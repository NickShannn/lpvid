#!/usr/bin/env python3
"""Validate explicitly supplied synthetic or generated v2 HTML reports.

Usage:
    python scripts/verify_report_v2.py /tmp/main_results_report_zh.html /tmp/main_results_report_en.html

At least one explicit HTML path is required. This script never scans runs directories.
"""

from __future__ import annotations

import html as html_lib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdc_pipeline.html_template_v2 import _tree_status_text

REPORT_DATA = re.compile(r"const REPORT_DATA = (\{.*?\});", re.S)
SECTION_IDS = ["overview", "risk", "sequences", "statistics", "abundance", "workflow", "links"]
ZH_TITLES = ["总体摘要", "风险提示", "代表序列", "详细序列统计", "相对丰度", "分析流程概览", "结果文件"]
EN_TITLES = [
    "Overview", "Prior-based risk notices", "Representative sequences", "Detailed sequence statistics",
    "Relative abundance", "Workflow overview", "Result files",
]
REQUIRED_IDS = [
    "chartTime", "chartSupergroup", "chartLength", "supergroupChartFrame", "lengthBinData",
    "stageTable", "sequenceTable", "abundanceTable", "globalFilter", "supergroupFilter",
]
RETIRED_TOKENS = ["chartFunnel", "chartCalls", "riskTable", "supergroupTable", "hostTable"]

FORBIDDEN_TERMS = [
    "管理层摘要", "结果解读", "Executive Summary",
    "Luca" + "Prot", "Luca" + "Virus", "Luca" + "Host",
    "embed" + "ding", "M" + "LP", "stage" + "-2", "top" + "-hit",
]
HAN = re.compile(r"[一-鿿]")


def _body(text: str) -> str:
    body = re.sub(r"<script.*?</script>", "", text, flags=re.S)
    return re.sub(r"<style.*?</style>", "", body, flags=re.S)


def _app_js(text: str) -> str:
    scripts = re.findall(r"<script[^>]*>(.*?)</script>", text, re.S)
    return scripts[-1] if scripts else ""


def _section(text: str, section_id: str) -> str:
    match = re.search(rf'<section id="{re.escape(section_id)}"[^>]*>(.*?)</section>', text, re.S)
    return match.group(1) if match else ""


def _div_container(text: str, marker: str) -> str:
    start = text.find(marker)
    if start < 0:
        return ""
    depth = 0
    for match in re.finditer(r"<div\b|</div>", text[start:]):
        if match.group(0) == "</div>":
            depth -= 1
            if depth == 0:
                return text[start:start + match.end()]
        else:
            depth += 1
    return ""


def _titles(text: str) -> list[str]:
    return [
        html_lib.unescape(re.sub(r"<[^>]+>", "", value)).strip()
        for value in re.findall(r"<h2[^>]*>(.*?)</h2>", text, re.S)
    ]


def _anchor_problems(fragment: str, label: str) -> list[str]:
    problems = []
    for anchor in re.findall(r"<a\b[^>]*>", fragment):
        if 'target="_blank"' not in anchor:
            problems.append(f"{label} 链接缺少 target=_blank: {anchor}")
        if 'rel="noopener"' not in anchor:
            problems.append(f"{label} 链接缺少 rel=noopener: {anchor}")
    return problems


def check_file(path: Path) -> list[str]:
    problems: list[str] = []
    text = path.read_text(encoding="utf-8")
    lang = "en" if re.search(r'<html[^>]*lang="en"', text) else "zh"
    body = _body(text)
    app_js = _app_js(text)
    match = REPORT_DATA.search(text)
    if not match:
        return [f"{path}: 未找到 REPORT_DATA"]
    try:
        data = json.loads(match.group(1))
    except json.JSONDecodeError as exc:
        return [f"{path}: REPORT_DATA 不是有效 JSON: {exc}"]

    if re.findall(r'<section id="([^"]+)"', text) != SECTION_IDS:
        problems.append("七个 section 的 DOM 顺序必须为 overview→risk→sequences→statistics→abundance→workflow→links")
    nav_match = re.search(r"<nav[^>]*>(.*?)</nav>", text, re.S)
    if not nav_match or re.findall(r'href="#([^"]+)"', nav_match.group(1)) != SECTION_IDS:
        problems.append("导航顺序与七段 DOM 顺序不一致")
    expected_titles = EN_TITLES if lang == "en" else ZH_TITLES
    if _titles(text) != expected_titles:
        problems.append(f"标题不匹配: {_titles(text)!r}")

    for element_id in REQUIRED_IDS:
        if f'id="{element_id}"' not in text:
            problems.append(f"缺少容器 #{element_id}")
    retired_region = body + "\n" + app_js
    for token in RETIRED_TOKENS:
        if token in retired_region:
            problems.append(f"存在旧 DOM/JS: {token}")
    if "echarts" in retired_region.lower():
        problems.append("存在 ECharts 残留")
    report_data_text = match.group(1)
    for term in FORBIDDEN_TERMS:
        lowered = term.lower()
        if any(lowered in value.lower() for value in (body, app_js, report_data_text)):
            problems.append(f"存在禁用词: {term}")
    for token in ("yAxes", "y2", "dualAxis"):
        if token in app_js:
            problems.append(f"图表疑似多轴配置: {token}")
    if "tooltip" not in app_js:
        problems.append("图表缺少 hover tooltip")
    if "#4f46e5" in text.lower() or "indigo" in body.lower():
        problems.append("存在旧营销紫色样式")

    stage_rows = data.get("stageRows") or []
    for row in stage_rows:
        for field in ("stage", "status", "input", "output", "minutes"):
            if field not in row:
                problems.append(f"stageRows 缺字段 {field}")
    if not stage_rows and not re.search(r"No stage metrics|没有阶段指标", _section(text, "workflow")):
        problems.append("空阶段指标缺少说明")
    if "chartFunnel" in text or "chartCalls" in text:
        problems.append("仍有重复 funnel 或结果构成图")

    counts = data.get("counts") or {}
    supergroup_rows = data.get("supergroupRows") or []
    if counts.get("representatives", 0) > 0 and not supergroup_rows:
        problems.append("存在代表序列但 supergroupRows 为空")
    for row in supergroup_rows:
        for field in ("representatives", "positive_orfs", "positive_contigs", "length_statistics", "classification_methods", "tree_statuses", "abundance"):
            if field not in row:
                problems.append(f"supergroupRows 缺字段 {field}: {row.get('supergroup')}")
        expected_status = _tree_status_text(row.get("tree_statuses") or {}, lang)
        if row.get("tree_status") != expected_status:
            problems.append(f"树状态不统一: {row.get('supergroup')}")
    height_match = re.search(r'id="supergroupChartFrame" style="height:(\d+)px"', text)
    expected_height = max(180, 72 + 44 * len(supergroup_rows))
    if not height_match or int(height_match.group(1)) != expected_height:
        problems.append(f"超群图高度未按类别数调整，应为 {expected_height}px")
    if not re.search(r'<details id="lengthBinData".*?<table', text, re.S):
        problems.append("长度图缺少可展开 bin 数据表")
    if "counted in more than one supergroup" not in body and "多个超群重复计入" not in body:
        problems.append("缺少跨超群 contig 计数说明")

    funnel = data.get("funnel") or []
    if [row.get("key") for row in funnel] != ["input", "orf", "prefilter", "positive", "representatives"]:
        problems.append("代表序列五节点流程不完整或顺序错误")
    for row in funnel:
        if not row.get("unit"):
            problems.append(f"流程步骤 {row.get('key')} 缺少单位（包括 0 数量）")
        node = re.search(rf'id="flow-{re.escape(str(row.get("key")))}">([^<]+)</strong>', text)
        if not node or node.group(1).replace(",", "") != str(row.get("count", 0)):
            problems.append(f"流程节点数量与数据不一致: {row.get('key')}")
    if len(re.findall(r'class="flow-node"', body)) != 5:
        problems.append("代表序列流程必须显示五个节点")

    sequence_columns = data.get("sequenceColumns") or []
    reference = next((row for row in sequence_columns if row.get("field") == "reference_id"), None)
    expected_reference = "Original reference-sequence ID" if lang == "en" else "原始参考序列 ID"
    if not reference or reference.get("title") != expected_reference:
        problems.append("reference_id 未明确标为原始参考序列 ID")
    for retired_field in ("host_human", "host_mammal", "host_vertebrate", "host_suppression_reason"):
        if any(column.get("field") == retired_field for column in sequence_columns):
            problems.append(f"序列表仍含旧列: {retired_field}")
    if "Object.values(row).some" not in app_js:
        problems.append("global search 不是任意字段 OR 匹配")
    all_groups = "All supergroups" if lang == "en" else "全部超群"
    if f'<option value="">{all_groups}</option>' not in text:
        problems.append("超群筛选缺少全部选项")
    for token, threshold in (("makeSequenceOptions", "rowCount > 20"), ("makeContigOptions", "rowCount > 10")):
        if token not in app_js or threshold not in app_js:
            problems.append(f"表格未按数据量动态设置高度/分页: {token}")
    if re.search(r"height\s*:\s*['\"]300px", app_js):
        problems.append("小数据表格仍使用固定 300px 高度")

    risk_summary = data.get("riskSummary") or {}
    for field in ("human_nn1", "human_nn2_5", "vertebrate_only", "no_placed_label"):
        if field not in risk_summary:
            problems.append(f"riskSummary 缺字段 {field}")
    risk_ids = {
        "human-nn1": "human_nn1", "human-nn2-5": "human_nn2_5",
        "vertebrate-only": "vertebrate_only", "no-placed-label": "no_placed_label",
    }
    for element_id, field in risk_ids.items():
        value_match = re.search(rf'id="risk-{element_id}">(\d+)</strong>', text)
        if not value_match or int(value_match.group(1)) != int(risk_summary.get(field, 0) or 0):
            problems.append(f"风险 KPI 未直接使用 risk_summary: {field}")
    risk_rows = data.get("riskRows") or []
    if risk_rows:
        if 'id="risk-banner"' not in text:
            problems.append("风险总横幅缺失")
        matches = re.findall(r'class="risk-family-card"[^>]*data-matched="(true|false)"', text)
        if "true" in matches and "false" in matches and matches.index("true") > matches.index("false"):
            problems.append("匹配名录的科卡片未优先")
        if not re.search(r"does not exclude unlisted|不排除名录外", body):
            problems.append("风险提示缺少非排除性说明")
    elif not re.search(r"No sequences were placed|没有序列完成系统发育放置|no usable nn1|没有可用的 nn1", _section(text, "risk")):
        problems.append("无放置/无 nn1 风险空状态不明确")
    for row in risk_rows:
        for neighbor in row.get("neighbors") or []:
            if neighbor.get("distance") is None:
                expected = "unavailable" if lang == "en" else "不可用"
                if expected not in body:
                    problems.append("null 距离未显示 unavailable")
                    break
    if risk_summary.get("human_nn2_5", 0):
        if not any(
            any(isinstance(n.get("rank"), int) and 2 <= n["rank"] <= 5 and n.get("human") is True for n in row.get("neighbors") or [])
            for row in risk_rows
        ):
            problems.append("human_nn2_5 没有真实 rank 2-5 证据")

    abundance_rows = data.get("abundanceRows") or []
    if data.get("route") == "reads" and counts.get("positive_contigs", 0) and not abundance_rows:
        problems.append("reads 路线有阳性 contigs 但 abundanceRows 为空")
    if data.get("route") != "reads" and not re.search(r"not applicable|不适用", _section(text, "abundance")):
        problems.append("非 reads 路线缺少丰度不适用说明")

    links = _section(text, "links")
    problems.extend(_anchor_problems(links, "结果文件"))
    for name, description in (("risk_warning.tsv", "Prior-based family notices" if lang == "en" else "先验风险提示"), ("neighbor_evidence.tsv", "Structured top-5 neighbor evidence" if lang == "en" else "结构化 top-5 近邻证据")):
        if name in links and description not in links:
            problems.append(f"新结果文件缺少描述: {name}")
    if counts.get("placed", 0):
        iframe = re.search(r'<iframe id="treePreview"[^>]*>', text)
        if not iframe or 'loading="lazy"' not in iframe.group(0):
            problems.append("placed>0 时缺少 lazy iframe 树预览")
        elif not re.search(r'src="(?:phylogenetic_explorer|tree_explorer)\.html"', iframe.group(0)):
            problems.append("树预览 iframe 不是同目录 explorer")
        standalone = re.search(r'<a id="treeStandalone"[^>]*>', text)
        if not standalone:
            problems.append("缺少独立树浏览器链接")
        else:
            problems.extend(_anchor_problems(standalone.group(0), "树浏览器"))
        preview_container = _div_container(text, '<div class="tree-preview no-print">')
        if not preview_container:
            problems.append("缺少 tree-preview 容器")
        elif 'id="treePreview"' not in preview_container:
            problems.append("树预览 iframe 未位于 tree-preview 容器内")
        if 'class="tree-preview no-print"' not in text or not re.search(r"@media print.*?\.tree-preview[^{}]*\{[^}]*display\s*:\s*none", text, re.S):
            problems.append("树预览未设置为 print 隐藏")
    elif 'id="treePreview"' in text or 'id="treeStandalone"' in text:
        problems.append("无放置结果不应显示树预览")
    if "per-query" in body.lower():
        problems.append("不应声称 per-query 树联动")

    if counts.get("representatives", 0) == 0 and not re.search(r"No representative sequences|没有可用的代表序列", _section(text, "sequences")):
        problems.append("空代表序列缺少说明")
    if not supergroup_rows and not re.search(r"No supergroup statistics|没有可用的超群统计", _section(text, "statistics")):
        problems.append("空超群统计缺少说明")

    if lang == "en":
        han = HAN.search(body)
        if han:
            snippet = body[max(0, han.start() - 30):han.start() + 30].replace("\n", " ")
            problems.append(f"英文报告正文含中文: …{snippet}…")
        for row in funnel:
            if HAN.search(str(row.get("unit", ""))) or HAN.search(str(row.get("stage", ""))):
                problems.append(f"英文流程未翻译: {row.get('key')}")
        for row in stage_rows:
            if any(HAN.search(str(row.get(field, ""))) for field in ("stage", "input", "output")):
                problems.append(f"英文 stageRows 未翻译: {row.get('stage')}")
        for row in supergroup_rows:
            if HAN.search(str(row.get("tree_status", ""))):
                problems.append(f"英文 tree_status 未翻译: {row.get('supergroup')}")
        for row in risk_rows:
            for field in ("pathogen_category", "statutory_category", "pathogen_examples", "statutory_examples", "hint"):
                if HAN.search(str(row.get(field, ""))):
                    problems.append(f"英文风险字段复用了中文: {field}")
        for column in sequence_columns:
            if HAN.search(str(column.get("title", ""))) or HAN.search(str(column.get("headerTooltip", ""))):
                problems.append(f"英文序列列未翻译: {column.get('field')}")
        for key, value in (data.get("labels") or {}).items():
            if HAN.search(str(value)):
                problems.append(f"英文前端标签未翻译: {key}")
        if HAN.search(json.dumps(risk_summary, ensure_ascii=False)):
            problems.append("英文 riskSummary 含中文类别")
    if "条条" in body or "Non-Probability" in body:
        problems.append("存在重复量词或退休指标")
    return problems


def main(argv: list[str]) -> int:
    paths = [Path(item) for item in argv[1:]]
    if not paths:
        print("请显式提供至少一个 HTML 路径 / explicit HTML paths are required; no directory scan is performed.", file=sys.stderr)
        return 2
    failures = 0
    for path in paths:
        if not path.is_file():
            problems = [f"文件不存在: {path}"]
        else:
            problems = check_file(path)
        if problems:
            failures += 1
            print(f"FAIL {path}")
            for problem in problems:
                print(f"     - {problem}")
        else:
            print(f"OK   {path}")
    print(f"\n{len(paths) - failures}/{len(paths)} 通过")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
