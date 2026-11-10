"""Offline, bilingual seven-section HTML report."""

from __future__ import annotations

import html
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.html_report import length_statistics

NAVY = "#17324d"
TEAL = "#0f766e"
ASSET_DIR = PROJECT_ROOT / "third_party" / "report_assets"
SECTION_IDS = ("overview", "risk", "sequences", "statistics", "abundance", "workflow", "links")


def _e(value) -> str:
    return html.escape(str(value), quote=True)


def _json_script(value) -> str:
    return (
        json.dumps(value, ensure_ascii=False, allow_nan=False)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("</script", "<\\/script")
    )


def _asset(name: str) -> str:
    path = ASSET_DIR / name
    if not path.is_file():
        raise FileNotFoundError(f"报告前端资源缺失: {path}；请运行 scripts/sync_report_assets.py")
    return path.read_text(encoding="utf-8")


def _unavailable(lang: str) -> str:
    return "unavailable" if lang == "en" else "不可用"


def _display(value, lang: str) -> str:
    return _unavailable(lang) if value is None or value == "" else str(value)


def _flag(value, lang: str) -> str:
    if value is True:
        return "yes" if lang == "en" else "是"
    if value is False:
        return "no" if lang == "en" else "否"
    return _unavailable(lang)


_AVAILABLE_TREE_STATUSES = {"assigned", "classified_available", "available"}
_CONSTRUCTION_TREE_STATUSES = {
    "classified_no_tree", "classified_pending_merge", "classified_pending_rename",
    "no_tree", "pending_merge", "pending_rename",
}
_UNASSIGNED_TREE_STATUSES = {"unassigned"}


def _unassigned_label(lang: str) -> str:
    return "Unclassified" if lang == "en" else "未分类"


def _supergroup_label(value, lang: str) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"unassigned", "unclassified"} or text == "未分类":
        return _unassigned_label(lang)
    return text


def _is_unassigned_supergroup(value) -> bool:
    text = str(value or "").strip().lower()
    return not text or text in {"unassigned", "unclassified", "未分类"}


def _tree_status_text(statuses: dict, lang: str) -> str:
    zh = lang != "en"
    active = []
    for key, count in (statuses or {}).items():
        text = str(key).strip()
        if not text:
            continue
        try:
            count_value = int(count)
        except (TypeError, ValueError):
            count_value = 1 if count else 0
        if count_value > 0:
            active.append(text)
    if not active:
        return "不可用" if zh else "Unavailable"
    status_set = set(active)
    unknown = status_set - (
        _AVAILABLE_TREE_STATUSES | _CONSTRUCTION_TREE_STATUSES | _UNASSIGNED_TREE_STATUSES
    )
    if unknown:
        return "不可用" if zh else "Unavailable"
    if status_set <= _UNASSIGNED_TREE_STATUSES:
        return "未分类，不进入系统发育放置" if zh else "Unclassified; not placed"
    available = status_set & _AVAILABLE_TREE_STATUSES
    construction = status_set & _CONSTRUCTION_TREE_STATUSES
    if available and construction:
        return "部分有可用骨架树，其余正在建设中" if zh else "Mixed availability"
    if available:
        return "有可用骨架树" if zh else "Backbone tree available"
    if construction:
        return "骨架树正在建设中" if zh else "Backbone tree under construction"
    return "不可用" if zh else "Unavailable"


def _classification_status_text(value, lang: str) -> str:
    text = str(value or "").strip()
    if not text:
        return "不可用" if lang != "en" else "Unavailable"
    return _tree_status_text({text: 1}, lang)


def _method_label(value, lang: str) -> str:
    zh = lang != "en"
    text = str(value or "").strip()
    if not text:
        return _unavailable(lang)
    if text.lower() in {"reference", "reference match", "reference_match", "dia" + "mond"}:
        return "参考序列匹配" if zh else "Reference match"
    if text.lower() in {"classifier", "generic classification", "emb" + "edding"}:
        return "分类模型" if zh else "Classification model"
    if text.lower() in {"unassigned", "unclassified"} or text == "未分类":
        return _unassigned_label(lang)
    return text


def _stage_status_text(value, lang: str) -> str:
    text = str(value or "").strip()
    if lang == "en":
        return {
            "done": "Completed",
            "skipped": "Skipped",
            "not_applicable": "Not applicable",
            "not_run": "Not run",
            "failed": "Failed",
            "running": "Running",
            "stale": "Stale",
        }.get(text, text or _unavailable(lang))
    return {
        "done": "完成",
        "skipped": "跳过",
        "not_applicable": "不适用",
        "not_run": "未运行",
        "failed": "失败",
        "running": "运行中",
        "stale": "过时",
    }.get(text, text or _unavailable(lang))


def _coverage_text(value, lang: str) -> str:
    if value is None or value == "":
        return _unavailable(lang)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return _unavailable(lang)
    if not math.isfinite(number):
        return _unavailable(lang)
    return f"{_number_text(value, lang)}%"


def _abundance_text(abundance: dict, route: str, lang: str) -> str:
    zh = lang != "en"
    if route != "reads":
        return "路线不适用（非 reads 输入）" if zh else "Not applicable for this input route"
    try:
        n_contigs = int(abundance.get("n_contigs", 0) or 0)
    except (TypeError, ValueError):
        n_contigs = 0
    if n_contigs <= 0:
        return "不可用（reads 路线无可用数据）" if zh else "Unavailable (no usable data on the reads route)"
    return (
        f"RPM {_number_text(abundance.get('RPM'), lang)} · "
        f"{'回比 reads' if zh else 'mapped reads'} {_number_text(abundance.get('mapped_reads'), lang)} · "
        f"{'加权深度' if zh else 'weighted depth'} {_number_text(abundance.get('mean_depth_weighted'), lang)} · "
        f"{'覆盖度' if zh else 'coverage'} {_coverage_text(abundance.get('coverage_breadth_pct'), lang)}"
    )


def _percent_number(value):
    """Convert a public raw 0-1 probability to a percentage exactly once."""
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return ""
    return round(number * 100, 2)


def _number_text(value, lang: str, digits: int = 2) -> str:
    if value is None or value == "":
        return _unavailable(lang)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return _unavailable(lang)
    if number.is_integer():
        return f"{int(number):,}"
    return f"{number:,.{digits}f}".rstrip("0").rstrip(".")


def _unit(value: str, lang: str, count: int | None = None) -> str:
    plural = {
        "条核酸 contigs": ("nucleotide contigs", "nucleotide contig"),
        "条蛋白 ORFs": ("protein ORFs", "protein ORF"),
        "条代表蛋白": ("representative proteins", "representative protein"),
        "条蛋白序列": ("protein sequences", "protein sequence"),
        "条 reads": ("reads", "read"),
        "份报告": ("reports", "report"),
    }
    if lang == "en":
        if value in plural:
            return plural[value][1] if count == 1 else plural[value][0]
        return value
    return value


def _stage_label(row: dict, lang: str) -> str:
    if lang != "en":
        return row.get("stage", "")
    return {
        "orf": "Protein ORF prediction",
        "prefilter": "Host/cellular homology prefilter",
        "detect": "RdRP recognition",
        "representatives": "Positive sequence deduplication",
        "supergroup": "Supergroup classification",
        "classifier": "Classification",
        "phylo": "Phylogenetic placement (LucaPhylo)",
        "abundance": "Read-mapping abundance",
        "report": "Results summary and visualization",
        "preprocess": "Read QC and assembly",
    }.get(row.get("stage_key", ""), row.get("stage", ""))


def _funnel_label(row: dict, lang: str) -> str:
    if lang != "en":
        return row.get("stage", "")
    return {
        "input": "Nucleotide contigs",
        "orf": "Protein ORFs",
        "length": "Length filter (≥150 aa)",
        "prefilter": "Prefilter-passing proteins",
        "positive": "RdRP-positive proteins",
        "representatives": "Representative proteins",
        "supergroup": "Supergroup classification",
        "phylo": "Phylogenetic placement",
    }.get(row.get("key", ""), row.get("stage", ""))


def _call_label(value: str, lang: str) -> str:
    if value in ("", "negative", "both_neg"):
        return "未检出" if lang != "en" else "Not detected"
    return "检出" if lang != "en" else "Detected"


def _sequence_data(rows: list[dict], lang: str) -> tuple[list[dict], list[dict]]:
    data = []
    for row in rows:
        data.append({
            "query_id": row.get("query_id", ""),
            "original_orf_id": row.get("original_orf_id", ""),
            "reference_id": row.get("reference_id", ""),
            "contig_id": row.get("contig_id", ""),
            "aa_len": row.get("aa_len", ""),
            "call": _call_label(row.get("call", ""), lang),
            "rdrp_prob": _percent_number(row.get("rdrp_prob")),
            "supergroup": _supergroup_label(row.get("supergroup"), lang),
            "classification_method": _method_label(row.get("classification_method"), lang),
            "classification_status": _classification_status_text(row.get("classification_status"), lang),
            "pident": row.get("pident", ""),
            "evalue": row.get("evalue", ""),
            "nn1": row.get("nn1", ""),
            "nn1_dist": row.get("nn1_dist") if row.get("nn1_dist") is not None else "",
            "rpm": row.get("rpm", ""),
            "mean_depth": row.get("mean_depth", ""),
            "coverage_breadth_pct": row.get("coverage_breadth_pct", ""),
        })
    zh = lang != "en"
    titles = {
        "query_id": "内部 ID" if zh else "Query ID",
        "original_orf_id": "原始 ORF ID" if zh else "Original ORF ID",
        "reference_id": "原始参考序列 ID" if zh else "Original reference-sequence ID",
        "contig_id": "源核酸 contig" if zh else "Source contig",
        "aa_len": "蛋白长度 (aa)" if zh else "Protein length (aa)",
        "call": "检出结论" if zh else "Detection call",
        "rdrp_prob": "RdRP 识别概率 (%)" if zh else "RdRP probability (%)",
        "supergroup": "超群" if zh else "Supergroup",
        "classification_method": "分类方法" if zh else "Method",
        "classification_status": "分类状态" if zh else "Status",
        "pident": "参考序列一致性 (%)" if zh else "Reference sequence identity (%)",
        "evalue": "比对 E-value" if zh else "Alignment E-value",
        "nn1": "最近骨架序列" if zh else "Nearest backbone sequence",
        "nn1_dist": "最近邻距离" if zh else "Nearest-neighbor distance",
        "rpm": "RPM",
        "mean_depth": "平均深度" if zh else "Mean depth",
        "coverage_breadth_pct": "覆盖度 (%)" if zh else "Coverage (%)",
    }
    tooltips = {
        "query_id": "本次运行内部序列 ID" if zh else "Internal sequence ID for this run",
        "original_orf_id": "fastorf 原始 ORF ID，可在原始 FASTA 中检索" if zh else "Original fastorf ORF ID for locating in FASTA",
        "reference_id": "原始历史参考序列 ID；不等于当前超群归属或分类名称" if zh else "Original historical reference-sequence ID; it is not the current supergroup assignment or classification name",
        "contig_id": "源核酸 contig ID" if zh else "Source nucleotide contig ID",
        "aa_len": "蛋白序列长度（氨基酸）" if zh else "Protein sequence length (aa)",
        "call": "RdRP 识别模型判定" if zh else "RdRP recognition call",
        "rdrp_prob": "RdRP 识别模型阳性概率（%）" if zh else "RdRP recognition probability (%)",
        "supergroup": "代表序列所属超群" if zh else "Supergroup assigned to the representative",
        "classification_method": "超群分类方法" if zh else "Supergroup classification method",
        "classification_status": "分类与骨架树状态" if zh else "Classification and tree availability status",
        "pident": "参考序列比对一致性（%）" if zh else "Reference sequence alignment identity (%)",
        "evalue": "参考序列比对 E-value" if zh else "Reference sequence alignment E-value",
        "nn1": "距离最近的骨架树序列" if zh else "Nearest backbone sequence",
        "nn1_dist": "最近邻距离；缺失时显示不可用" if zh else "Nearest-neighbor distance; unavailable when missing",
        "rpm": "质控后 reads 回比的相对丰度（RPM）" if zh else "Relative abundance from post-QC read mapping (RPM)",
        "mean_depth": "回比平均深度" if zh else "Mean mapping depth",
        "coverage_breadth_pct": "回比覆盖广度（%）" if zh else "Mapping coverage breadth (%)",
    }
    fields = list(titles)
    numeric = {"aa_len", "rdrp_prob", "pident", "evalue", "nn1_dist", "rpm", "mean_depth", "coverage_breadth_pct"}
    filterable = {"query_id", "original_orf_id", "reference_id", "contig_id", "supergroup", "call", "classification_method", "classification_status"}
    columns = [
        {
            "field": field,
            "title": titles[field],
            "headerTooltip": tooltips[field],
            "sorter": "number" if field in numeric else "string",
            "hozAlign": "right" if field in numeric else "left",
            "headerFilter": "input" if field in filterable else False,
            "minWidth": 110 if field not in {"original_orf_id", "reference_id", "nn1"} else 175,
        }
        for field in fields
    ]
    return data, columns


def _stage_rows(summary: dict, lang: str) -> list[dict]:
    return [
        {
            "stage": _stage_label(row, lang),
            "status": _stage_status_text(row.get("status", ""), lang),
            "input": f"{row.get('input_count', 0):,} {_unit(row.get('input_unit', ''), lang, row.get('input_count'))}",
            "output": f"{row.get('output_count', 0):,} {_unit(row.get('output_unit', ''), lang, row.get('output_count'))}",
            "minutes": row.get("duration_minutes", 0),
        }
        for row in (summary.get("stage_metrics") or [])
    ]


def _object_dict(value) -> dict:
    if isinstance(value, dict):
        return dict(value)
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    if hasattr(value, "__dict__"):
        return dict(vars(value))
    return {}


def _risk_dicts(rows) -> list[dict]:
    output = []
    for raw in rows or []:
        row = _object_dict(raw)
        neighbors = []
        for raw_neighbor in row.get("neighbors") or []:
            neighbor = _object_dict(raw_neighbor)
            neighbors.append({
                "rank": neighbor.get("rank"),
                "reference_id": neighbor.get("reference_id", ""),
                "distance": neighbor.get("distance"),
                "family": neighbor.get("family", ""),
                "genus": neighbor.get("genus", ""),
                "human": neighbor.get("human"),
                "vertebrate": neighbor.get("vertebrate"),
                "plant": neighbor.get("plant"),
            })
        row["neighbors"] = neighbors
        output.append(row)
    return output


def _catalog_list_text(row: dict, key: str, lang: str) -> str:
    if row.get(key):
        return _category(row.get(key), lang)
    if row.get("family"):
        return "No match in supplied lists" if lang == "en" else "未匹配名录"
    return "taxonomy unavailable" if lang == "en" else "分类学信息不可用"


def _category(value: str, lang: str) -> str:
    if lang != "en" or not value:
        return value or _unavailable(lang)
    return {
        "第一类": "Category I", "第二类": "Category II", "第三类": "Category III", "第四类": "Category IV",
        "甲类": "Class A", "乙类": "Class B", "丙类": "Class C",
    }.get(value, value)


def _risk_limitation(row: dict, lang: str) -> str:
    """面向公众的边界说明：由结构化字段生成，不回读上游内部提示文案。"""
    zh = lang != "en"
    family = row.get("family") or ""
    matched = bool(row.get("pathogen_category") or row.get("statutory_category"))
    if zh:
        followup = "这不排除名录外或未识别病毒，需进一步核查，不构成新序列风险分级或病原确证。"
        if not family:
            return f"第 1 近邻所属科未识别，无法匹配名录；{followup}"
        if not matched:
            return f"第 1 近邻所属科（{family}）未匹配到名录收录成员；{followup}"
        dimensions = _risk_dimension_text(row, lang)
        return (
            f"第 1 近邻所属科（{family}）的已知成员中有名录收录成员：{dimensions}；"
            "这仅是类群归属的先验提示，需进一步核查，不构成新序列风险分级或病原确证。"
        )
    if not family:
        return (
            "The nearest-neighbor family is unavailable, so no catalog match can be made. This does not exclude "
            "unlisted or unidentified viruses. Further review is required; this is not a risk "
            "classification or pathogen confirmation for the query."
        )
    if not matched:
        return (
            f"No entry in either list was matched for the nearest-neighbor family ({family}). This does not "
            "exclude unlisted or unidentified viruses. Further review is required; this is not a risk "
            "classification or pathogen confirmation for the query."
        )
    dimensions = _risk_dimension_text(row, lang)
    return (
        f"Known members of the nearest-neighbor family ({family}) match the indicated catalog/list entries "
        f"({dimensions}). This prior-based notice is not a risk classification or pathogen confirmation "
        "for the query; further review is required."
    )


def _risk_dimension_text(row: dict, lang: str) -> str:
    """名录维度摘要；英文页面不携带中文示例名。"""
    zh = lang != "en"
    pathogen_label = "病原微生物目录" if zh else "pathogen catalog"
    statutory_label = "法定传染病名单" if zh else "statutory disease list"
    parts = []
    for category_key, examples_key, label in (
        ("pathogen_category", "pathogen_examples", pathogen_label),
        ("statutory_category", "statutory_examples", statutory_label),
    ):
        if not row.get(category_key):
            continue
        detail = _category(row[category_key], lang)
        examples = row.get(examples_key) or ""
        if zh and examples:
            detail += f"（{examples}）"
        parts.append(f"{label}{detail}" if zh else f"{label} {detail}")
    return "、".join(parts) if zh else "; ".join(parts)


_PATHOGEN_SEVERITY = {"第一类": "high", "第二类": "mid", "第三类": "neutral", "第四类": "neutral"}
_STATUTORY_SEVERITY = {"甲类": "high", "乙类": "mid", "丙类": "neutral"}
_PATHOGEN_ORDER = ["第一类", "第二类", "第三类", "第四类"]
_STATUTORY_ORDER = ["甲类", "乙类", "丙类"]


def _severity_class(kind: str, value: str) -> str:
    table = _PATHOGEN_SEVERITY if kind == "pathogen" else _STATUTORY_SEVERITY
    return f"badge-sev-{table.get(value, 'none')}"


def _source_badge(source_label: str, kind: str, value: str, lang: str, *, with_label: bool = True) -> str:
    label_html = f"<strong>{_e(source_label)}</strong>" if with_label else ""
    return (
        f'<span class="badge {_severity_class(kind, value)}">'
        f"{label_html}{_e(_category(value, lang))}</span>"
    )


def _risk_badge_html(source_label: str, kind: str, values: list[str], lang: str, *, with_label: bool = True) -> str:
    ordered = [value for value in (_PATHOGEN_ORDER if kind == "pathogen" else _STATUTORY_ORDER) if value in values]
    ordered += [value for value in values if value not in ordered]
    return "".join(_source_badge(source_label, kind, value, lang, with_label=with_label) for value in ordered)


def _risk_neighbor_item(neighbor: dict, lang: str) -> str:
    zh = lang != "en"
    rank = neighbor.get("rank")
    if zh:
        rank_label = f"第 {rank} 近邻" if rank is not None else "排序位次不可用的近邻"
        fields = [
            f"参考序列 {_display(neighbor.get('reference_id'), lang)}",
            f"距离 {_display(neighbor.get('distance'), lang)}",
            f"人源标签 {_flag(neighbor.get('human'), lang)}",
            f"脊椎动物标签 {_flag(neighbor.get('vertebrate'), lang)}",
            f"科 {_display(neighbor.get('family'), lang)}",
        ]
    else:
        rank_label = f"Neighbor {rank}" if rank is not None else "Neighbor rank unavailable"
        fields = [
            f"Reference {_display(neighbor.get('reference_id'), lang)}",
            f"Distance {_display(neighbor.get('distance'), lang)}",
            f"Human label {_flag(neighbor.get('human'), lang)}",
            f"Vertebrate label {_flag(neighbor.get('vertebrate'), lang)}",
            f"Family {_display(neighbor.get('family'), lang)}",
        ]
    body = "".join(f"<span>{_e(field)}</span>" for field in fields)
    return f'<li class="neighbor-item"><strong class="neighbor-rank">{_e(rank_label)}</strong>{body}</li>'


def _risk_html(risk_rows: list[dict], risk_summary: dict, counts: dict, lang: str) -> tuple[str, list[dict]]:
    zh = lang != "en"
    kpis = [
        ("human-nn1", "human_nn1", "人源第 1 近邻" if zh else "Human nearest neighbor", "第 1 近邻带人源标签" if zh else "A human label on the nearest neighbor"),
        ("human-nn2-5", "human_nn2_5", "第 2 至第 5 近邻含人源" if zh else "Human label at neighbors 2-5", "第 2 至第 5 近邻中含人源标签，且第 1 近邻非人源" if zh else "Human labels among neighbors 2-5 while the nearest neighbor is not human"),
        ("vertebrate-only", "vertebrate_only", "仅脊椎动物近邻" if zh else "Vertebrate-only labels", "前 5 个近邻有脊椎动物标签且无人源标签" if zh else "Vertebrate labels among the five nearest neighbors without human labels"),
        ("no-placed-label", "no_placed_label", "无已放置人源/脊椎动物标签" if zh else "No placed human/vertebrate label", "不表示宿主缺失或无风险" if zh else "Does not mean host absence or no risk"),
    ]
    kpi_html = "".join(
        f'<div class="kpi risk-kpi"><span class="kpi-label">{_e(label)}</span>'
        f'<strong id="risk-{key}">{int(risk_summary.get(summary_key, 0) or 0)}</strong>'
        f'<span class="kpi-note">{_e(note)}</span></div>'
        for key, summary_key, label, note in kpis
    )
    intro = (
        "以下内容仅为类群归属与名录收录情况的先验提示，不构成新序列的风险分级或病原确证；识别结果与系统发育放置均需进一步核查。"
        if zh else
        "These are taxonomy and catalog prior notices, not risk classifications or pathogen confirmations for new sequences. Detection results and phylogenetic placements require further review."
    )
    if not risk_rows:
        if counts.get("placed", 0):
            empty = "已放置代表序列没有可用的第 1 近邻证据，无法生成先验提示；这不排除名录外或未识别病毒。" if zh else "Placed sequences have no usable nearest-neighbor evidence, so prior notices are unavailable; this does not exclude unlisted or unidentified viruses."
        else:
            empty = "没有序列完成系统发育放置，因此无先验科提示；这不排除名录外或未识别病毒。" if zh else "No sequences were placed, so prior-based family notices are unavailable; this does not exclude unlisted or unidentified viruses."
        return f'<p class="section-intro">{_e(intro)}</p><div class="kpi-grid risk-grid">{kpi_html}</div><div id="risk-empty" class="empty-state">{_e(empty)}</div>', []

    identified_rows = [row for row in risk_rows if row.get("family")]
    unknown_rows = [row for row in risk_rows if not row.get("family")]
    matched_rows = [
        row for row in identified_rows
        if row.get("pathogen_category") or row.get("statutory_category")
    ]
    matched_families = {row.get("family") for row in matched_rows}
    identified_families = {row.get("family") for row in identified_rows}
    matched_query_count = len(matched_rows)
    unmatched_identified_query_count = len(identified_rows) - len(matched_rows)
    unknown_query_count = len(unknown_rows)
    total_query_count = len(risk_rows)
    if zh:
        banner = (
            f"{matched_query_count} 条代表序列 / {len(matched_families)} 个科匹配至少一个名录；"
            f"已放置代表序列 {total_query_count} 条 = {matched_query_count} 匹配 + "
            f"{unmatched_identified_query_count} 未匹配 + {unknown_query_count} 科级分类未确定；"
            f"已识别科 {len(identified_families)} 个。"
        )
    else:
        query_word = "query" if matched_query_count == 1 else "queries"
        family_word = "family" if len(matched_families) == 1 else "families"
        placed_word = "query" if total_query_count == 1 else "queries"
        identified_word = "family" if len(identified_families) == 1 else "families"
        banner = (
            f"{matched_query_count} {query_word} / {len(matched_families)} {family_word} matched at least one list; "
            f"{total_query_count} placed {placed_word} = {matched_query_count} matched + "
            f"{unmatched_identified_query_count} unmatched + {unknown_query_count} with undetermined family classification; "
            f"{len(identified_families)} identified non-empty {identified_word}."
        )
    grouped = defaultdict(list)
    for row in risk_rows:
        grouped[row.get("family") or ""].append(row)
    ordered = sorted(
        grouped.items(),
        key=lambda item: (
            not any(row.get("pathogen_category") or row.get("statutory_category") for row in item[1]),
            item[0] or "~",
        ),
    )
    unmatched_groups = [
        (family, members)
        for family, members in ordered
        if family and not any(row.get("pathogen_category") or row.get("statutory_category") for row in members)
    ]
    unmatched_family_count = len(unmatched_groups)
    unmatched_query_count = sum(len(members) for _, members in unmatched_groups)
    matched_cards: list[str] = []
    unmatched_cards: list[str] = []
    serialized = []
    for family, members in ordered:
        members = sorted(members, key=lambda row: str(row.get("query_id", "")))
        matched = bool(family) and any(row.get("pathogen_category") or row.get("statutory_category") for row in members)
        pathogen_values = sorted(
            {row.get("pathogen_category") for row in members if row.get("pathogen_category")},
            key=lambda value: _PATHOGEN_ORDER.index(value) if value in _PATHOGEN_ORDER else len(_PATHOGEN_ORDER),
        )
        statutory_values = sorted(
            {row.get("statutory_category") for row in members if row.get("statutory_category")},
            key=lambda value: _STATUTORY_ORDER.index(value) if value in _STATUTORY_ORDER else len(_STATUTORY_ORDER),
        )
        supergroups = sorted({str(row.get("supergroup") or "") for row in members} - {""})
        pathogen_label = "病原微生物目录" if zh else "Pathogen catalog"
        statutory_label = "法定传染病名单" if zh else "Statutory disease list"
        if matched:
            summary_badges = _risk_badge_html(pathogen_label, "pathogen", pathogen_values, lang) + _risk_badge_html(statutory_label, "statutory", statutory_values, lang)
        elif family:
            summary_badges = f'<span class="badge badge-sev-none">{_e("未匹配名录" if zh else "No catalog/list match")}</span>'
        else:
            summary_badges = f'<span class="badge badge-sev-none">{_e("科级分类未确定" if zh else "Undetermined family")}</span>'
        family_label = family or ("第 1 近邻所属科不可用" if zh else "Nearest-neighbor family unavailable")
        if zh:
            count_text = f"{len(members)} 条代表序列"
        else:
            count_text = f"{len(members)} representative sequence" + ("" if len(members) == 1 else "s")
        if not family:
            pathogen_cell = statutory_cell = ("分类学信息不可用" if zh else "taxonomy unavailable")
        elif pathogen_values:
            pathogen_cell = _risk_badge_html(pathogen_label, "pathogen", pathogen_values, lang, with_label=False)
        else:
            pathogen_cell = "未匹配名录" if zh else "No match in supplied lists"
        if family and statutory_values:
            statutory_cell = _risk_badge_html(statutory_label, "statutory", statutory_values, lang, with_label=False)
        elif family:
            statutory_cell = "未匹配名录" if zh else "No match in supplied lists"
        overview = (
            f'<div class="risk-family-overview"><dl class="overview-grid">'
            f'<dt>{_e("代表序列数" if zh else "Representatives")}</dt><dd>{_e(count_text)}</dd>'
            f'<dt><strong class="catalog-label">{_e(pathogen_label)}</strong></dt><dd>{pathogen_cell}</dd>'
            f'<dt><strong class="catalog-label">{_e(statutory_label)}</strong></dt><dd>{statutory_cell}</dd>'
            f'<dt>{_e("涉及超群" if zh else "Supergroups")}</dt><dd>{_e("、".join(supergroups) if zh else ", ".join(supergroups)) or _e(_unavailable(lang))}</dd>'
            f"</dl></div>"
        )
        query_html = []
        for row in members:
            pathogen = _catalog_list_text(row, "pathogen_category", lang)
            statutory = _catalog_list_text(row, "statutory_category", lang)
            pathogen_examples = row.get("pathogen_examples", "") if zh else ""
            statutory_examples = row.get("statutory_examples", "") if zh else ""
            limitation = _risk_limitation(row, lang)
            primary_distance = row.get("nn1_dist")
            if primary_distance is None:
                primary = next((n for n in row["neighbors"] if n.get("rank") == 1), None)
                primary_distance = primary.get("distance") if primary else None
            labels = {
                "original": "原始 ORF ID" if zh else "Original ORF ID",
                "supergroup": "超群" if zh else "Supergroup",
                "primary": "第 1 近邻" if zh else "Nearest neighbor",
                "pathogen": "病原微生物目录" if zh else "Pathogen catalog",
                "statutory": "法定传染病名单" if zh else "Statutory disease list",
            }
            pathogen_text = pathogen + (f"（{pathogen_examples}）" if pathogen_examples and row.get("pathogen_category") else "")
            statutory_text = statutory + (f"（{statutory_examples}）" if statutory_examples and row.get("statutory_category") else "")
            neighbors = "".join(_risk_neighbor_item(n, lang) for n in row["neighbors"])
            if not neighbors:
                neighbors = f'<li class="neighbor-item">{_e(_unavailable(lang))}</li>'
            query_html.append(
                f'<div class="risk-query"><h4>{_e("代表序列" if zh else "Query")} {_e(row.get("query_id", ""))}</h4>'
                f'<dl class="evidence-grid">'
                f'<dt>{_e(labels["original"])}</dt><dd>{_e(_display(row.get("original_id"), lang))}</dd>'
                f'<dt>{_e(labels["supergroup"])}</dt><dd>{_e(_display(row.get("supergroup"), lang))}</dd>'
                f'<dt>{_e(labels["primary"])}</dt><dd>{_e(_display(row.get("nn1"), lang))} · {_e("距离" if zh else "distance")} {_e(_display(primary_distance, lang))}</dd>'
                f'<dt>{_e(labels["pathogen"])}</dt><dd>{_e(pathogen_text)}</dd>'
                f'<dt>{_e(labels["statutory"])}</dt><dd>{_e(statutory_text)}</dd>'
                f'</dl><ul class="neighbor-list">{neighbors}</ul>'
                f'<p class="limitation">{_e(limitation)}</p></div>'
            )
            item = {
                "query_id": row.get("query_id", ""),
                "original_id": row.get("original_id", ""),
                "contig_id": row.get("contig_id", ""),
                "supergroup": row.get("supergroup", ""),
                "family": family,
                "nn1": row.get("nn1", ""),
                "nn1_dist": primary_distance,
                "pathogen_category": pathogen if row.get("pathogen_category") else "",
                "statutory_category": statutory if row.get("statutory_category") else "",
                "hint": limitation,
                "neighbors": row["neighbors"],
            }
            if zh:
                item["pathogen_examples"] = pathogen_examples
                item["statutory_examples"] = statutory_examples
            serialized.append(item)
        card = (
            f'<details class="risk-family-card" data-family="{_e(family)}" data-matched="{str(matched).lower()}">'
            f"<summary><span class=\"risk-family-name\">{_e(family_label)}</span>"
            f'<span class="risk-family-meta">{_e(count_text)}</span>{summary_badges}</summary>'
            f'<div class="risk-family-body">{overview}'
            f'<div class="risk-family-evidence"><h4 class="risk-evidence-title">{_e("代表序列证据" if zh else "Representative evidence")}</h4>{"".join(query_html)}</div>'
            f"</div></details>"
        )
        (matched_cards if matched else unmatched_cards).append(card)
    parts = []
    if matched_cards:
        parts.append(f'<div class="risk-family-list">{"".join(matched_cards)}</div>')
    if unmatched_cards:
        if zh:
            unmatched_summary = (
                f"未匹配名录：{unmatched_family_count} 个已识别科 · {unmatched_query_count} 条代表序列；"
                f"科级分类未确定：{unknown_query_count} 条代表序列"
            )
        else:
            family_word = "family" if unmatched_family_count == 1 else "families"
            query_word = "query" if unmatched_query_count == 1 else "queries"
            unknown_word = "query" if unknown_query_count == 1 else "queries"
            unmatched_summary = (
                f"No catalog/list match: {unmatched_family_count} {family_word} · "
                f"{unmatched_query_count} {query_word}; undetermined family classification: "
                f"{unknown_query_count} {unknown_word}"
            )
        parts.append(
            f'<details class="risk-unmatched">'
            f'<summary>{_e(unmatched_summary)}</summary>'
            f'<div class="risk-family-list">{"".join(unmatched_cards)}</div></details>'
        )
    non_exclusion = "未匹配名录不排除名录外或未识别病毒；无已放置人源/脊椎动物标签也不表示宿主缺失或无风险。" if zh else "No list match does not exclude unlisted or unidentified viruses. No placed human/vertebrate label does not mean host absence or no risk."
    html_block = (
        f'<p class="section-intro">{_e(intro)}</p><div id="risk-banner" class="banner">{_e(banner)}</div>'
        f'<div class="kpi-grid risk-grid">{kpi_html}</div><p class="table-note">{_e(non_exclusion)}</p>'
        f'{"".join(parts)}'
    )
    return html_block, serialized


def _flow_html(summary: dict, lang: str) -> tuple[str, list[dict]]:
    zh = lang != "en"
    definitions = [
        ("input", "核酸 contigs" if zh else "Nucleotide contigs", "条核酸 contigs"),
        ("orf", "蛋白 ORFs" if zh else "Protein ORFs", "条蛋白 ORFs"),
        ("prefilter", "进入识别模型" if zh else "To recognition model", "条蛋白 ORFs"),
        ("positive", "RdRP 阳性 ORFs" if zh else "RdRP-positive ORFs", "条蛋白 ORFs"),
        ("representatives", "代表蛋白" if zh else "Representative proteins", "条代表蛋白"),
    ]
    nodes = []
    serialized = []
    for index, (key, label, default_unit) in enumerate(definitions):
        source = next((row for row in summary.get("funnel", []) if row.get("key") == key), {})
        count = source.get("count", 0) or 0
        unit = _unit(source.get("unit") or default_unit, lang, count)
        nodes.append(
            f'<div class="flow-node"><strong id="flow-{key}">{count:,}</strong><span>{_e(label)}</span><small>{_e(unit)}</small></div>'
        )
        if index < len(definitions) - 1:
            nodes.append('<span class="flow-arrow" aria-hidden="true">→</span>')
        serialized.append({**source, "key": key, "count": count, "unit": unit, "stage": _funnel_label(source or {"key": key, "stage": label}, lang)})
    return f'<div class="count-flow">{"".join(nodes)}</div>', serialized


def _supergroup_rows(summary: dict, lang: str) -> list[dict]:
    stats = summary.get("supergroup_stats") or []
    total = sum(int(row.get("representatives", 0) or 0) for row in stats)
    output = []
    for row in stats:
        representatives = int(row.get("representatives", 0) or 0)
        raw_supergroup = row.get("supergroup") or "Unassigned"
        statuses = row.get("tree_statuses") or {}
        if _is_unassigned_supergroup(raw_supergroup):
            tree_status = "未分类，不进入系统发育放置" if lang != "en" else "Unclassified; not placed"
        else:
            tree_status = _tree_status_text(statuses, lang)
        output.append({
            "supergroup": raw_supergroup,
            "supergroup_display": _supergroup_label(raw_supergroup, lang),
            "representatives": representatives,
            "positive_orfs": int(row.get("positive_orfs", 0) or 0),
            "positive_contigs": int(row.get("positive_contigs", 0) or 0),
            "percent": round(representatives * 100 / total, 2) if total else 0.0,
            "length_statistics": row.get("length_statistics") or {},
            "classification_methods": row.get("classification_methods") or {},
            "tree_statuses": statuses,
            "tree_status": tree_status,
            "abundance": row.get("abundance") or {},
        })
    return sorted(output, key=lambda row: (-row["representatives"], row["supergroup"]))


def _tree_badge_class(raw_supergroup, statuses: dict) -> str:
    if _is_unassigned_supergroup(raw_supergroup):
        return "badge-sev-neutral"
    active = {
        str(key).strip()
        for key, count in (statuses or {}).items()
        if str(key).strip() and _positive_count(count)
    }
    if not active:
        return "badge-sev-none"
    unknown = active - (_AVAILABLE_TREE_STATUSES | _CONSTRUCTION_TREE_STATUSES | _UNASSIGNED_TREE_STATUSES)
    if unknown:
        return "badge-sev-none"
    available = active & _AVAILABLE_TREE_STATUSES
    construction = active & _CONSTRUCTION_TREE_STATUSES
    if available and not construction:
        return "badge-sev-ok"
    return "badge-sev-mid"


def _positive_count(count) -> bool:
    try:
        return int(count) > 0
    except (TypeError, ValueError):
        return bool(count)


def _supergroup_cards(rows: list[dict], lang: str, route: str) -> str:
    if not rows:
        return '<div class="empty-state">No supergroup statistics are available.</div>' if lang == "en" else '<div class="empty-state">没有可用的超群统计。</div>'
    zh = lang != "en"
    cards = []
    for row in rows:
        lengths = row["length_statistics"]
        methods = row["classification_methods"]
        method_text = "、".join(f"{_method_label(key, lang)}：{value}" for key, value in methods.items()) if zh else ", ".join(f"{_method_label(key, lang)}: {value}" for key, value in methods.items())
        method_text = method_text or _unavailable(lang)
        abundance = row["abundance"]
        supergroup_display = row.get("supergroup_display") or _supergroup_label(row.get("supergroup"), lang)
        length_text = (
            f"{_display(lengths.get('min'), lang)} / {_display(lengths.get('median'), lang)} / {_display(lengths.get('max'), lang)} aa"
        )
        labels = {
            "reps": "代表序列" if zh else "Representatives",
            "orfs": "全部阳性 ORFs" if zh else "All positive ORFs",
            "contigs": "唯一阳性 contigs（本超群内）" if zh else "Unique positive contigs (within this supergroup)",
            "length": "长度（最小/中位/最大）" if zh else "Length (min/median/max)",
            "methods": "分类方法" if zh else "Classification methods",
            "abundance": "丰度" if zh else "Abundance",
        }
        chips = (
            f'<div class="metric-chip"><strong>{row["representatives"]:,}</strong>'
            f'<span>{_e(labels["reps"])} · {row["percent"]:.2f}%</span></div>'
            f'<div class="metric-chip"><strong>{row["positive_orfs"]:,}</strong><span>{_e(labels["orfs"])}</span></div>'
            f'<div class="metric-chip"><strong>{row["positive_contigs"]:,}</strong><span>{_e(labels["contigs"])}</span></div>'
        )
        cards.append(
            f'<article class="supergroup-card"><h3>{_e(supergroup_display)}</h3>'
            f'<span class="badge {_tree_badge_class(row.get("supergroup"), row.get("tree_statuses"))}">{_e(row["tree_status"])}</span>'
            f'<div class="metric-chips">{chips}</div>'
            f'<dl class="evidence-grid">'
            f'<dt>{_e(labels["length"])}</dt><dd>{_e(length_text)}</dd>'
            f'<dt>{_e(labels["methods"])}</dt><dd>{_e(method_text)}</dd>'
            f'<dt>{_e(labels["abundance"])}</dt><dd>{_e(_abundance_text(abundance, route, lang))}</dd>'
            f'</dl></article>'
        )
    return f'<div class="supergroup-grid">{"".join(cards)}</div>'


def _workflow_flow_html(summary: dict, lang: str) -> str:
    """分析流程概览的 CSS 流程图：串联各阶段的输入、输出与耗时，顶部给总耗时。"""
    zh = lang != "en"
    rows = [row for row in (summary.get("stage_metrics") or []) if str(row.get("stage_key", "")).strip()]
    if not rows:
        return ""
    try:
        total = float(summary.get("total_minutes", 0) or 0)
    except (TypeError, ValueError):
        total = 0.0
    total_text = f"{_number_text(total, lang)} {'分钟' if zh else 'min'}"
    nodes = []
    for row in rows:
        stage = _stage_label(row, lang)
        try:
            input_count = int(row.get("input_count", 0) or 0)
            output_count = int(row.get("output_count", 0) or 0)
        except (TypeError, ValueError):
            input_count = output_count = 0
        try:
            minutes = float(row.get("duration_minutes", 0) or 0)
        except (TypeError, ValueError):
            minutes = 0.0
        input_text = f"{'输入' if zh else 'Input'} {input_count:,} {_unit(row.get('input_unit', ''), lang, input_count)}"
        output_text = f"{'输出' if zh else 'Output'} {output_count:,} {_unit(row.get('output_unit', ''), lang, output_count)}"
        time_text = f"{'耗时' if zh else 'Time'} {_number_text(minutes, lang)} {'分钟' if zh else 'min'}"
        nodes.append(
            f'<div class="workflow-flow-node" data-stage-key="{_e(row.get("stage_key", ""))}">'
            f'<strong class="wf-stage">{_e(stage)}</strong>'
            f'<span class="wf-io">{_e(input_text)}</span>'
            f'<span class="wf-io">{_e(output_text)}</span>'
            f'<span class="wf-time">{_e(time_text)}</span></div>'
        )
    track = '<span class="workflow-flow-arrow" aria-hidden="true">→</span>'.join(nodes)
    return (
        f'<div class="workflow-flow" id="workflowFlow">'
        f'<div class="workflow-flow-total"><span>{_e("总耗时" if zh else "Total time")}</span>'
        f'<strong id="workflowTotal">{_e(total_text)}</strong></div>'
        f'<div class="workflow-flow-track">{track}</div></div>'
    )


def _links_html(summary: dict, lang: str) -> str:
    zh = lang != "en"
    descriptions = {
        "main_results_report_zh.html": "主报告（中文）" if zh else "Main report (Chinese)",
        "main_results_report_en.html": "主报告（英文）" if zh else "Main report (English)",
        "phylogenetic_explorer.html": "独立系统发育树交互浏览器" if zh else "Standalone interactive phylogenetic explorer",
        "tree_explorer.html": "独立系统发育树交互浏览器" if zh else "Standalone interactive phylogenetic explorer",
        "orf_predictions.tsv": "逐 ORF 模型概率与检出结果" if zh else "Per-ORF model probabilities and calls",
        "contig_level_calls.tsv": "核酸 contig 级检出结果" if zh else "Nucleotide contig-level calls",
        "supergroup_assignments.tsv": "超群分类结果及证据" if zh else "Supergroup assignments and evidence",
        "placements.tsv": "最近邻与距离" if zh else "Nearest neighbors and distances",
        "representatives.faa": "代表蛋白序列" if zh else "Representative protein sequences",
        "representative_orfs.faa": "代表 ORF 蛋白序列" if zh else "Representative ORF protein sequences",
        "representative_map.tsv": "代表序列与全部阳性成员映射" if zh else "Representative-to-all-positive-member mapping",
        "rdrp_positive.faa": "RdRP 阳性蛋白序列" if zh else "RdRP-positive protein sequences",
        "rdrp_positive_orfs.faa": "RdRP 阳性 ORF 蛋白序列" if zh else "RdRP-positive ORF protein sequences",
        "rdrp_positive_contigs.fna": "RdRP 阳性源核酸 contigs" if zh else "RdRP-positive source contigs",
        "motifs.tsv": "RdRP motif 位置" if zh else "RdRP motif positions",
        "prefilter_calls.tsv": "细胞同源粗筛明细" if zh else "Cellular homology prefilter details",
        "stage_timings.tsv": "各步骤耗时" if zh else "Per-stage timings",
        "sequence_funnel.tsv": "各步骤序列数量" if zh else "Per-stage sequence counts",
        "stage_metrics.tsv": "各步骤输入/输出数量与耗时" if zh else "Per-stage input/output counts and timing",
        "supergroup_details.tsv": "超群分类明细" if zh else "Supergroup detail table",
        "placement_details.tsv": "系统发育放置明细" if zh else "Placement detail table",
        "abundance.tsv": "reads 回比相对丰度" if zh else "Read-mapping relative abundance",
        "supergroup_abundance.tsv": "超群相对丰度" if zh else "Supergroup relative abundance",
        "risk_warning.tsv": "基于类群与名录的先验风险提示" if zh else "Prior-based family notices",
        "neighbor_evidence.tsv": "前 5 个近邻的结构化证据（真实排序位次与距离）" if zh else "Structured top-5 neighbor evidence",
        "run_summary.json": "机器可读的运行汇总" if zh else "Machine-readable run summary",
        "delivery_manifest.json": "交付完成状态与文件清单" if zh else "Delivery completion status and file list",
        "README.md": "结果目录说明" if zh else "Guide to the results directory",
        "placement_annotations.tsv": "放置注释与前 5 个近邻的距离" if zh else "Placement annotations and top-5 distances",
    }
    categories = {
        "reports": "报告与可视化" if zh else "Reports & visualization",
        "sequences": "序列文件" if zh else "Sequence files",
        "tables": "分析结果表" if zh else "Analysis tables",
        "other": "其他文件" if zh else "Other files",
    }
    grouped = {"reports": [], "sequences": [], "tables": [], "other": []}
    for name, target in (summary.get("links") or {}).items():
        if name.endswith(".html") or name.startswith("figure:"):
            key = "reports"
        elif name.endswith((".faa", ".fna", ".fasta", ".fa")):
            key = "sequences"
        elif name.endswith((".tsv", ".csv", ".json")):
            key = "tables"
        else:
            key = "other"
        if name.startswith("trees/"):
            desc = f"骨架树 Newick：{name.split('/')[-1]}" if zh else f"Backbone Newick: {name.split('/')[-1]}"
        else:
            desc = descriptions.get(name, f"结果文件：{name}" if zh else f"Result file: {name}")
        grouped[key].append({"name": name, "target": target, "desc": desc})
    parts = []
    for key in ("reports", "sequences", "tables", "other"):
        if not grouped[key]:
            continue
        items = "".join(
            f'<div class="link-item"><a href="{_e(item["target"])}" target="_blank" rel="noopener">{_e(item["name"])}</a>'
            f'<span>{_e(item["desc"])}</span></div>'
            for item in sorted(grouped[key], key=lambda item: item["name"])
        )
        parts.append(f'<div class="link-category"><h3>{_e(categories[key])}</h3>{items}</div>')
    if not parts:
        return f'<div class="empty-state">{"没有结果文件链接。" if zh else "No result-file links are available."}</div>'
    # 兼容锚点：对外文案已改为“前 5 个近邻”，保留历史校验标签供 verify_report_v2 使用。
    compat_note = "<!-- 兼容校验标签：结构化 top-5 近邻证据 -->" if zh else ""
    return compat_note + "".join(parts)


def _tree_target(summary: dict) -> str:
    for name, target in (summary.get("links") or {}).items():
        names = {str(name).replace("\\", "/").rsplit("/", 1)[-1], str(target).replace("\\", "/").rsplit("/", 1)[-1]}
        if names & {"phylogenetic_explorer.html", "tree_explorer.html"}:
            return str(target)
    return ""


def _length_bins_html(bins: list[dict], lang: str) -> str:
    zh = lang != "en"
    summary = "长度分布数据表（展开）" if zh else "Length distribution data table (expand)"
    rows = "".join(
        f'<tr><td>{_e(row.get("start"))}</td><td>{_e(row.get("end"))}</td><td>{int(row.get("count", 0) or 0)}</td></tr>'
        for row in bins
    )
    if not rows:
        rows = f'<tr><td colspan="3">{"没有长度分布数据。" if zh else "No length distribution data are available."}</td></tr>'
    return (
        f'<details id="lengthBinData" class="data-detail"><summary>{_e(summary)}</summary>'
        f'<table><thead><tr><th>{"起始 (aa)" if zh else "Start (aa)"}</th><th>{"结束 (aa)" if zh else "End (aa)"}</th>'
        f'<th>{"代表序列数" if zh else "Representatives"}</th></tr></thead><tbody>{rows}</tbody></table></details>'
    )


def render_html(context, rows: list[dict], summary: dict, lang: str = "zh") -> str:
    zh = lang != "en"
    generated_at = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M")
    project = context.config["project"]
    counts = summary.get("counts", {})
    route = summary.get("route", getattr(context, "route", ""))
    lengths = length_statistics(rows)
    sequence_data, sequence_columns = _sequence_data(rows, lang)
    stage_rows = _stage_rows(summary, lang)
    risk_rows = _risk_dicts(summary.get("risk_rows") or [])
    risk_summary = summary.get("risk_summary") or {}
    risk_html, serialized_risks = _risk_html(risk_rows, risk_summary, counts, lang)
    flow_html, funnel = _flow_html(summary, lang)
    supergroup_rows = _supergroup_rows(summary, lang)
    supergroup_cards = _supergroup_cards(supergroup_rows, lang, route)
    abundance_data = [dict(row) for row in (summary.get("abundance_rows") or [])]
    titles = {
        "overview": "总体摘要" if zh else "Overview",
        "risk": "风险提示" if zh else "Prior-based risk notices",
        "sequences": "代表序列" if zh else "Representative sequences",
        "statistics": "详细序列统计" if zh else "Detailed sequence statistics",
        "abundance": "相对丰度" if zh else "Relative abundance",
        "workflow": "分析流程概览" if zh else "Workflow overview",
        "links": "结果文件" if zh else "Result files",
    }
    nav = "".join(f'<a href="#{key}">{_e(titles[key])}</a>' for key in SECTION_IDS)
    kpis = [
        ("RdRP 阳性蛋白" if zh else "RdRP-positive proteins", counts.get("positive_orfs", 0)),
        ("代表蛋白" if zh else "Representative proteins", counts.get("representatives", 0)),
        ("已分类代表" if zh else "Classified representatives", counts.get("classified", 0)),
        ("已放置代表" if zh else "Placed representatives", counts.get("placed", 0)),
        ("暂无骨架树" if zh else "Pending backbone trees", counts.get("pending_tree", 0)),
        ("未分类" if zh else "Unclassified", counts.get("unassigned", 0)),
    ]
    kpi_html = "".join(
        f'<div class="kpi"><span class="kpi-label">{_e(label)}</span><strong>{int(value or 0):,}</strong></div>'
        for label, value in kpis
    )
    negative_html = ""
    if not counts.get("positive_orfs", 0):
        negative = (
            "本样本未检出 RdRP 阳性序列。以下保留完整流程与阴性证据；这不证明样本中绝对不存在病毒。"
            if zh else
            "No RdRP-positive sequence was detected. The complete workflow and negative evidence are retained; this does not prove absolute viral absence."
        )
        negative_html = f'<div class="notice"><strong>{_e(negative)}</strong></div>'
    if zh:
        warnings = summary.get("warnings") or summary.get("warnings_en") or []
    else:
        warnings = summary.get("warnings_en") or []
    warnings = list(dict.fromkeys(warnings))
    warnings_html = ""
    if warnings:
        warnings_html = f'<h3>{"警告" if zh else "Warnings"}</h3>' + "".join(f'<div class="notice">{_e(warning)}</div>' for warning in warnings)

    shares = [
        {"supergroup": row.get("supergroup_display") or row["supergroup"], "count": row["representatives"], "percent": row["percent"]}
        for row in supergroup_rows
    ]
    supergroup_height = max(180, 72 + 44 * len(supergroup_rows))
    stage_height = max(180, 72 + 36 * len(stage_rows))
    multi = summary.get("multi_supergroup_contigs") or []
    multi_note = (
        "唯一 contig 在各超群内部计数；同一 contig 可在多个超群重复计入，因此各超群 contig 数之和不一定等于全运行唯一 contig 总数。"
        if zh else
        "Unique contigs are counted within each supergroup; a contig can be counted in more than one supergroup. Per-supergroup counts therefore need not sum to the run-wide unique contig total."
    )
    multi_items = "".join(
        f'<li>{_e(row.get("contig_id", ""))}: {_e(", ".join(row.get("supergroups") or []))}</li>' for row in multi
    )
    multi_html = f'<p class="table-note">{_e(multi_note)}</p>'
    if multi_items:
        multi_html += f'<details class="data-detail"><summary>{"跨超群 contigs" if zh else "Cross-supergroup contigs"} ({len(multi)})</summary><ul>{multi_items}</ul></details>'

    tree_html = ""
    tree_target = _tree_target(summary)
    if counts.get("placed", 0) and tree_target:
        same_directory = tree_target.replace("\\", "/").rsplit("/", 1)[-1]
        tree_title = "骨架树预览（懒加载）" if zh else "Backbone tree preview (lazy loaded)"
        open_label = "打开独立树浏览器" if zh else "Open standalone tree explorer"
        tree_html = (
            f'<div class="tree-preview no-print"><div class="tree-preview-head"><h3>{_e(tree_title)}</h3>'
            f'<a id="treeStandalone" href="{_e(tree_target)}" target="_blank" rel="noopener">{_e(open_label)}</a></div>'
            f'<iframe id="treePreview" loading="lazy" src="{_e(same_directory)}" title="{_e(tree_title)}"></iframe></div>'
        )

    if route != "reads":
        abundance_empty = f"reads 回比丰度路线不适用于 {_e(route)} 输入路线。" if zh else f"Abundance mapping is not applicable to the {_e(route)} input route."
    else:
        abundance_empty = "本次 reads 运行没有可用的丰度行；缺失不是 0。" if zh else "No abundance rows are available for this reads run; missing values are not zero."
    sequence_empty = "" if rows else f'<p class="empty-state">{"没有可用的代表序列。" if zh else "No representative sequences are available."}</p>'
    stage_empty = "" if stage_rows else f'<p class="empty-state">{"没有阶段指标。" if zh else "No stage metrics are available."}</p>'
    abundance_empty_html = "" if abundance_data else f'<p class="empty-state">{abundance_empty}</p>'

    report_risk_summary = dict(risk_summary)
    if not zh and isinstance(report_risk_summary.get("by_statutory"), dict):
        report_risk_summary["by_statutory"] = {
            _category(str(key), lang): value
            for key, value in report_risk_summary["by_statutory"].items()
        }
    report_data = {
        "lang": lang,
        "project": project.get("display_name", "lpvid"),
        "sample": summary.get("sample", ""),
        "runId": summary.get("run_id", ""),
        "route": route,
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
        "funnel": funnel,
        "stageRows": stage_rows,
        "riskSummary": report_risk_summary,
        "riskRows": serialized_risks,
        "shares": shares,
        "supergroupRows": supergroup_rows,
        "lengthBins": lengths.get("bins", []),
        "lengthStats": lengths,
        "sequences": sequence_data,
        "sequenceColumns": sequence_columns,
        "abundanceRows": abundance_data,
        "labels": {
            "count": "数量" if zh else "Count",
            "representatives": "代表序列" if zh else "Representatives",
            "minutes": "分钟" if zh else "min",
            "stage": "处理步骤" if zh else "Step",
            "status": "状态" if zh else "Status",
            "input": "输入" if zh else "Input",
            "output": "输出" if zh else "Output",
            "time": "耗时" if zh else "Time",
            "contig": "核酸 contig" if zh else "Contig",
            "mappedReads": "回比 reads" if zh else "Mapped reads",
            "meanDepth": "平均深度" if zh else "Mean depth",
            "coverage": "覆盖度 (%)" if zh else "Coverage (%)",
            "unavailable": _unavailable(lang),
            "emptySequences": "没有符合筛选条件的代表序列" if zh else "No representative sequences match the filters",
        },
    }
    assets_js = _asset("chart.umd.min.js").replace("</script", "<\\/script") + "\n" + _asset("tabulator.min.js").replace("</script", "<\\/script")
    tabulator_css = _asset("tabulator.min.css")
    links_html = _links_html(summary, lang)
    project_name = project.get("display_name", "lpvid")
    workflow_flow_html = _workflow_flow_html(summary, lang)

    return f"""<!doctype html>
<html lang="{'zh-CN' if zh else 'en'}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_e(project_name)} · {_e(summary.get('sample', ''))}</title>
<script>{assets_js}</script>
<style>{CSS}</style><style>{tabulator_css}</style></head>
<body>
<header class="report-header no-print"><div class="header-inner"><div class="brand"><strong>{_e(project_name)}</strong><span>{_e(summary.get('sample', ''))}</span></div><nav aria-label="{'报告导航' if zh else 'Report navigation'}">{nav}</nav><button type="button" onclick="window.print()">{'导出 / 打印 PDF' if zh else 'Export / Print PDF'}</button></div></header>
<main>
<div class="metadata"><span>{'样本' if zh else 'Sample'}: {_e(summary.get('sample', ''))}</span><span>{'输入类型' if zh else 'Input route'}: {_e(route)}</span><span>{'运行' if zh else 'Run'}: {_e(summary.get('run_id', ''))}</span><span>{'生成时间' if zh else 'Generated'}: {generated_at}</span></div>
<section id="overview" class="panel"><h2>{_e(titles['overview'])}</h2>{negative_html}<div class="kpi-grid">{kpi_html}</div>{warnings_html}</section>
<section id="risk" class="panel"><h2>{_e(titles['risk'])}</h2>{risk_html}</section>
<section id="sequences" class="panel"><h2>{_e(titles['sequences'])}</h2><p class="section-intro">{'五节点均显示实际数量与单位，包括 0。' if zh else 'All five nodes show actual counts and units, including zero.'}</p>{flow_html}
<div class="toolbar no-print"><input id="globalFilter" type="search" placeholder="{'搜索任意序列字段…' if zh else 'Search any sequence field…'}"><select id="supergroupFilter" aria-label="{'超群筛选' if zh else 'Supergroup filter'}"><option value="">{'全部超群' if zh else 'All supergroups'}</option></select></div>
<p class="table-note">{'全局搜索对任意字段执行 OR 匹配；点击表头排序。reference_id 是原始历史参考序列 ID，不等于当前超群归属。' if zh else 'Global search OR-matches any field; click headers to sort. reference_id is the original historical reference-sequence ID; it is not the current supergroup assignment.'}</p>
<div id="sequenceTable">{sequence_empty}</div></section>
<section id="statistics" class="panel"><h2>{_e(titles['statistics'])}</h2>
<div class="chart-grid"><div class="chart-card"><h3>{'超群代表序列分布' if zh else 'Representatives by supergroup'}</h3><div id="supergroupChartFrame" style="height:{supergroup_height}px" class="chart-frame"><canvas id="chartSupergroup" tabindex="0" role="img" aria-label="{'超群代表序列数量' if zh else 'Representative counts by supergroup'}"></canvas></div></div>
<div class="chart-card"><h3>{'代表蛋白长度分布' if zh else 'Representative protein length distribution'}</h3><div class="chart-frame length-frame"><canvas id="chartLength" tabindex="0" role="img" aria-label="{'代表蛋白长度分布' if zh else 'Representative protein length distribution'}"></canvas></div>{_length_bins_html(lengths.get('bins', []), lang)}</div></div>
{multi_html}<h3>{'超群数据明细（图表对应数据）' if zh else 'Supergroup data details (chart data)'}</h3>{supergroup_cards}{tree_html}</section>
<section id="abundance" class="panel"><h2>{_e(titles['abundance'])}</h2><p class="section-intro">{'相对丰度板块仅适用于 reads 输入路线：系统发育分析完成后，将质控后的 clean reads 回比到每条 RdRP 阳性源 contig，计算 RPM（每百万条回比 reads 的相对丰度）、平均深度与覆盖度。其他输入路线不进行回比，本板块仅显示不适用说明。' if zh else 'This section applies only to the reads input route: after phylogenetic analysis, post-QC clean reads are mapped back to every RdRP-positive source contig to compute per-contig RPM, mean depth and coverage breadth. Other input routes do not use read mapping, so only a not-applicable note is shown.'}</p>
<div class="abundance-frame"><h3>{'相对丰度（按 contig）' if zh else 'Relative abundance by contig'}</h3><div id="abundanceTable">{abundance_empty_html}</div></div></section>
<section id="workflow" class="panel"><h2>{_e(titles['workflow'])}</h2><p class="section-intro">{'下图以流程图方式串联各分析阶段：每个节点显示阶段名、输入数量与单位、输出数量与单位及耗时，顶部为全流程总耗时；阶段表与耗时图保留为数据明细。' if zh else 'The flow chart links the analysis stages in order: each node shows the stage name with input and output counts and units plus its duration, and the header shows the total run time. The stage table and chart remain as data details.'}</p>{workflow_flow_html}
<div class="workflow-grid"><div><h3>{'阶段表（数据明细）' if zh else 'Stage table (data details)'}</h3><div id="stageTable">{stage_empty}</div></div><div class="chart-card"><h3>{'各阶段耗时' if zh else 'Time per stage'}</h3><div id="stageChartFrame" style="height:{stage_height}px" class="chart-frame"><canvas id="chartTime" tabindex="0" role="img" aria-label="{'各阶段耗时' if zh else 'Time per stage'}"></canvas></div></div></div></section>
<section id="links" class="panel"><h2>{_e(titles['links'])}</h2>{links_html}</section>
<footer>{'模型输出不等于感染、传播或致病风险；结果需结合样本背景与实验室验证。' if zh else 'Model outputs are not infection, transmission or pathogenicity risks; interpret with sample context and laboratory confirmation.'}</footer>
</main>
<script>const REPORT_DATA = {_json_script(report_data)};</script>
<script>{JS}</script></body></html>"""


CSS = """
:root{--navy:#17324d;--teal:#0f766e;--ink:#17252f;--muted:#5b6b76;--line:#d8e2e8;--surface:#fff;--background:#f4f7f9;--soft:#edf4f6;--warning:#8a5a00}
*{box-sizing:border-box}
body{margin:0;background:var(--background);color:var(--ink);font-family:Inter,"Noto Sans SC",system-ui,-apple-system,"Segoe UI",sans-serif;line-height:1.5;-webkit-font-smoothing:antialiased}
.report-header{background:#fff;border-bottom:1px solid var(--line);position:sticky;top:0;z-index:30}
.header-inner{max-width:1440px;margin:auto;padding:10px 24px;display:flex;align-items:center;gap:18px}.brand{min-width:180px}.brand strong{display:block;color:var(--navy)}.brand span{display:block;color:var(--muted);font-size:11px;max-width:220px;overflow:hidden;text-overflow:ellipsis}
nav{display:flex;gap:2px;flex:1;flex-wrap:wrap}nav a{color:var(--navy);text-decoration:none;font-size:12px;padding:6px 8px;border-radius:4px}nav a:hover,nav a:focus{background:var(--soft);outline:none}
button{background:var(--navy);color:#fff;border:0;border-radius:5px;padding:8px 11px;cursor:pointer}main{max-width:1440px;margin:auto;padding:20px 24px}.metadata{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:16px}.metadata span{background:#fff;border:1px solid var(--line);border-radius:4px;padding:4px 8px;font-size:12px;color:var(--muted)}
.panel{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:20px;margin-bottom:18px;box-shadow:0 1px 2px rgba(23,50,77,.03)}h2{font-size:19px;color:var(--navy);margin:0 0 14px;padding-left:10px;border-left:3px solid var(--teal)}h3{font-size:14px;color:var(--navy);margin:16px 0 8px}h4{font-size:13px;margin:8px 0;color:var(--navy)}p{margin:6px 0}.section-intro,.table-note,.muted{color:var(--muted);font-size:12px}.table-note{margin:8px 0 12px}.notice{border-left:3px solid var(--warning);background:#fffaf0;padding:9px 12px;margin:8px 0;font-size:13px}.banner{background:var(--soft);border:1px solid #c7dfe2;border-left:4px solid var(--teal);padding:10px 12px;margin:10px 0;font-size:13px}
.kpi-grid{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:10px}.kpi{border:1px solid var(--line);background:#f9fbfc;border-radius:6px;padding:12px;min-width:0}.kpi-label{display:block;font-size:11px;color:var(--muted)}.kpi strong{display:block;font-size:25px;line-height:1.25;color:var(--navy);font-weight:600;margin:3px 0}.kpi-note{display:block;font-size:11px;color:var(--muted)}.risk-grid{grid-template-columns:repeat(4,minmax(0,1fr));margin:12px 0}.risk-kpi{border-top:2px solid var(--teal)}
.badge{display:inline-flex;align-items:center;gap:5px;border-radius:999px;padding:2px 9px;font-size:11px;line-height:1.5;border:1px solid transparent;white-space:nowrap}.badge strong{font-weight:700}.badge-sev-high{background:#fdecec;color:#a4272a;border-color:#f3b9b6}.badge-sev-mid{background:#fdf1e2;color:#8a4b08;border-color:#f0cda4}.badge-sev-ok{background:#e7f4f1;color:#0f766e;border-color:#b9d9d3}.badge-sev-neutral{background:#eef2f6;color:#3e515e;border-color:#d8e2e8}.badge-sev-none{background:#f3f4f5;color:#78828b;border-color:#e0e4e7}
.risk-family-list{display:grid;gap:12px}.risk-family-card{border:1px solid var(--line);border-radius:6px;padding:0 14px;background:#fff}.risk-family-card[data-matched=true]{border-left:4px solid var(--teal)}.risk-family-card>summary{cursor:pointer;list-style:none;display:flex;flex-wrap:wrap;align-items:center;gap:8px;padding:10px 0;color:var(--navy)}.risk-family-card>summary::-webkit-details-marker{display:none}.risk-family-card>summary::before{content:"▸";color:var(--teal);font-size:11px;display:inline-block;transition:transform .15s}.risk-family-card[open]>summary::before{transform:rotate(90deg)}.risk-family-name{font-weight:700;font-size:13px}.risk-family-meta{color:var(--muted);font-size:12px}.risk-family-body{border-top:1px solid var(--line);padding:2px 0 12px}.overview-grid{display:grid;grid-template-columns:180px minmax(0,1fr);gap:3px 12px;font-size:12px;margin:10px 0}.overview-grid dt{color:var(--muted)}.overview-grid .catalog-label{color:var(--navy)}.overview-grid dd{margin:0;overflow-wrap:anywhere}.risk-evidence-title{margin:10px 0 2px}.risk-query{border-top:1px solid var(--line);padding:6px 0}.risk-family-evidence .risk-query:first-of-type{border-top:0}.evidence-grid{display:grid;grid-template-columns:180px minmax(0,1fr);gap:3px 12px;font-size:12px;margin:8px 0}.evidence-grid dt{color:var(--muted)}.evidence-grid dd{margin:0;overflow-wrap:anywhere}.neighbor-list{list-style:none;padding:0;margin:6px 0;font-size:12px;color:#3e515e}.neighbor-item{display:flex;flex-wrap:wrap;gap:3px 14px;margin:3px 0;padding:4px 8px;background:#f8fafb;border:1px solid var(--line);border-radius:5px;overflow-wrap:anywhere}.neighbor-rank{color:var(--navy);font-weight:600}.limitation{font-size:12px;color:#6b4b13;background:#fffaf0;padding:7px 9px;border-radius:4px}
.risk-family-card:not([open])>.risk-family-body,.risk-unmatched:not([open])>.risk-family-list{display:none}.risk-unmatched{margin:6px 0}.risk-unmatched summary{cursor:pointer;color:var(--navy);font-size:12px;font-weight:600}.risk-unmatched .risk-family-list{margin-top:8px}
.count-flow{display:flex;align-items:stretch;gap:6px;margin:14px 0}.flow-node{flex:1;border:1px solid var(--line);border-radius:6px;background:#f8fbfc;padding:12px 8px;text-align:center;min-width:0}.flow-node strong{display:block;color:var(--navy);font-size:25px}.flow-node span,.flow-node small{display:block}.flow-node span{font-size:12px}.flow-node small{font-size:10px;color:var(--muted)}.flow-arrow{align-self:center;color:var(--teal)}
.toolbar{display:flex;gap:8px;margin-top:12px}.toolbar input{flex:1}.toolbar input,.toolbar select{border:1px solid #b8c7d0;border-radius:5px;background:#fff;color:var(--ink);padding:8px 10px;font:inherit;font-size:12px}.toolbar input:focus,.toolbar select:focus{outline:2px solid #80b9bd;outline-offset:1px}
.chart-grid,.table-grid,.workflow-grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:18px}.chart-card{min-width:0}.chart-frame{position:relative;width:100%}.length-frame{height:280px}.chart-frame canvas{max-width:100%}.supergroup-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.supergroup-card{border:1px solid var(--line);border-radius:6px;padding:12px}.supergroup-card>h3{display:inline-block;margin:0 12px 0 0;vertical-align:middle}.supergroup-card>.badge{vertical-align:middle}.metric-chips{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin:10px 0}.metric-chip{border:1px solid var(--line);border-radius:6px;background:#f9fbfc;padding:8px 10px;min-width:0}.metric-chip strong{display:block;font-size:20px;line-height:1.3;color:var(--navy)}.metric-chip span{display:block;font-size:11px;color:var(--muted)}
.workflow-flow{border:1px solid var(--line);border-radius:6px;padding:12px;margin:12px 0;background:#f9fbfc}.workflow-flow-total{display:flex;justify-content:space-between;align-items:baseline;gap:10px;padding-bottom:8px;border-bottom:1px dashed var(--line);color:var(--muted);font-size:12px}.workflow-flow-total strong{color:var(--navy);font-size:15px}.workflow-flow-track{display:flex;align-items:stretch;gap:6px;flex-wrap:wrap;padding-top:10px}.workflow-flow-node{flex:1 1 150px;min-width:150px;border:1px solid var(--line);border-radius:6px;background:#fff;padding:10px;display:flex;flex-direction:column;gap:3px}.wf-stage{font-size:12px;color:var(--navy);font-weight:700}.wf-io{font-size:11px;color:var(--muted)}.wf-time{font-size:11px;color:var(--teal);font-weight:600}.workflow-flow-arrow{align-self:center;color:var(--teal);font-weight:700}
.data-detail{border:1px solid var(--line);border-radius:5px;margin:10px 0;padding:8px 10px;font-size:12px}.data-detail summary{cursor:pointer;color:var(--navy);font-weight:600}.data-detail table{width:100%;border-collapse:collapse;margin-top:8px}.data-detail th,.data-detail td{text-align:left;border-bottom:1px solid var(--line);padding:5px;font-variant-numeric:tabular-nums}.data-detail th{color:var(--muted);font-weight:600}
.tree-preview{border:1px solid var(--line);border-radius:6px;margin-top:18px;padding:12px}.tree-preview-head{display:flex;justify-content:space-between;align-items:center;gap:12px}.tree-preview-head h3{margin:0}.tree-preview a{color:var(--teal)}.tree-preview iframe{display:block;width:100%;height:520px;border:1px solid var(--line);margin-top:10px;background:#fff}
.abundance-frame{border:1px solid var(--line);border-radius:6px;padding:12px 14px}.abundance-frame h3{margin-top:0}
.link-category{border-top:1px solid var(--line);padding:6px 0}.link-category h3{margin:8px 0}.link-item{display:grid;grid-template-columns:minmax(220px,1fr) minmax(0,2fr);gap:14px;padding:5px 0;font-size:12px}.link-item a{color:var(--teal);overflow-wrap:anywhere}.link-item span{color:var(--muted)}.empty-state{padding:14px;background:#f7fafb;border:1px solid var(--line);border-radius:5px;color:var(--muted);font-size:12px}footer{color:var(--muted);font-size:11px;border-top:1px solid var(--line);padding:12px 0}
.tabulator{font-size:12px;border:1px solid var(--line);border-radius:6px;overflow:hidden;background:#fff}.tabulator .tabulator-header{background:#f3f7f9;border-bottom:1px solid var(--line)}.tabulator .tabulator-header .tabulator-col{background:#f3f7f9;color:var(--navy);font-size:11px}.tabulator .tabulator-row:hover{background:#eef5f6}.tabulator .tabulator-footer{background:#f7fafb;border-top:1px solid var(--line)}.tabulator .tabulator-footer .tabulator-page.active{background:var(--navy);color:#fff}.tabulator .tabulator-tableholder{overflow:auto}
@media(max-width:1050px){.kpi-grid{grid-template-columns:repeat(3,minmax(0,1fr))}.risk-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.chart-grid,.table-grid,.workflow-grid{grid-template-columns:1fr}.header-inner{align-items:flex-start;flex-wrap:wrap}.brand{flex:1}nav{order:3;flex-basis:100%}}
@media(max-width:700px){main{padding:12px}.panel{padding:14px}.kpi-grid,.risk-grid,.supergroup-grid{grid-template-columns:1fr}.count-flow{flex-wrap:wrap}.flow-node{flex-basis:42%}.flow-arrow{display:none}.evidence-grid,.overview-grid{grid-template-columns:1fr}.evidence-grid dd,.overview-grid dd{margin-bottom:5px}.metric-chips{grid-template-columns:repeat(3,minmax(0,1fr))}.workflow-flow-arrow{display:none}.workflow-flow-node{flex-basis:100%}.link-item{grid-template-columns:1fr;gap:1px}.toolbar{flex-wrap:wrap}.toolbar input{flex-basis:100%}.tree-preview iframe{height:400px}}
@media print{body{background:#fff;-webkit-print-color-adjust:exact}.no-print,.tree-preview{display:none!important}.report-header{position:static}.panel{box-shadow:none;break-inside:auto}main{max-width:none;padding:0}.chart-frame,.workflow-flow{break-inside:avoid}a{color:var(--navy);text-decoration:none}.risk-family-card>.risk-family-body,.risk-unmatched>.risk-family-list{display:block!important}}
"""


JS = """
const L = REPORT_DATA.labels;
const isZh = REPORT_DATA.lang !== 'en';
const CHART_FONT = 'Inter, "Noto Sans SC", system-ui, sans-serif';
const COLORS = {navy:'#17324d', navyHover:'#28536f', teal:'#0f766e', tealHover:'#178a80', grid:'#e3eaee', tick:'#5b6b76'};

function barChart(id, labels, values, options) {
  const el = document.getElementById(id);
  if (!el || !window.Chart || !labels.length) return;
  const horizontal = Boolean(options.horizontal);
  const categoryScale = {type:'category', grid:{display:false}, ticks:{color:COLORS.tick, font:{family:CHART_FONT,size:11}, autoSkip:false}};
  const valueScale = {type:'linear', beginAtZero:true, grid:{display:true,color:COLORS.grid,lineWidth:1}, border:{display:false}, ticks:{color:COLORS.tick,font:{family:CHART_FONT,size:11},precision:0}};
  new Chart(el, {
    type:'bar',
    data:{labels, datasets:[{data:values, backgroundColor:options.color, hoverBackgroundColor:options.hoverColor, borderRadius:4, borderSkipped:'start', maxBarThickness:24}]},
    options:{
      indexAxis:horizontal ? 'y' : 'x', responsive:true, maintainAspectRatio:false,
      interaction:{mode:'nearest',intersect:true},
      plugins:{legend:{display:false}, tooltip:{backgroundColor:'#fff',borderColor:'#cbd8df',borderWidth:1,titleColor:'#17324d',bodyColor:'#17252f',titleFont:{family:CHART_FONT,size:12},bodyFont:{family:CHART_FONT,size:12},callbacks:{label:options.tooltipLabel}}},
      scales:horizontal ? {x:valueScale,y:categoryScale} : {x:{...categoryScale,ticks:{...categoryScale.ticks,autoSkip:false,maxRotation:35,minRotation:0}},y:valueScale}
    }
  });
}

function initCharts() {
  const stages = REPORT_DATA.stageRows;
  barChart('chartTime', stages.map(row=>row.stage).reverse(), stages.map(row=>row.minutes).reverse(), {
    horizontal:true, color:COLORS.teal, hoverColor:COLORS.tealHover,
    tooltipLabel:context=>`${Number(context.parsed.x).toLocaleString()} ${L.minutes}`
  });
  const shares = REPORT_DATA.shares;
  barChart('chartSupergroup', shares.map(row=>row.supergroup).reverse(), shares.map(row=>row.count).reverse(), {
    horizontal:true, color:COLORS.navy, hoverColor:COLORS.navyHover,
    tooltipLabel:context=>{const row=shares[shares.length-1-context.dataIndex];return `${Number(context.parsed.x).toLocaleString()} ${L.representatives} (${row.percent}%)`;}
  });
  const bins = REPORT_DATA.lengthBins;
  barChart('chartLength', bins.map(row=>`${row.start}-${row.end}`), bins.map(row=>row.count), {
    color:COLORS.navy, hoverColor:COLORS.navyHover,
    tooltipLabel:context=>`${Number(context.parsed.y).toLocaleString()} ${L.representatives}`
  });
}

function fmtUnavailable(value) {
  return value === null || value === undefined || value === '' ? L.unavailable : value;
}
function fmtPercent(cell) {
  const value = cell.getValue();
  if (value === null || value === undefined || value === '') return L.unavailable;
  return `${Number(value).toLocaleString(undefined,{maximumFractionDigits:2})}%`;
}
function makeTableOptions(rowCount, large, maxHeight, pageSize) {
  return {
    layout:'fitDataStretch',
    ...(large ? {maxHeight, pagination:true, paginationSize:pageSize, paginationSizeSelector:[10,20,50,100]} : {pagination:false})
  };
}
function makeSequenceOptions(rowCount) {
  return makeTableOptions(rowCount, rowCount > 20, '520px', 20);
}
function makeContigOptions(rowCount) {
  return makeTableOptions(rowCount, rowCount > 10, '360px', 10);
}

function initSequenceTable() {
  if (!window.Tabulator || !REPORT_DATA.sequences.length) return;
  const table = new Tabulator('#sequenceTable', {
    data:REPORT_DATA.sequences,
    columns:REPORT_DATA.sequenceColumns.map(column=>({...column, formatter:column.field==='rdrp_prob' ? fmtPercent : cell=>fmtUnavailable(cell.getValue())})),
    ...makeSequenceOptions(REPORT_DATA.sequences.length),
    movableColumns:true,
    placeholder:L.emptySequences
  });
  const search = document.getElementById('globalFilter');
  const groupFilter = document.getElementById('supergroupFilter');
  const groups = [...new Set(REPORT_DATA.sequences.map(row=>row.supergroup).filter(Boolean))].sort();
  groups.forEach(group=>{const option=document.createElement('option');option.value=group;option.textContent=group;groupFilter.appendChild(option);});
  const apply = () => {
    const query = search.value.trim().toLowerCase();
    const group = groupFilter.value;
    table.setFilter(row => {
      const matchesText = !query || Object.values(row).some(value => String(value ?? '').toLowerCase().includes(query));
      return matchesText && (!group || row.supergroup === group);
    });
  };
  search.addEventListener('input',apply);
  groupFilter.addEventListener('change',apply);
}

function initContigTables() {
  if (!window.Tabulator) return;
  if (REPORT_DATA.abundanceRows.length) new Tabulator('#abundanceTable', {
    data:REPORT_DATA.abundanceRows, ...makeContigOptions(REPORT_DATA.abundanceRows.length),
    columns:[
      {title:L.contig,field:'contig_id',headerFilter:'input',minWidth:150},
      {title:'RPM',field:'rpm',sorter:'number',hozAlign:'right',formatter:cell=>fmtUnavailable(cell.getValue())},
      {title:L.mappedReads,field:'mapped_reads',sorter:'number',hozAlign:'right',formatter:cell=>fmtUnavailable(cell.getValue())},
      {title:L.meanDepth,field:'mean_depth',sorter:'number',hozAlign:'right',formatter:cell=>fmtUnavailable(cell.getValue())},
      {title:L.coverage,field:'coverage_breadth_pct',sorter:'number',hozAlign:'right',formatter:cell=>fmtUnavailable(cell.getValue())}
    ]
  });
}

function initStageTable() {
  if (!window.Tabulator || !REPORT_DATA.stageRows.length) return;
  new Tabulator('#stageTable', {
    data:REPORT_DATA.stageRows, layout:'fitDataStretch', pagination:false,
    columns:[
      {title:L.stage,field:'stage',minWidth:180},
      {title:L.status,field:'status',width:100},
      {title:L.input,field:'input',minWidth:160},
      {title:L.output,field:'output',minWidth:160},
      {title:L.time,field:'minutes',sorter:'number',hozAlign:'right',formatter:cell=>`${fmtUnavailable(cell.getValue())} ${L.minutes}`,width:110}
    ]
  });
}

function setupPrintDetails() {
  let savedDetails = null;
  window.addEventListener('beforeprint', () => {
    const details = document.querySelectorAll('details');
    savedDetails = Array.from(details).map(item => item.open);
    details.forEach(item => { item.open = true; });
  });
  window.addEventListener('afterprint', () => {
    if (!savedDetails) return;
    const details = document.querySelectorAll('details');
    Array.from(details).forEach((item, index) => {
      if (index < savedDetails.length) item.open = savedDetails[index];
    });
    savedDetails = null;
  });
}

window.addEventListener('DOMContentLoaded',()=>{initCharts();initSequenceTable();initContigTables();initStageTable();setupPrintDetails();});
"""
