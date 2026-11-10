#!/usr/bin/env python3
"""校验 phylogenetic_explorer.html（树浏览器）的关键内容与清理项。

用法：
    python scripts/verify_tree_explorer.py runs/real_data/*/results/phylogenetic_explorer.html
"""

from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

FORBIDDEN = [
    "LWR", "confidence", "pendant", "VirID 分类", "骨架树库", "🧬",
    "HOST_CLASSES", "hostClassForTip", "宿主大类",

    "host_suppression_reason",
]
REQUIRED = ["exportDlg", "globalNotice", "original_id", "svg2pdf", "jsPDF", "report-meta"]


def check_file(path: Path, require_original_ids: bool = True) -> list[str]:
    """校验树浏览器契约。

    require_original_ids 默认 True：contigs/reads 路线有 orf_contig_map，
    每条 placement 都应带回原始 ORF ID。proteins 路线没有映射表，
    original_id 允许为空，但模板残留与 placement schema 契约仍执行。
    """
    problems: list[str] = []
    text = path.read_text(encoding="utf-8")

    for token in FORBIDDEN:
        if token in text:
            problems.append(f"残留禁止内容: {token}")
    for token in REQUIRED:
        if token not in text:
            problems.append(f"缺少必需组件: {token}")

    match = re.search(r'<script id="virid-bundles" type="application/json">\s*(.*?)\s*</script>', text, re.S)
    if not match:
        problems.append("未找到 virid-bundles 数据")
        return problems
    try:
        data = json.loads(match.group(1))
    except json.JSONDecodeError as exc:
        problems.append(f"bundles JSON 解析失败: {exc}")
        return problems

    if not data:

        meta = re.search(r'<script id="report-meta" type="application/json">\s*(.*?)\s*</script>', text, re.S)
        try:
            notice = json.loads(meta.group(1)).get("notice", "") if meta else ""
        except json.JSONDecodeError:
            notice = ""
        if not notice:
            problems.append("bundles 为空且无说明提示")
        return problems
    missing_original = 0
    legacy_fields = {"host", "host_suppression_reason"}
    for supergroup, bundle in data.items():
        for placement in bundle.get("placements", []):
            if "original_id" not in placement:
                problems.append(f"{supergroup} 的 placement 缺 original_id 字段")
                break
            if not placement["original_id"]:
                missing_original += 1
            if placement.get("neighbors"):
                neighbor = placement["neighbors"][0]
                if not isinstance(neighbor, list) or len(neighbor) != 2:
                    problems.append(f"{supergroup} 的邻居不是 [id, dist] 两元组")
                    break
            leaked = sorted(legacy_fields & set(placement))
            if leaked:
                problems.append(
                    f"{supergroup} 的 {placement.get('query_id')} 含旧版宿主概率字段 {leaked}；"
                    "树为旧产物，请重跑 visualization"
                )
    if require_original_ids and missing_original:
        problems.append(f"{missing_original} 条 placement 的 original_id 为空")

    meta = re.search(r'<script id="report-meta" type="application/json">\s*(.*?)\s*</script>', text, re.S)
    if meta:
        try:
            json.loads(meta.group(1))
        except json.JSONDecodeError:
            problems.append("report-meta JSON 解析失败")

    return problems


def main(argv: list[str]) -> int:
    paths = [Path(item) for item in argv[1:]] or [
        Path(item) for item in sorted(glob.glob("runs/real_data/**/results/phylogenetic_explorer.html", recursive=True))
    ]
    if not paths:
        print("未找到 phylogenetic_explorer.html", file=sys.stderr)
        return 2
    failures = 0
    for path in paths:
        problems = check_file(path)
        sample = path.parent.parent.name
        if problems:
            failures += 1
            print(f"FAIL {sample}/{path.name}")
            for problem in problems[:6]:
                print(f"     - {problem}")
        else:
            print(f"OK   {sample}/{path.name}")
    print(f"\n{len(paths) - failures}/{len(paths)} 通过")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
