"""从运行产物中提取各阶段耗时、序列漏斗和超群详情。"""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from cdc_pipeline.stages.common import read_tsv


def _parse_time(s: str) -> datetime | None:
    if not s:
        return None
    try:
        return datetime.fromisoformat(s)
    except (ValueError, TypeError):
        return None


def extract_stage_timings(run_dir: str | Path) -> list[dict]:
    """从 state.json 提取各阶段耗时（分钟）。"""
    state = json.loads((Path(run_dir) / "state.json").read_text())
    stages = state.get("stages", {})
    total_start = None
    total_end = None
    timings = []

    for name in [
        "preprocess", "orf", "prefilter", "detect", "representatives",
        "supergroup", "classifier", "phylo", "visualization",
        "abundance", "report",
    ]:
        stage = stages.get(name)
        if not stage:
            continue
        started = _parse_time(stage.get("started_at"))
        ended = _parse_time(stage.get("ended_at"))
        duration = 0.0
        if started and ended:
            duration = (ended - started).total_seconds() / 60.0
            if total_start is None or started < total_start:
                total_start = started
            if total_end is None or ended > total_end:
                total_end = ended
        timings.append({
            "stage": name,
            "status": stage.get("status", ""),
            "duration_minutes": round(duration, 2),
        })

    total = 0.0
    if total_start and total_end:
        total = (total_end - total_start).total_seconds() / 60.0

    return timings, round(total, 2)


def _find(run_dir: Path, pattern: str) -> Path | None:
    matches = sorted(run_dir.glob(pattern))
    return matches[0] if matches else None


def _find_tiered(
    run_dir: Path,
    filename: str,
    stage: str,
    *,
    delivery_dir: str = "tables",
    delivery_names: tuple[str, ...] | None = None,
) -> Path | None:
    """按阶段原始产物 → 新交付目录 → 旧 results 根目录查找文件。

    原始阶段产物优先，报告重建不依赖交付副本；同时兼容 tables/、sequences/
    和早期平铺布局，避免交付文件移动后统计归零。
    """
    for pattern in (f"step*_{stage}/{filename}", f"0*_{stage}/{filename}"):
        found = _find(run_dir, pattern)
        if found:
            return found
    names = delivery_names or (filename,)
    for name in names:
        found = _find(run_dir, f"results/{delivery_dir}/{name}")
        if found:
            return found
    for name in names:
        found = _find(run_dir, f"results/{name}")
        if found:
            return found
    return None


def _count_fasta(path: Path | None) -> int:
    return sum(1 for line in open(path) if line.startswith(">")) if path and path.is_file() else 0


def _count_tsv(path: Path | None) -> int:
    return len(read_tsv(path)) if path and path.is_file() else 0


def extract_supergroup_outcome(run_dir: Path | str) -> dict[str, int]:
    run_dir = Path(run_dir)
    """合并参考匹配与分类的分类结果，避免把未分类和已分类混为一谈。"""
    assignment_path = _find_tiered(run_dir, "supergroup_assignments.tsv", "supergroup")
    assignments = _count_tsv(assignment_path)
    assignment_rows = read_tsv(assignment_path) if assignment_path else []
    statuses = Counter(row.get("status", "") for row in assignment_rows)
    classifier_path = _find_tiered(run_dir, "classifier_assignments.tsv", "classifier")
    classifier_rows = read_tsv(classifier_path) if classifier_path and classifier_path.is_file() else []
    remaining_path = _find_tiered(run_dir, "remaining_unclassified.faa", "classifier")
    unresolved = _count_fasta(remaining_path) if remaining_path else statuses.get("unassigned", 0)
    return {
        "representatives": assignments,
        "reference_classified": assignments - statuses.get("unassigned", 0),
        "reference_unassigned": statuses.get("unassigned", 0),
        "classifier_input": len(classifier_rows),
        "classifier_classified": sum(1 for row in classifier_rows if row.get("supergroup")),
        "classifier_tree_available": sum(1 for row in classifier_rows if row.get("tree_status") == "available"),
        "tree_available": statuses.get("assigned", 0) + sum(1 for row in classifier_rows if row.get("tree_status") == "available"),
        "pending_merge": statuses.get("classified_pending_merge", 0) + sum(1 for row in classifier_rows if row.get("tree_status") == "pending_merge"),
        "no_tree": statuses.get("classified_no_tree", 0) + sum(1 for row in classifier_rows if row.get("tree_status") == "no_tree"),
        "unresolved": unresolved,
    }


def extract_sequence_funnel(run_dir: str | Path) -> list[dict]:
    """构建带明确序列类型的数量漏斗，避免核酸序列与蛋白序列混淆。"""
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    route = manifest.get("route", "")
    funnel = []

    def add(key: str, stage: str, count: int, unit: str, note: str = ""):
        funnel.append({"key": key, "stage": stage, "count": count, "unit": unit, "note": note})


    if route in {"contigs", "proteins"}:
        input_path = Path(manifest.get("inputs", {}).get("input", ""))
        unit = "条核酸 contigs" if route == "contigs" else "条蛋白序列"
        add("input", "原始输入", _count_fasta(input_path), unit,
            "组装后的核酸 contigs" if route == "contigs" else "用户提供的蛋白序列")
    else:
        add("input", "原始输入", 0, "条 reads", "原始测序 reads 经质控和组装后转为核酸 contigs")


    proteins = _find(run_dir, "step*_orf/proteins.faa") or _find(run_dir, "0*_orf/proteins.faa")
    orf_count = _count_fasta(proteins)
    if orf_count:
        add("orf", "fastorf ORF 预测", orf_count, "条蛋白 ORFs", "核酸 contigs 翻译得到的 ORF 蛋白，≥75 nt")


    main_faa = _find(run_dir, "step*_prefilter/main_analysis.faa") or _find(run_dir, "0*_prefilter/main_analysis.faa")
    main_count = _count_fasta(main_faa)
    if main_count:
        add("length", "长度过滤（≥150 aa）", main_count, "条蛋白 ORFs", "进入粗筛和模型分析")


    model_faa = _find(run_dir, "step*_prefilter/model_input.faa") or _find(run_dir, "0*_prefilter/model_input.faa")
    model_count = _count_fasta(model_faa)
    if model_count:
        kv = _count_fasta(_find(run_dir, "step*_prefilter/KEEP_viral.faa") or _find(run_dir, "0*_prefilter/KEEP_viral.faa"))
        ku = _count_fasta(_find(run_dir, "step*_prefilter/KEEP_unknown.faa") or _find(run_dir, "0*_prefilter/KEEP_unknown.faa"))
        add("prefilter", "RefSeq 粗筛后送模型", model_count, "条蛋白 ORFs", f"病毒命中 KEEP_viral={kv}，未知保守 KEEP_unknown={ku}")


    pred_path = _find_tiered(run_dir, "orf_predictions.tsv", "detect")
    if pred_path:
        rows = [r for r in read_tsv(pred_path) if r.get("analysis_set") == "main"]
        positive = sum(1 for r in rows if r.get("call") == "positive")
        add("positive", "RdRP 阳性蛋白", positive, "条蛋白 ORFs", "RdRP 识别模型判定")


    reps = _find_tiered(
        run_dir, "representatives.faa", "representatives",
        delivery_dir="sequences", delivery_names=("representative_orfs.faa", "representatives.faa"),
    )
    rep_count = _count_fasta(reps)
    if rep_count:
        add("representatives", "去冗余代表序列", rep_count, "条代表蛋白", "CD-HIT 90% identity / 80% coverage 后")


    outcome = extract_supergroup_outcome(run_dir)
    if outcome["representatives"]:
        add("supergroup", "超群分类", outcome["representatives"], "条代表蛋白",
            f"参考匹配={outcome['reference_classified']}，分类={outcome['classifier_classified']}；"
            f"骨架树可用={outcome['tree_available']}，骨架树正在建设中={outcome['pending_merge'] + outcome['no_tree']}")


    placements = _find_tiered(run_dir, "placements.tsv", "phylo")
    placed = _count_tsv(placements)
    if placed:
        add("phylo", "系统发育放置", placed, "条代表蛋白", "有可用骨架树的超群")

    return funnel


def extract_stage_metrics(run_dir: str | Path, timings: list[dict] | None = None) -> list[dict]:
    """按业务执行顺序返回每步输入/输出序列数量、单位和耗时。"""
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    route = manifest.get("route", "")
    state = json.loads((run_dir / "state.json").read_text())
    stages = state.get("stages", {})
    timing_map = {row["stage"]: row for row in (timings or extract_stage_timings(run_dir)[0])}
    funnel_rows = extract_sequence_funnel(run_dir)
    by_key = {row["key"]: row for row in funnel_rows}
    outcome = extract_supergroup_outcome(run_dir)

    def get(key, default_unit=None):
        row = by_key.get(key)
        if row is None:
            return 0, default_unit
        return row.get("count", 0), row.get("unit") or default_unit

    rows = []

    def add(key, label, input_count, input_unit, output_count, output_unit, note=""):
        stage = stages.get(key, {})
        timing = timing_map.get(key, {})
        status = stage.get("status", "not_run")
        if key == "report" and status == "running":
            status = "done"
        rows.append({
            "stage_key": key,
            "stage": label,
            "status": status,
            "input_count": input_count,
            "input_unit": input_unit,
            "output_count": output_count,
            "output_unit": output_unit,
            "duration_minutes": timing.get("duration_minutes", 0.0),
            "note": note,
        })

    if route == "reads":
        preprocess = stages.get("preprocess", {}).get("summary", {}) or {}
        qc = preprocess.get("qc", {})
        after = qc.get("after_filtering", {})
        add("preprocess", "测序数据质控与组装", 0, "条 reads",
            0, "条核酸 contigs", f"质控后 reads={after.get('total_reads', '')}")

    input_count, input_unit = get("input", "条核酸 contigs" if route == "contigs" else "条蛋白序列")
    orf_count, orf_unit = get("orf", "条蛋白 ORFs")
    main_count, _ = get("length", "条蛋白 ORFs")
    model_count, model_unit = get("prefilter", "条蛋白 ORFs")
    positive_count, positive_unit = get("positive")
    rep_count, rep_unit = get("representatives", "条代表蛋白")
    placed_count, placed_unit = get("phylo", "条代表蛋白")

    add("orf", "蛋白 ORF 预测", input_count, input_unit, orf_count, orf_unit, "fastorf，核酸 contigs → 蛋白 ORFs")
    add("prefilter", "宿主/细胞同源粗筛", orf_count, orf_unit, model_count, model_unit,
        f"长度≥150 aa 的蛋白 ORFs={main_count}；RefSeq 粗筛后进入模型")
    add("detect", "RdRP 识别", model_count, model_unit, positive_count, positive_unit or "条蛋白 ORFs",
        "RdRP 识别模型给出阳性概率并按阈值判定")
    add("representatives", "阳性序列去冗余", positive_count, positive_unit or "条蛋白 ORFs", rep_count, rep_unit,
        "CD-HIT 聚类并选择代表序列")
    add("supergroup", "超群分类", rep_count, rep_unit, outcome["reference_classified"], rep_unit,
        f"参考匹配未覆盖 {outcome['reference_unassigned']} 条进入分类")
    if outcome["reference_unassigned"]:
        add("classifier", "分类", outcome["reference_unassigned"], rep_unit,
            outcome["classifier_classified"], rep_unit, "仅处理前述未匹配序列")
    add("phylo", "系统发育放置（LucaPhylo）", outcome["tree_available"], rep_unit, placed_count, placed_unit,
        f"骨架树正在建设中：{outcome['pending_merge'] + outcome['no_tree']} 条暂未放置")
    add("report", "结果汇总与可视化", rep_count, rep_unit, 1, "份报告", "生成双语主报告和离线系统发育浏览器")
    return rows


def extract_supergroup_details(run_dir: str | Path) -> list[dict]:
    """提取超群分配详情，含参考匹配与分类结果。"""
    run_dir = Path(run_dir)
    sg_tsv = _find_tiered(run_dir, "supergroup_assignments.tsv", "supergroup")
    if not sg_tsv or not sg_tsv.is_file():
        return []

    rows = read_tsv(sg_tsv)
    details = []
    for r in rows:
        detail = {
            "query_id": r.get("query_id", ""),
            "supergroup": r.get("supergroup", ""),
            "status": r.get("status", ""),
            "method": r.get("classification_method", ""),
            "reference_id": r.get("reference_id", ""),
            "evalue": r.get("evalue", ""),
            "pident": r.get("pident", ""),
            "qcovhsp": r.get("qcovhsp", ""),
        }
        details.append(detail)


    classifier_tsv = _find_tiered(run_dir, "classifier_assignments.tsv", "classifier")
    if classifier_tsv and classifier_tsv.is_file():
        classifier_rows = read_tsv(classifier_tsv)
        classifier_map = {r["query_id"]: r for r in classifier_rows}
        for detail in details:
            if detail["query_id"] in classifier_map:
                row = classifier_map[detail["query_id"]]
                detail["method"] = "classifier"
                detail["supergroup"] = row.get("supergroup", "")
                detail["confidence"] = row.get("confidence", "")

    return details


def extract_placement_details(run_dir: str | Path) -> list[dict]:
    """从 placements.tsv 提取系统发育放置信息。"""
    run_dir = Path(run_dir)
    pl_tsv = _find_tiered(run_dir, "placements.tsv", "phylo")
    if not pl_tsv or not pl_tsv.is_file():
        return []

    rows = read_tsv(pl_tsv)
    details = []
    for r in rows:
        details.append({
            "query_id": r.get("query_id", ""),
            "supergroup": r.get("supergroup", ""),
            "nearest_neighbor": r.get("nn1", ""),
            "nn1_distance": r.get("nn1_dist", ""),
        })
    return details


def build_comprehensive_metrics(run_dir: str | Path) -> dict:
    """汇总所有指标。"""
    timings, total_minutes = extract_stage_timings(run_dir)
    funnel = extract_sequence_funnel(run_dir)
    supergroup_details = extract_supergroup_details(run_dir)
    placement_details = extract_placement_details(run_dir)


    from collections import Counter
    sg_counter = Counter(d["supergroup"] for d in supergroup_details if d["supergroup"])
    method_counter = Counter(d.get("method", "") for d in supergroup_details)
    status_counter = Counter(d["status"] for d in supergroup_details)

    return {
        "timings": timings,
        "total_minutes": total_minutes,
        "funnel": funnel,
        "supergroup_summary": {
            "by_supergroup": dict(sg_counter.most_common()),
            "by_method": dict(method_counter),
            "by_status": dict(status_counter),
        },
        "supergroup_details": supergroup_details,
        "placement_details": placement_details,
        "n_supergroup_details": len(supergroup_details),
        "n_placement_details": len(placement_details),
    }
