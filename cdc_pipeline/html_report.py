"""独立 HTML 报告的数据连接、统计和自包含页面生成。"""

from __future__ import annotations

import csv
import math
import statistics
from collections import Counter
from pathlib import Path

from cdc_pipeline.stages.common import read_tsv


def _rows(data: dict, key: str) -> list[dict[str, str]]:
    path = data.get(key)
    return read_tsv(path) if path and Path(path).is_file() else []


def _iter_tsv(path: str | Path):
    """逐行读取 TSV，调用方可在遍历时过滤，避免大表整体入内存。"""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _filtered_rows(data: dict, key: str, required_ids: set[str], id_key: str) -> list[dict]:
    path = data.get(key)
    if not path or not required_ids or not Path(path).is_file():
        return []
    return [row for row in _iter_tsv(path) if row.get(id_key) in required_ids]


def _by(rows: list[dict], key: str) -> dict[str, dict]:
    return {row[key]: row for row in rows if row.get(key)}


def _number(value):
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return ""
    return number if math.isfinite(number) else ""


def _integer(value):
    if value is None or value == "":
        return ""
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return ""


def _finite_sum(values: list[float]) -> float | None:
    if not values:
        return None
    try:
        result = math.fsum(values)
    except (OverflowError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _text(value) -> str:
    return str(value).strip() if value is not None else ""


def supergroup_shares(rows: list[dict]) -> list[dict]:
    if not rows:
        return []
    counts = Counter(row.get("supergroup") or "Unassigned" for row in rows)
    ordered = sorted(((sg, count) for sg, count in counts.items() if sg != "Unassigned"), key=lambda item: (-item[1], item[0]))
    if "Unassigned" in counts:
        ordered.append(("Unassigned", counts["Unassigned"]))
    total = len(rows)
    result = []
    remaining = 100.0
    for index, (supergroup, count) in enumerate(ordered):
        percent = round(remaining, 2) if index == len(ordered) - 1 else round(count / total * 100, 2)
        remaining = round(remaining - percent, 2)
        result.append({"supergroup": supergroup, "count": count, "percent": percent})
    return result


def _distribution(values: list[float]) -> dict:
    values = sorted(values)
    if not values:
        return {"count": 0, "min": None, "median": None, "max": None, "bins": []}
    bins_count = min(10, max(1, math.ceil(math.sqrt(len(values)))))
    minimum, maximum = values[0], values[-1]
    if minimum == maximum:
        bins = [len(values)]
        edges = [(minimum, maximum)]
    else:
        bins = [0] * bins_count
        width = (maximum - minimum) / bins_count
        for value in values:
            index = min(int((value - minimum) / width), bins_count - 1)
            bins[index] += 1
        edges = [(minimum + index * width, minimum + (index + 1) * width) for index in range(bins_count)]
    return {
        "count": len(values),
        "min": minimum,
        "median": statistics.median(values),
        "max": maximum,
        "bins": [
            {"start": round(start, 4), "end": round(end, 4), "count": count}
            for (start, end), count in zip(edges, bins)
        ],
    }


def length_statistics(rows: list[dict]) -> dict:
    values = [_integer(row.get("aa_len")) for row in rows]
    return _distribution([value for value in values if value != ""])


def categorical_counts(rows: list[dict], key: str, empty_label: str = "Unassigned") -> list[dict]:
    counts = Counter(row.get(key) or empty_label for row in rows)
    return [{"name": name, "count": count} for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))]


def motif_statistics(rows: list[dict]) -> list[dict]:
    return [
        {"name": motif, "count": sum(motif in (row.get("motif_classes") or "").split(";") for row in rows)}
        for motif in ["F", "A", "B", "C"]
    ]


def confidence_statistics(rows: list[dict]) -> dict:
    """已废弃：不再计算 confidence（非 LucaPhylo 原生输出）。保留空实现以兼容。"""
    return {"count": 0, "min": None, "median": None, "max": None, "bins": []}


def _classifier_rows(context) -> list[dict]:
    data = context.data
    value = data.get("classifier_assignments")
    if isinstance(value, (list, tuple)):
        return list(value)
    path = Path(value) if value else None
    if not path or not path.is_file():
        run_dir = Path(getattr(context, "run_dir", "")) if getattr(context, "run_dir", None) else None
        candidates = sorted(run_dir.glob("step*_classifier/classifier_assignments.tsv")) if run_dir else []
        path = candidates[0] if candidates else None
    return read_tsv(path) if path and path.is_file() else []


def _load_report_data(context) -> dict:
    """加载一次报告数据；大 ORF 表仅保留所需 ID，并缓存到 context。"""
    cached = getattr(context, "_html_report_data_cache", None)
    if cached is not None:
        return cached

    data = context.data
    representative_rows = _rows(data, "representative_map")
    member_ids = {row.get("query_id", "") for row in representative_rows if row.get("query_id")}
    representative_ids = {
        row["query_id"] for row in representative_rows
        if str(row.get("is_representative")) == "1" and row.get("query_id")
    }
    detection = _by(_filtered_rows(data, "orf_predictions", representative_ids, "query_id"), "query_id")
    orf_map = _by(_filtered_rows(data, "orf_contig_map", member_ids, "query_id"), "query_id")
    positive_contigs = {row["contig_id"] for row in orf_map.values() if row.get("contig_id")}
    abundance = _by(_filtered_rows(data, "abundance_tsv", positive_contigs, "contig_id"), "contig_id")
    loaded = {
        "representative_rows": representative_rows,
        "detection": detection,
        "assignments": _by(_rows(data, "supergroup_assignments"), "query_id"),
        "placements": _by(_rows(data, "placements_tsv"), "query_id"),
        "classifier": _by(_classifier_rows(context), "query_id"),
        "orf_map": orf_map,
        "positive_contigs": positive_contigs,
        "abundance": abundance,
    }
    try:
        setattr(context, "_html_report_data_cache", loaded)
    except (AttributeError, TypeError):
        pass
    return loaded


def join_sequence_rows(context) -> list[dict]:
    loaded = _load_report_data(context)
    representatives = [
        row for row in loaded["representative_rows"]
        if str(row.get("is_representative")) == "1"
    ]
    detection = loaded["detection"]
    assignments = loaded["assignments"]
    placements = loaded["placements"]
    orf_map = loaded["orf_map"]
    abundance = loaded["abundance"]
    classifier = loaded["classifier"]
    output = []
    seen = set()
    for representative in representatives:
        query_id = representative["query_id"]
        if query_id in seen:
            raise ValueError(f"代表序列 ID 重复: {query_id}")
        seen.add(query_id)
        det = detection.get(query_id, {})
        assignment = assignments.get(query_id, {})
        placement = placements.get(query_id, {})
        classifier_row = classifier.get(query_id, {})
        orf_row = orf_map.get(query_id, {})
        contig_id = orf_row.get("contig_id", "")
        abundance_row = abundance.get(contig_id, {})
        assignment_status = _text(assignment.get("status"))
        supergroup = _text(assignment.get("supergroup"))
        method = _text(assignment.get("classification_method"))
        if supergroup:
            method = method or "reference"
        elif _text(classifier_row.get("supergroup")):
            supergroup = _text(classifier_row.get("supergroup"))
            method = "classifier"
            assignment_status = f"classified_{_text(classifier_row.get('tree_status')) or 'unknown'}"
        output.append({
            "query_id": query_id,
            "original_orf_id": orf_row.get("orf_id", ""),
            "contig_id": contig_id,
            "aa_len": _integer(representative.get("aa_len")),
            "call": det.get("call", ""),
            "rdrp_prob": _number(det.get("rdrp_prob")),
            "supergroup": supergroup,
            "classification_method": method,
            "classification_status": assignment_status,
            "reference_id": assignment.get("reference_id", ""),
            "evalue": _number(assignment.get("evalue")),
            "bitscore": _number(assignment.get("bitscore")),
            "pident": _number(assignment.get("pident")),
            "nn1": placement.get("nn1", ""),
            "nn1_dist": _number(placement.get("nn1_dist")),
            "rpm": _number(abundance_row.get("RPM")),
            "mean_depth": _number(abundance_row.get("mean_depth")),
            "coverage_breadth_pct": _number(abundance_row.get("coverage_breadth_pct")),
        })
    return output


def _weighted_mean(rows: list[dict], value_key: str) -> float | None:
    weighted = [
        (row[value_key], row["length"]) for row in rows
        if row[value_key] != "" and row["length"] != "" and row["length"] > 0
    ]
    length = sum(length for _, length in weighted)
    numerator = _finite_sum([value * length for value, length in weighted])
    if not length or numerator is None:
        return None
    result = numerator / length
    return result if math.isfinite(result) else None


def _abundance_summary(contigs: set[str], ambiguous: set[str], abundance_rows: dict[str, dict]) -> dict:
    rows = [abundance_rows[contig] for contig in sorted(contigs - ambiguous) if contig in abundance_rows]
    rpm_values = [row["rpm"] for row in rows if row["rpm"] != ""]
    mapped_values = [row["mapped_reads"] for row in rows if row["mapped_reads"] != ""]
    return {
        "n_contigs": len(rows),
        "mapped_reads": sum(mapped_values) if mapped_values else None,
        "RPM": _finite_sum(rpm_values),
        "mean_depth_weighted": _weighted_mean(rows, "mean_depth"),
        "coverage_breadth_pct": _weighted_mean(rows, "coverage_breadth_pct"),
    }


def enrich_report_summary(context, rows: list[dict], summary: dict) -> dict:
    """原地补齐全阳性成员口径的超群和丰度汇总。"""
    loaded = _load_report_data(context)
    representative_rows = loaded["representative_rows"]
    orf_map = loaded["orf_map"]

    positive_contigs = loaded["positive_contigs"]

    abundance_rows = []
    for contig in sorted(loaded["positive_contigs"]):
        source = loaded["abundance"].get(contig)
        if source is None:
            continue
        abundance_rows.append({
            "contig_id": contig,
            "length": _integer(source.get("length")),
            "mapped_reads": _integer(source.get("mapped_reads")),
            "rpm": _number(source.get("RPM")),
            "mean_depth": _number(source.get("mean_depth")),
            "coverage_breadth_pct": _number(source.get("coverage_breadth_pct")),
            "total_clean_reads": _integer(source.get("total_clean_reads")),
        })
    abundance_by_contig = {row["contig_id"]: row for row in abundance_rows}

    groups: dict[str, dict] = {}
    representative_groups = {}
    for row in rows:
        supergroup = row.get("supergroup") or "Unassigned"
        representative_groups[row["query_id"]] = supergroup
        group = groups.setdefault(supergroup, {"rows": [], "orf_ids": set(), "contigs": set()})
        group["rows"].append(row)
    for member in representative_rows:
        query_id = member.get("query_id", "")
        representative = member.get("representative_id", "")
        supergroup = representative_groups.get(representative)
        if supergroup is None:
            continue
        group = groups[supergroup]
        if query_id:
            group["orf_ids"].add(query_id)
        orf_row = orf_map.get(query_id, {})
        if orf_row.get("contig_id"):
            group["contigs"].add(orf_row["contig_id"])

    contig_groups: dict[str, set[str]] = {}
    for supergroup, group in groups.items():
        if supergroup == "Unassigned":
            continue
        for contig in group["contigs"]:
            contig_groups.setdefault(contig, set()).add(supergroup)
    multi_supergroup_contigs = [
        {"contig_id": contig, "supergroups": sorted(supergroups)}
        for contig, supergroups in sorted(contig_groups.items()) if len(supergroups) > 1
    ]
    ambiguous = {row["contig_id"] for row in multi_supergroup_contigs}

    supergroup_stats = []
    for supergroup in sorted(groups, key=lambda name: (name == "Unassigned", name)):
        group = groups[supergroup]
        methods = Counter(row.get("classification_method") or "Unassigned" for row in group["rows"])
        statuses = Counter(row.get("classification_status") or "unassigned" for row in group["rows"])
        supergroup_stats.append({
            "supergroup": supergroup,
            "representatives": len(group["rows"]),
            "positive_orfs": len(group["orf_ids"]),
            "positive_contigs": len(group["contigs"]),
            "representative_ids": sorted(row["query_id"] for row in group["rows"]),
            "positive_orf_ids": sorted(group["orf_ids"]),
            "contig_ids": sorted(group["contigs"]),
            "length_statistics": length_statistics(group["rows"]),
            "classification_methods": dict(sorted(methods.items())),
            "tree_statuses": dict(sorted(statuses.items())),
            "abundance": _abundance_summary(group["contigs"], ambiguous, abundance_by_contig),
        })

    summary["supergroup_stats"] = supergroup_stats
    summary["multi_supergroup_contigs"] = multi_supergroup_contigs
    summary["multi_supergroup_contig_count"] = len(multi_supergroup_contigs)
    summary["abundance_rows"] = abundance_rows
    summary["_abundance_rows"] = abundance_rows
    return summary


def build_html_report(context, output_path: str | Path, summary: dict | None = None, lang: str = "zh") -> Path:
    """将当前运行的代表序列和统计写成独立自包含 HTML。lang='zh'|'en'"""
    from cdc_pipeline.html_template_v2 import render_html

    rows = join_sequence_rows(context)
    if summary is None:
        summary = {
            "run_id": context.run_dir.name,
            "route": context.route,
            "sample": context.sample,
            "counts": {},
            "warnings": [],
            "links": {},
            "stages": context.state.snapshot.get("stages", {}),
        }
    enrich_report_summary(context, rows, summary)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_html(context, rows, summary, lang=lang), encoding="utf-8")
    return output_path


def build_bilingual_reports(context, output_dir: str | Path, summary: dict | None = None) -> tuple[Path, Path]:
    """生成中文和英文两版 main_results_report。"""
    output_dir = Path(output_dir)
    zh = build_html_report(context, output_dir / "main_results_report_zh.html", summary, lang="zh")
    en = build_html_report(context, output_dir / "main_results_report_en.html", summary, lang="en")
    return zh, en
