"""真实 query bundle 和当前运行专属离线系统发育树。"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from ete3 import Tree

from cdc_pipeline.assets import load_supergroup_map
from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.placement import Neighbor, Placement
from cdc_pipeline.risk_assessment import parse_taxonomy_flag
from cdc_pipeline.stages.common import artifact_path, read_tsv
from cdc_pipeline.state import atomic_write_json

RANK_ORDER = ["Realm", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Suborder", "Family", "Subfamily", "Genus", "Species"]


def _taxonomy(path: Path) -> tuple[dict, dict]:
    taxonomy, hosts = {}, {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            tip = row["Cluster_sseqid"]
            taxonomy[tip] = {rank: row.get(rank, "").strip() for rank in RANK_ORDER if row.get(rank, "").strip()}
            marks = [int(parse_taxonomy_flag(row.get(key)) is True) for key in ["Human", "Vertebrates", "Plant"]]
            if any(marks):
                hosts[tip] = marks
    return taxonomy, hosts


def build_bundle(
    supergroup: str,
    tree_path: str | Path,
    taxonomy_path: str | Path,
    placements: list[Placement],
    provenance: dict,
    original_ids: dict[str, str] | None = None,
) -> dict:
    tree_path = Path(tree_path)
    tree = Tree(str(tree_path), format=1)
    tips = {tip.name for tip in tree}
    taxonomy, host_marks = _taxonomy(Path(taxonomy_path))
    missing = tips - set(taxonomy)
    if missing:
        raise ValueError(f"{supergroup} 骨架叶缺少分类信息: {sorted(missing)[:5]}")
    ranks = [rank for rank in RANK_ORDER if any(taxonomy[tip].get(rank) for tip in tips)]
    value_counts = {rank: {} for rank in ranks}
    for tip in tips:
        for rank in ranks:
            value = taxonomy[tip].get(rank)
            if value:
                value_counts[rank][value] = value_counts[rank].get(value, 0) + 1
    placement_rows = []
    for placement in placements:
        neighbors = []
        for neighbor in placement.neighbors:
            neighbors.append([neighbor.tip_id, round(neighbor.distance, 6)])
        placement_rows.append({
            "query_id": placement.query_id,

            "original_id": (original_ids or {}).get(placement.query_id, ""),
            "supergroup": placement.supergroup,
            "neighbors": neighbors,
        })
    return {
        "supergroup": supergroup,
        "kind": "virid",
        "backbone_newick": tree_path.read_text().strip(),
        "n_backbone": len(tips),
        "n_query": len(placement_rows),
        "placements": placement_rows,
        "truth": None,
        "metrics": {},
        "taxonomy": taxonomy,
        "host_marks": host_marks,
        "ranks": ranks,
        "value_counts": value_counts,
        "provenance": provenance,
    }


def build_tree_html(bundles: list[dict], output_path: str | Path, notice: str = "") -> Path:
    source = PROJECT_ROOT / "third_party" / "placement_view"
    template = (source / "tree_explorer_template.html").read_text()
    core = (source / "phylo_core.js").read_text()
    core = re.sub(r"^/\*.*?\*/", "", core, count=1, flags=re.S).strip()

    def safe(value: str) -> str:
        return value.replace("</script", "<\\/script")

    virids = {bundle["supergroup"]: bundle for bundle in bundles}
    vendor = source / "vendor"
    jspdf = (vendor / "jspdf.umd.min.js").read_text() if (vendor / "jspdf.umd.min.js").is_file() else ""
    svg2pdf = (vendor / "svg2pdf.umd.min.js").read_text() if (vendor / "svg2pdf.umd.min.js").is_file() else ""
    output = template.replace("/*__PHYLO_CORE__*/", safe(core))
    output = output.replace("/*__BACKBONE_BUNDLES__*/", "{}")
    output = output.replace("/*__FITCH_BUNDLES__*/", "{}")
    output = output.replace("/*__VIRID_BUNDLES__*/", safe(json.dumps(virids, ensure_ascii=False)))
    output = output.replace("/*__REPORT_META__*/", safe(json.dumps({"notice": notice}, ensure_ascii=False)))
    output = output.replace("/*__JSPDF__*/", safe(jspdf))
    output = output.replace("/*__SVG2PDF__*/", safe(svg2pdf))
    markers = ["__PHYLO_CORE__", "__BACKBONE_BUNDLES__", "__FITCH_BUNDLES__", "__VIRID_BUNDLES__", "__REPORT_META__", "__JSPDF__", "__SVG2PDF__"]
    if any(marker in output for marker in markers):
        raise ValueError("离线树模板仍有未替换占位符")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(output, encoding="utf-8")
    return output_path


def _placements(path: str | Path, top_k: int) -> list[Placement]:
    placements = []
    for row in read_tsv(path):
        neighbors = []
        for index in range(1, top_k + 1):
            if row.get(f"nn{index}"):
                neighbors.append(Neighbor(row[f"nn{index}"], float(row[f"nn{index}_dist"])))
        placements.append(Placement(row["query_id"], row["supergroup"], neighbors))
    return placements


def _original_ids(path: str | Path | None, wanted: set[str] | None = None) -> dict[str, str]:
    """query_id → 原始 fastorf ORF ID。缺失时返回空映射（旧运行兼容）。

    orf_contig_map 在大样本上可达数百 MB 到 GB，按行流式过滤，只保留已放置的 ID。
    """
    if not path or not Path(path).is_file():
        return {}
    mapping: dict[str, str] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query_id = row.get("query_id")
            if not query_id or (wanted is not None and query_id not in wanted):
                continue
            original = row.get("orf_id") or row.get("original_orf_id") or ""
            if original:
                mapping[query_id] = original
    return mapping


def _iter_tsv(path: str | Path):
    """逐行读取 TSV，避免把大 orf_contig_map 整体载入内存。"""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def run(context: PipelineContext) -> StageResult:
    from cdc_pipeline.supergroup_mapping import lucaphylo_to_asset_name
    placements = _placements(context.data["placements_tsv"], context.config["phylo"]["top_k"])
    original_ids = _original_ids(context.data.get("orf_contig_map"), {p.query_id for p in placements})
    assets = load_supergroup_map(context.config)
    bundles = []
    bundle_paths = {}
    for sg, matrix in context.data.get("phylo_matrices", {}).items():
        current = [placement for placement in placements if placement.supergroup == sg]
        if not current:
            continue
        asset_name = lucaphylo_to_asset_name(sg)
        if asset_name not in assets:
            continue
        asset = assets[asset_name]
        bundle = build_bundle(sg, asset.tree, asset.taxonomy, current, {
            "run_id": context.run_dir.name,
            "matrix": matrix,
            "tree": str(asset.tree),
            "adapter": str(asset.adapter),
            "query_source": "real_pipeline_query",
        }, original_ids=original_ids)
        path = context.stage_dir / "bundles" / f"{sg.replace('/', '_')}.json"
        atomic_write_json(path, bundle)
        bundle_paths[sg] = path
        bundles.append(bundle)
    no_tree_groups = sorted(
        set(context.data.get("supergroup_no_tree_groups", {}))
        | set(context.data.get("classifier_groups_no_tree", {}))
    )
    no_tree_count = int(context.data.get("supergroup_counts", {}).get("classified_no_tree", 0)) + int(
        context.data.get("classifier_counts", {}).get("classified_no_tree", 0)
    )
    notice = ""
    if not bundles:
        notice = "本样本没有完成系统发育放置的序列（RdRP 阴性或无可放置序列）。具体检测过程与阴性证据见主报告。"
    elif no_tree_groups or no_tree_count:
        names = "、".join(no_tree_groups[:8]) + (" 等" if len(no_tree_groups) > 8 else "")
        notice = (
            f"部分识别序列的骨架树正在建设中：{names or '若干超群'}（共 {no_tree_count} 条）。"
            "这些序列已完成超群分类，待骨架树建设完成后纳入系统发育放置，详见主报告。"
        )
    output = build_tree_html(bundles, context.stage_dir / "tree_explorer.html", notice=notice)


    result_dir = context.run_dir / "results"
    tree_dir = result_dir / "trees"
    tree_dir.mkdir(parents=True, exist_ok=True)
    tree_files: dict[str, Path] = {}
    annotation_rows = [[
        "query_id", "original_id", "supergroup",
        "nn1", "nn1_dist", "nn2", "nn2_dist", "nn3", "nn3_dist", "nn4", "nn4_dist", "nn5", "nn5_dist",
    ]]
    for bundle in bundles:
        safe_name = re.sub(r"[^\w.-]+", "_", bundle["supergroup"])
        tree_file = tree_dir / f"{safe_name}.nwk"
        tree_file.write_text(bundle["backbone_newick"].strip() + "\n", encoding="utf-8")
        tree_files[bundle["supergroup"]] = tree_file
        for placement in bundle["placements"]:
            row = [placement["query_id"], placement.get("original_id", ""), placement.get("supergroup", bundle["supergroup"])]
            for index in range(5):
                neighbor = placement["neighbors"][index] if index < len(placement["neighbors"]) else None
                row += [neighbor[0], f"{neighbor[1]:.6f}"] if neighbor else ["", ""]
            annotation_rows.append(row)
    annotation_path = result_dir / "tables" / "placement_annotations.tsv"
    annotation_path.parent.mkdir(parents=True, exist_ok=True)
    annotation_path.write_text("\n".join("\t".join(row) for row in annotation_rows) + "\n", encoding="utf-8")

    return StageResult(
        artifacts={
            "tree_html": artifact_path(context, output),
            "bundles": {sg: artifact_path(context, path) for sg, path in bundle_paths.items()},
            "trees": {sg: artifact_path(context, path) for sg, path in tree_files.items()},
            "placement_annotations": artifact_path(context, annotation_path),
        },
        summary={"bundle_count": len(bundles), "query_count": len(placements), "reason": "" if bundles else "no_placements"},
        data={
            "tree_html": str(output),
            "bundle_paths": {sg: str(path) for sg, path in bundle_paths.items()},
            "tree_files": {sg: str(path) for sg, path in tree_files.items()},
            "placement_annotations": str(annotation_path),
        },
        status="done" if bundles else "skipped",
    )
