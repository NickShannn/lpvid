"""按超群运行 LucaPhylo stage-4 并汇总 placements.tsv。"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.assets import load_supergroup_map
from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.fasta import FastaRecord, read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.placement import Placement, placements_from_matrix, write_placements_tsv
from cdc_pipeline.process import pick_free_gpu, run_command
from cdc_pipeline.stages.common import artifact_path
from cdc_pipeline.state import atomic_write_json, compute_signature


@dataclass
class PhyloResult:
    placements_tsv: Path
    matrices: dict[str, Path]
    counts: dict[str, int]
    errors: dict[str, str]
    status: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def phylo_cache_key(base_model: str | Path, adapter: str | Path) -> str:
    base_model = Path(base_model)
    adapter = Path(adapter)
    adapter_weight = adapter / "adapter_model.safetensors"
    base_files = []
    for path in sorted(base_model.glob("*.safetensors")):
        stat = path.stat()
        base_files.append({"name": path.name, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    return compute_signature({
        "base_path": str(base_model.resolve()),
        "base_files": base_files,
        "base_config": _sha256(base_model / "config.json"),
        "adapter_sha256": _sha256(adapter_weight),
    })[:24]


def _run_stage4_one(sg: str, query_faa: Path, asset, config: dict, output_dir: Path, gpu: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    query_ids = {record.record_id for record in read_fasta(query_faa)}
    backbone_ids = {record.record_id for record in read_fasta(asset.fasta)}
    collision = query_ids & backbone_ids
    if collision:
        raise ValueError(f"{sg} query/backbone ID 冲突: {sorted(collision)[:5]}")
    if not query_ids:
        raise ValueError(f"{sg} query FASTA 为空")
    cache_key = phylo_cache_key(config["assets"]["phylo_base_model"], asset.adapter)
    cache_dir = Path(config["paths"]["cache_dir"]) / "lucaphylo" / f"{sg}_{cache_key}"
    matrix = output_dir / "stage4_full_distance_matrix.csv"
    if matrix.exists():
        matrix.unlink()
    run_command(
        [
            config["environments"]["lucaphylo"], PROJECT_ROOT / "third_party" / "lucaphylo" / "infer_stage_4.py",
            "--model_name_or_path", config["assets"]["phylo_base_model"],
            "--stage4_adapter_path", asset.adapter,
            "--backbone_tree", asset.tree, "--backbone_seqs", asset.fasta,
            "--query_seqs", query_faa, "--output_dir", output_dir,
            "--max_seq_length", config["phylo"]["max_seq_length"],
            "--batch_size", config["phylo"]["batch_size"],
            "--distance_batch_size", config["phylo"]["distance_batch_size"],
            "--tree_distance_backend", config["phylo"]["tree_distance_backend"],
            "--" + "emb" + "edding_cache_dir", cache_dir,
            "--device", "cuda:0",
        ],
        env_name="lucaphylo",
        cwd=PROJECT_ROOT / "third_party" / "lucaphylo",
        log_path=output_dir / "stage4.log",
        env={"CUDA_VISIBLE_DEVICES": str(gpu)},
    )
    if not matrix.is_file():
        raise FileNotFoundError(f"{sg} 未产出距离矩阵: {matrix}")
    return matrix


def run_phylo(supergroup_groups: dict[str, str | Path], config: dict, workdir: str | Path, gpu: str | int) -> PhyloResult:
    from cdc_pipeline.supergroup_mapping import lucaphylo_to_asset_name
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    assets = load_supergroup_map(config)
    placements: list[Placement] = []
    matrices = {}
    errors = {}
    selected_gpu = pick_free_gpu(gpu) if supergroup_groups else str(gpu)
    for sg, query_path in sorted(supergroup_groups.items()):
        try:

            asset_name = lucaphylo_to_asset_name(sg)
            if asset_name not in assets:
                raise KeyError(f"Supergroup '{sg}' (asset '{asset_name}') has no backbone/adapter assets")
            output_dir = workdir / sg.replace("/", "_")
            matrix = _run_stage4_one(sg, Path(query_path), assets[asset_name], config, output_dir, selected_gpu)
            backbone_ids = {record.record_id for record in read_fasta(assets[asset_name].fasta)}
            current = placements_from_matrix(matrix, backbone_ids, config["phylo"]["top_k"], sg)
            placements.extend(current)
            matrices[sg] = matrix
            write_placements_tsv(output_dir / "placements.tsv", current, config["phylo"]["top_k"])
        except Exception as exc:
            errors[sg] = str(exc)
            (workdir / f"{sg.replace('/', '_')}_error.txt").write_text(str(exc) + "\n")
    if supergroup_groups and not matrices:
        raise RuntimeError("所有超群 LucaPhylo 均失败: " + json.dumps(errors, ensure_ascii=False))
    output = workdir / "placements.tsv"
    write_placements_tsv(output, placements, config["phylo"]["top_k"])
    counts = {sg: sum(placement.supergroup == sg for placement in placements) for sg in matrices}
    counts["placed"] = len(placements)
    counts["failed_supergroups"] = len(errors)
    atomic_write_json(workdir / "phylo_summary.json", {"counts": counts, "errors": errors, "matrices": {key: str(value) for key, value in matrices.items()}})
    return PhyloResult(output, matrices, counts, errors, "done" if supergroup_groups else "skipped")


def prewarm_backbones(config: dict) -> dict[str, str]:
    assets = load_supergroup_map(config)
    gpu = pick_free_gpu(config.get("runtime", {}).get("gpu", "auto"))
    outputs = {}
    for sg, asset in assets.items():
        first = next(read_fasta(asset.fasta))
        output_dir = Path(config["paths"]["cache_dir"]) / "prewarm_runs" / sg
        query = output_dir / "prewarm_query.faa"
        write_fasta(query, [FastaRecord("PREWARM_QUERY", first.sequence, "PREWARM_QUERY")])
        outputs[sg] = str(_run_stage4_one(sg, query, asset, config, output_dir, gpu))
    return outputs


def run(context: PipelineContext) -> StageResult:

    groups = dict(context.data.get("supergroup_groups", {}))
    for sg, path in context.data.get("classifier_groups_with_tree", {}).items():
        if sg not in groups:
            groups[sg] = path
        else:

            from cdc_pipeline.fasta import read_fasta as _read
            existing = list(_read(groups[sg]))
            new_seqs = list(_read(path))
            seen = {r.header_id for r in existing}
            merged = existing + [r for r in new_seqs if r.header_id not in seen]
            merged_path = context.stage_dir / f"merged_{sg.replace('/', '_')}.faa"
            from cdc_pipeline.fasta import write_fasta as _write
            _write(merged_path, merged)
            groups[sg] = merged_path

    result = run_phylo(groups, context.config, context.stage_dir, context.gpu)
    artifacts = {
        "placements": artifact_path(context, result.placements_tsv),
        "summary": artifact_path(context, context.stage_dir / "phylo_summary.json"),
        "matrices": {sg: artifact_path(context, path) for sg, path in result.matrices.items()},
    }
    for sg in result.errors:
        artifacts[f"error_{sg}"] = artifact_path(context, context.stage_dir / f"{sg.replace('/', '_')}_error.txt")
    return StageResult(
        artifacts=artifacts,
        summary={**result.counts, "errors": result.errors, "reason": "" if groups else "no_assigned_supergroup"},
        data={
            "placements_tsv": str(result.placements_tsv),
            "phylo_matrices": {sg: str(path) for sg, path in result.matrices.items()},
            "phylo_counts": result.counts,
            "phylo_errors": result.errors,
        },
        status=result.status,
    )
