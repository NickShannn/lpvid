"""24 超群资产校验、参考库构建和 prepare 入口。"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ete3 import Tree

from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.fasta import FastaRecord, read_fasta, write_fasta
from cdc_pipeline.process import run_command
from cdc_pipeline.state import atomic_write_json, compute_signature, file_fingerprint


@dataclass(frozen=True)
class SupergroupAsset:
    name: str
    tree: Path
    fasta: Path
    taxonomy: Path
    adapter: Path


def load_supergroup_map(config: dict[str, Any]) -> dict[str, SupergroupAsset]:
    root = Path(config["assets"]["phylo_root"])
    backbone = root / "backbone_dir"
    if not backbone.is_dir():
        return {}
    assets = {}
    for directory in sorted(path for path in backbone.iterdir() if path.is_dir()):
        trees = sorted(directory.glob("*.nwk"))
        fastas = sorted(directory.glob("*.fas"))
        assets[directory.name] = SupergroupAsset(
            name=directory.name,
            tree=trees[0] if len(trees) == 1 else directory / "MISSING_OR_MULTIPLE.nwk",
            fasta=fastas[0] if len(fastas) == 1 else directory / "MISSING_OR_MULTIPLE.fas",
            taxonomy=root / "taxonomy" / f"{directory.name}.csv",
            adapter=root / "adapter_dir" / directory.name,
        )
    return assets


def resolve_adapter_base(recorded: str | Path) -> Path:
    """Resolve the base-model path recorded in an adapter configuration.

    Adapters shipped in the resource record use repository-relative paths so
    that a downloaded package works from any installation directory. Absolute
    paths remain supported for locally trained adapters.
    """
    path = Path(recorded)
    if path.is_absolute():
        return path.resolve()
    return (PROJECT_ROOT / path).resolve()


def _taxonomy_ids(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {row["Cluster_sseqid"] for row in csv.DictReader(handle)}


def normalize_diamond_sequence(sequence: str) -> tuple[str, dict[str, int]]:
    """仅规范化 DIAMOND 参考副本；原 backbone 序列不修改。"""
    counts = {
        "gaps_removed": sequence.count("-"),
        "unknown_replaced": sequence.count("?"),
        "terminal_stops_removed": len(sequence) - len(sequence.rstrip("*")),
        "internal_stops_replaced": 0,
    }
    normalized = sequence.replace("-", "").replace("?", "X").rstrip("*")
    counts["internal_stops_replaced"] = normalized.count("*")
    normalized = normalized.replace("*", "X")
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    invalid = sorted(set(normalized) - allowed)
    if invalid:
        raise ValueError(f"DIAMOND 参考序列含未登记字符: {invalid}")
    return normalized, counts


def _write_mapping(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["reference_id", "supergroup", "original_id", "source_fasta"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare_supergroup_assets(config: dict[str, Any]) -> dict[str, Any]:
    assets = load_supergroup_map(config)
    if len(assets) != 24:
        raise RuntimeError(f"需要 24 个超群资产，实际发现 {len(assets)}；请按 assets/RESOURCE_LAYOUT.md 配置资源包")
    output = Path(config["paths"]["assets_dir"]) / "supergroup24"
    output.mkdir(parents=True, exist_ok=True)
    reference_records = []
    mapping_rows = []
    mismatches = []
    taxonomy_missing = []
    adapter_mismatches = []
    normalization = {}
    expected_base = str(Path(config["assets"]["phylo_base_model"]).resolve())
    counter = 0
    for sg, asset in assets.items():
        for required in [asset.tree, asset.fasta, asset.taxonomy, asset.adapter / "adapter_config.json", asset.adapter / "adapter_model.safetensors"]:
            if not required.is_file():
                raise FileNotFoundError(f"{sg} 缺少资产: {required}")
        records = list(read_fasta(asset.fasta))
        fasta_ids = {record.record_id for record in records}
        tree_tips = {tip.name for tip in Tree(str(asset.tree), format=1)}
        if tree_tips != fasta_ids:
            mismatches.append({
                "supergroup": sg,
                "tree_only": sorted(tree_tips - fasta_ids),
                "fasta_only": sorted(fasta_ids - tree_tips),
            })
        taxonomy_ids = _taxonomy_ids(asset.taxonomy)
        if not tree_tips <= taxonomy_ids:
            taxonomy_missing.append({"supergroup": sg, "tip_ids": sorted(tree_tips - taxonomy_ids)})
        adapter_config = json.loads((asset.adapter / "adapter_config.json").read_text())
        recorded_base = str(resolve_adapter_base(adapter_config["base_model_name_or_path"]))
        if recorded_base != expected_base:
            adapter_mismatches.append({"supergroup": sg, "adapter_base": recorded_base, "configured_base": expected_base})
        normalization[sg] = {"gaps_removed": 0, "unknown_replaced": 0, "terminal_stops_removed": 0, "internal_stops_replaced": 0}
        for record in records:
            counter += 1
            reference_id = f"R{counter:07d}"
            normalized, counts = normalize_diamond_sequence(record.sequence)
            for key, value in counts.items():
                normalization[sg][key] += value
            reference_records.append(FastaRecord(
                record_id=reference_id,
                sequence=normalized,
                description=f"{reference_id} original_id={record.record_id} supergroup={sg}",
            ))
            mapping_rows.append({
                "reference_id": reference_id,
                "supergroup": sg,
                "original_id": record.record_id,
                "source_fasta": str(asset.fasta),
            })
    if mismatches or taxonomy_missing or adapter_mismatches:
        raise ValueError(json.dumps({
            "tree_fasta_id_mismatches": mismatches,
            "taxonomy_missing": taxonomy_missing,
            "adapter_base_mismatches": adapter_mismatches,
        }, ensure_ascii=False))

    reference = output / "reference.faa"
    mapping = output / "mapping.tsv"
    write_fasta(reference, reference_records)
    _write_mapping(mapping, mapping_rows)
    signature = compute_signature({
        "reference_sha256": file_fingerprint(reference)["sha256"],
        "mapping_sha256": file_fingerprint(mapping)["sha256"],
        "diamond": config["tools"]["diamond"],
    })
    database = output / "reference.dmnd"
    build_path = output / "build.json"
    previous = json.loads(build_path.read_text()) if build_path.exists() else {}
    built = False
    if not database.exists() or previous.get("signature") != signature:
        run_command(
            [config["tools"]["diamond"], "makedb", "--in", reference, "-d", output / "reference", "--threads", "16"],
            env_name="lucaphylo",
            cwd=output,
            log_path=output / "prepare.log",
        )
        built = True
    atomic_write_json(build_path, {"signature": signature, "reference_count": counter})
    validation = {
        "supergroup_count": len(assets),
        "reference_count": counter,
        "tree_fasta_id_mismatches": mismatches,
        "taxonomy_missing": taxonomy_missing,
        "adapter_base_mismatches": adapter_mismatches,
        "sequence_normalization": normalization,
        "database_built": built,
        "database": str(database),
    }
    atomic_write_json(output / "validation.json", validation)
    return {
        "supergroups": len(assets),
        "reference_sequences": counter,
        "database_built": built,
        "reference": str(reference),
        "mapping": str(mapping),
        "database": str(database),
    }


def prepare_assets(config: dict[str, Any], prewarm: bool = False) -> dict[str, Any]:
    result = prepare_supergroup_assets(config)
    if prewarm:
        from cdc_pipeline.stages.phylo import prewarm_backbones
        result["prewarm"] = prewarm_backbones(config)
    return result
