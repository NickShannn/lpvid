import csv
import json
from pathlib import Path

import pytest

from cdc_pipeline.assets import load_supergroup_map, normalize_diamond_sequence, resolve_adapter_base
from cdc_pipeline.config import PROJECT_ROOT, load_config


def test_diamond_sequence_normalization_is_explicit():
    normalized, counts = normalize_diamond_sequence("AC-?DEFG**")
    assert normalized == "ACXDEFG"
    assert counts == {"gaps_removed": 1, "unknown_replaced": 1, "terminal_stops_removed": 2, "internal_stops_replaced": 0}


def test_adapter_base_path_is_portable():
    relative = "resources/models/lucaphylo/checkpoint-81000"
    assert resolve_adapter_base(relative) == (PROJECT_ROOT / relative).resolve()
    absolute = "/opt/lpvid/checkpoint-81000"
    assert resolve_adapter_base(absolute) == Path(absolute).resolve()


@pytest.mark.skipif(not (PROJECT_ROOT / "assets" / "phylo").is_dir(), reason="外部 24 超群资源包未配置")
def test_24_supergroups_have_tree_fasta_taxonomy_and_adapter():
    assets = load_supergroup_map(load_config(PROJECT_ROOT / "config" / "local.yaml"))
    # 24 原始 + Quenyavirus 重命名副本 = 25
    assert len(assets) >= 24
    for sg, asset in assets.items():
        assert asset.tree.is_file(), sg
        assert asset.fasta.is_file(), sg
        assert asset.taxonomy.is_file(), sg
        assert (asset.adapter / "adapter_model.safetensors").is_file(), sg
        assert (asset.adapter / "adapter_config.json").is_file(), sg


@pytest.mark.skipif(not (PROJECT_ROOT / "assets" / "supergroup24" / "reference.dmnd").is_file(), reason="本地 24 超群 DIAMOND 库尚未 prepare")
def test_prepared_reference_has_24_supergroups_and_unique_ids():
    root = PROJECT_ROOT / "assets" / "supergroup24"
    mapping = root / "mapping.tsv"
    assert mapping.is_file()
    with mapping.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len({row["supergroup"] for row in rows}) == 24
    assert len({row["reference_id"] for row in rows}) == len(rows)
    assert (root / "reference.dmnd").stat().st_size > 0
    validation = json.loads((root / "validation.json").read_text())
    assert validation["supergroup_count"] == 24
    assert validation["tree_fasta_id_mismatches"] == []
    assert validation["taxonomy_missing"] == []
