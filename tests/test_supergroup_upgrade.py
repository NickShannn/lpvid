"""180 超群分配测试：参考库、名称映射与分类模型。"""

import csv
from pathlib import Path

import pytest

from cdc_pipeline.config import PROJECT_ROOT, load_config
from cdc_pipeline.supergroup_mapping import (
    get_tree_availability,
    load_supergroup_mapping,
)


def test_supergroups_mapping_loads_correctly():
    """用户提供的 24→21 映射表可加载，关键字段完整。"""
    mapping = load_supergroup_mapping()
    assert len(mapping) == 24  # 24 行原始 VirID 超群
    # 检查合并映射
    assert mapping["Qinvirus"] == "Qin-Yue"
    assert mapping["Yuevirus"] == "Qin-Yue"
    assert mapping["Zhaovirus"] == "Tombus-Noda"
    assert mapping["na-daxi"] == "Quenyavirus"
    # other-negative 应为空
    assert mapping.get("other-negative", "") == ""


REF_DB = Path(load_config()["assets"]["rdrp_180_db"])
REF_CSV = Path(load_config()["assets"]["rdrp_180_csv"])
PHYLO_ROOT = Path(load_config()["assets"]["phylo_root"])
CLASSIFIER_DIR = Path(load_config()["assets"]["classifier_pipeline"])


@pytest.mark.skipif(not REF_DB.is_file(), reason="外部 RdRP 参考库资源包未配置")
def test_ref_database_exists():
    """180 超群参考库存在。"""
    config = load_config()
    db = Path(config["assets"]["rdrp_180_db"])
    assert db.is_file(), f"Missing: {db}"


@pytest.mark.skipif(not REF_CSV.is_file(), reason="外部 RdRP 参考标签表未配置")
def test_ref_mapping_covers_all_named_supergroups():
    """参考库 CSV 的 super_group 列可映射到 LucaPhylo_SGs。"""
    mapping = load_supergroup_mapping()
    # 构建 sseqid → LucaPhylo_SGs 映射
    ref_csv = REF_CSV
    with ref_csv.open() as f:
        rows = list(csv.DictReader(f))
    # 检查参考库中有已命名超群
    named = {row["super_group"] for row in rows if not row["super_group"].startswith("Supergroup")}
    assert "Astro-Poty" in named
    assert "Bunya-Arena" in named
    assert "Picorna" in named
    assert "Qin-Yue" in named
    assert "Tombus-Noda" in named
    # na-daxi 需要映射到 Quenyavirus
    assert "na-daxi" in named
    # other-negative 不应出现在映射后
    assert "other-negative" not in {mapping.get(v, v) for v in named if v != "other-negative"}


@pytest.mark.skipif(not PHYLO_ROOT.is_dir(), reason="外部超群骨架树资源包未配置")
def test_tree_availability_check():
    """检查 21 个 LucaPhylo_SGs 中哪些有骨架树。"""
    availability = get_tree_availability()
    # 18 个直接映射的超群应有树
    assert availability["Astro-Poty"] == "available"
    assert availability["Picorna"] == "available"
    assert availability["Flavi"] == "available"
    # 合并超群暂时无树
    assert availability["Qin-Yue"] == "pending_merge"
    assert availability["Tombus-Noda"] == "pending_merge"
    # Quenyavirus 已从 na-daxi 重命名
    assert availability["Quenyavirus"] == "available"
    # 无树的新超群
    assert availability.get("Supergroup022") == "no_tree"
    assert availability.get("Supergroup100") == "no_tree"


@pytest.mark.skipif(not (CLASSIFIER_DIR / "classify.py").is_file(), reason="外部分类器资源未配置")
def test_classifier_paths():
    """分类模型所需文件路径存在。"""
    config = load_config()
    classifier_dir = Path(config["assets"]["classifier_pipeline"])
    assert (classifier_dir / "classify.py").is_file()
    assert (classifier_dir / "labels_supergroups.tsv").is_file()
