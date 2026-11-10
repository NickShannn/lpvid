"""历史 24 超群到 180 超群体系的名称映射和骨架树可用性检查。"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from cdc_pipeline.assets import load_supergroup_map


def load_supergroup_mapping() -> dict[str, str]:
    """加载历史 24 超群到 180 超群名称的映射。"""
    from cdc_pipeline.config import load_config
    config = load_config()
    path = Path(config["assets"]["supergroups_mapping"])
    mapping = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            legacy_name = row["VirID_SGs"].strip()
            lucaphylo = row["LucaPhylo_SGs"].strip()
            if legacy_name:
                mapping[legacy_name] = lucaphylo
    return mapping


def load_reference_supergroup_map() -> dict[str, str]:
    """从 180 超群参考库 CSV 构建 sseqid→supergroup 映射。"""
    from cdc_pipeline.config import load_config
    config = load_config()
    path = Path(config["assets"]["rdrp_180_csv"])
    result = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):

            sseqid = row.get("RdRp_sseqid", row.get("RdRP_sseqid", "")).strip()
            supergroup = row.get("super_group", row.get("Super_group", "")).strip()
            if sseqid and supergroup:
                result[sseqid] = supergroup
    return result


def normalize_supergroup_name(name: str, mapping: dict[str, str] | None = None) -> str:
    """将参考库名称标准化为 LucaPhylo_SGs 名称。

    - na-daxi → Quenyavirus
    - other-negative → 空字符串（删除）
    - 已命名的直接返回
    """
    if mapping is None:
        mapping = load_supergroup_mapping()

    if name in mapping:
        return mapping[name]

    if name == "na-daxi":
        return "Quenyavirus"

    if name == "other-negative":
        return ""

    return name


def get_tree_availability() -> dict[str, str]:
    """检查每个超群的骨架树可用性。

    返回值：
    - available: 有现成树和 adapter
    - pending_merge: 需要合并建树（Qin-Yue、Tombus-Noda 含 Zhaovirus）
    - pending_rename: 需要重命名（Quenyavirus ← na-daxi）
    - no_tree: 无骨架树（大部分新超群）
    """
    from cdc_pipeline.config import load_config
    config = load_config()
    existing = load_supergroup_map(config)
    existing_names = set(existing.keys())
    mapping = load_supergroup_mapping()


    named_supergroups = {v for v in mapping.values() if v}


    merge_sources = {
        "Qin-Yue": {"Qinvirus", "Yuevirus"},
        "Tombus-Noda": {"Tombus-Noda", "Zhaovirus"},
    }


    reverse = {}
    for legacy_name, lucaphylo in mapping.items():
        if lucaphylo and legacy_name in existing_names:
            reverse.setdefault(lucaphylo, set()).add(legacy_name)

    availability = {}
    for sg in sorted(named_supergroups):
        if sg in merge_sources:
            sources = merge_sources[sg]
            if sources.issubset(existing_names):
                availability[sg] = "pending_merge"
            else:
                availability[sg] = "no_tree"
        elif sg == "Quenyavirus":

            if "Quenyavirus" in existing_names or "na-daxi" in existing_names:
                availability[sg] = "available"
            else:
                availability[sg] = "no_tree"
        elif sg in existing_names or sg in reverse:

            availability[sg] = "available"
        else:
            availability[sg] = "no_tree"



    ref_map = load_reference_supergroup_map()
    ref_supergroups = {normalize_supergroup_name(sg, mapping) for sg in set(ref_map.values())}
    for sg in sorted(ref_supergroups):
        if sg and sg.startswith("Supergroup") and sg not in availability:
            availability[sg] = "no_tree"

    return availability


def get_placement_eligible_supergroups() -> set[str]:
    """返回当前可以直接进行 LucaPhylo 放置的超群名。"""
    availability = get_tree_availability()
    return {sg for sg, status in availability.items() if status == "available"}


def lucaphylo_to_asset_name(lucaphylo_sg: str) -> str:
    """将 180 超群名称映射回资产目录使用的历史名称。

    例如 Picorna → Picorna-Calici, Bunya-Arena → Negative_Bunya-Arena。
    Quenyavirus 已重命名，直接匹配。
    """

    if lucaphylo_sg == "Quenyavirus":
        return "Quenyavirus"
    mapping = load_supergroup_mapping()
    for legacy_name, lucaphylo in mapping.items():
        if lucaphylo == lucaphylo_sg:
            return legacy_name
    return lucaphylo_sg
