"""风险提示板块：类群归属 + 先验风险信息（非风险判定）。

设计原则（与用户对齐）：
- 不判定"这条序列有风险"，而是"这条序列属于某个类群，该类群先验上含有风险病毒，提示关注"。
- 模型可能出错、系统发育放置可能失准，因此只做类群归属展示 + 先验风险提示。
- 两个名单并列：病原微生物目录（危害程度）与法定传染病（甲乙丙类）。
- 不做物种级判定（不把普通冠状病毒判为新冠）；只做科/属级归属 + 提示。
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class TaxonomyEntry:
    tip_id: str
    family: str
    genus: str
    species: str
    cluster_taxon: str
    human: bool | None
    vertebrate: bool | None
    plant: bool | None


@dataclass(frozen=True)
class PathogenEntry:
    seq: str
    name_cn: str
    name_en: str
    family_cn: str
    family_en: str
    category: str
    statutory: str
    covered: str


@dataclass(frozen=True)
class NeighborEvidence:
    """一个真实存在的 nn1-nn5 列的结构化近邻证据。"""

    rank: int
    reference_id: str
    distance: float | None
    family: str
    genus: str
    human: bool | None
    vertebrate: bool | None
    plant: bool | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "rank": self.rank,
            "reference_id": self.reference_id,
            "distance": self.distance,
            "family": self.family,
            "genus": self.genus,
            "human": self.human,
            "vertebrate": self.vertebrate,
            "plant": self.plant,
        }


@dataclass(frozen=True)
class RiskRow:
    query_id: str
    original_id: str
    contig_id: str
    supergroup: str
    nn1: str
    nn1_dist: float | None
    family: str
    genus: str
    host_human: bool | None
    host_vertebrate: bool | None
    host_plant: bool | None
    pathogen_category: str
    pathogen_examples: str
    statutory_category: str
    statutory_examples: str
    hint: str
    neighbors: list[NeighborEvidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "query_id": self.query_id,
            "original_id": self.original_id,
            "contig_id": self.contig_id,
            "supergroup": self.supergroup,
            "nn1": self.nn1,
            "nn1_dist": self.nn1_dist,
            "family": self.family,
            "genus": self.genus,
            "host_human": self.host_human,
            "host_vertebrate": self.host_vertebrate,
            "host_plant": self.host_plant,
            "pathogen_category": self.pathogen_category,
            "pathogen_examples": self.pathogen_examples,
            "statutory_category": self.statutory_category,
            "statutory_examples": self.statutory_examples,
            "hint": self.hint,
            "neighbors": [neighbor.to_dict() for neighbor in self.neighbors],
        }


def parse_taxonomy_flag(value: Any) -> bool | None:
    """taxonomy 宿主标记：仅显式 0/1 映射 False/True；缺值、空值、非法值为 None。"""
    if value is None:
        return None
    text = str(value).strip()
    if text == "1":
        return True
    if text == "0":
        return False
    return None


def load_taxonomy(taxonomy_dir: Path) -> dict[str, TaxonomyEntry]:
    """加载骨架树 taxonomy：Cluster_sseqid → TaxonomyEntry。"""
    result: dict[str, TaxonomyEntry] = {}
    for path in sorted(taxonomy_dir.glob("*.csv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                tip_id = row.get("Cluster_sseqid", "").strip()
                if not tip_id:
                    continue
                result[tip_id] = TaxonomyEntry(
                    tip_id=tip_id,
                    family=row.get("Family", "").strip(),
                    genus=row.get("Genus", "").strip(),
                    species=row.get("Species", "").strip(),
                    cluster_taxon=row.get("Cluster_taxon", "").strip(),
                    human=parse_taxonomy_flag(row.get("Human")),
                    vertebrate=parse_taxonomy_flag(row.get("Vertebrates")),
                    plant=parse_taxonomy_flag(row.get("Plant")),
                )
    return result


def load_pathogen_catalog(path: Path) -> list[PathogenEntry]:
    """加载病原微生物目录 + 法定传染病清单。"""
    entries: list[PathogenEntry] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            entries.append(PathogenEntry(
                seq=row.get("seq", ""),
                name_cn=row.get("name_cn", ""),
                name_en=row.get("name_en", ""),
                family_cn=row.get("family_cn", ""),
                family_en=row.get("family_en", ""),
                category=row.get("category", ""),
                statutory=row.get("statutory", ""),
                covered=row.get("covered", ""),
            ))
    return entries


def _catalog_by_family(entries: list[PathogenEntry]) -> dict[str, list[PathogenEntry]]:
    grouped: dict[str, list[PathogenEntry]] = defaultdict(list)
    for entry in entries:
        grouped[entry.family_en].append(entry)
    return grouped


def _pick_category(entries: list[PathogenEntry]) -> tuple[str, str, str, str]:
    """从同科条目中取最高危害程度与法定类别，并给出代表病毒/疾病。"""
    if not entries:
        return "", "", "", ""
    order = {"第一类": 0, "第二类": 1, "第三类": 2, "第四类": 3, "": 9}
    statutory_order = {"甲类": 0, "乙类": 1, "丙类": 2, "": 9}
    by_cat = sorted(entries, key=lambda e: order.get(e.category, 9))
    top_cat = by_cat[0].category
    cat_examples = [e.name_cn for e in by_cat if e.category == top_cat][:3]
    stat_entries = [e for e in entries if e.statutory]
    if stat_entries:
        stat_sorted = sorted(stat_entries, key=lambda e: statutory_order.get(e.statutory, 9))
        top_stat = stat_sorted[0].statutory
        stat_examples = [e.name_cn for e in stat_sorted if e.statutory == top_stat][:3]
    else:
        stat_sorted = []
        top_stat, stat_examples = "", []
    return (
        top_cat,
        "、".join(cat_examples) + ("等" if len([e for e in by_cat if e.category == top_cat]) > 3 else ""),
        top_stat,
        "、".join(stat_examples) + ("等" if len([e for e in stat_sorted if e.statutory == top_stat]) > 3 else ""),
    )


_FOLLOWUP = "需进一步核查，不构成新序列风险分级或病原确证。"


def _nn_distance(placement: dict[str, Any], rank: int) -> float | None:
    """读取 nnN_dist；缺失、空串、非法、NaN、Infinity 均返回 None。"""
    value = placement.get(f"nn{rank}_dist")
    if value is None:
        return None
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
    try:
        distance = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(distance):
        return None
    return distance


def _collect_neighbors(
    placement: dict[str, Any],
    taxonomy: dict[str, TaxonomyEntry],
) -> list[NeighborEvidence]:
    """按 nn1..nn5 真实列顺序收集非空近邻；空白列跳过但不重排序号。

    taxonomy 缺失的近邻仍保留 reference_id/distance，family/genus 为空，
    host flags 为 None 表示 unknown，而不是阴性。
    """
    neighbors: list[NeighborEvidence] = []
    for rank in range(1, 6):
        reference_id = str(placement.get(f"nn{rank}", "") or "").strip()
        if not reference_id:
            continue
        tax = taxonomy.get(reference_id)
        neighbors.append(NeighborEvidence(
            rank=rank,
            reference_id=reference_id,
            distance=_nn_distance(placement, rank),
            family=tax.family if tax else "",
            genus=tax.genus if tax else "",
            human=tax.human if tax else None,
            vertebrate=tax.vertebrate if tax else None,
            plant=tax.plant if tax else None,
        ))
    return neighbors


def _hint(
    family: str,
    category: str,
    category_examples: str,
    statutory: str,
    statutory_examples: str,
) -> str:
    """提示仅描述 nn1 所属科的已知成员/名录收录成员，不把整个科定级。

    同时展示病原微生物目录与法定传染病两个维度（若都存在）。
    """
    subject = f"{family}（nn1 所属科）" if family else "nn1 所属科未识别"
    dimensions: list[str] = []
    if category:
        dimensions.append(f"《病原微生物目录》{category}（{category_examples}）")
    if statutory:
        dimensions.append(f"《传染病防治法》{statutory}（{statutory_examples}）")
    if not dimensions:
        separator = "，" if not family else ""
        return (
            f"{subject}{separator}未匹配到上述名录收录成员；"
            "不排除名录外或未识别病毒，"
            f"{_FOLLOWUP}"
        )
    return f"{subject}的已知成员中，有名录收录成员涉及" + "和".join(dimensions) + f"；{_FOLLOWUP}"


def build_risk_rows(
    placements: list[dict[str, Any]],
    taxonomy: dict[str, TaxonomyEntry],
    catalog: list[PathogenEntry],
    *,
    original_ids: dict[str, str] | None = None,
    contig_map: dict[str, str] | None = None,
) -> list[RiskRow]:
    """为每条 placed query 生成类群归属 + 先验风险提示。

    每条 RiskRow 的 neighbors 按真实列序保存 nn1-nn5 的结构化近邻证据；
    病原目录/法定名单的科级匹配仅锚定 nn1 的科，不受 nn2-nn5 影响。

    参数：
    - placements: results/placements.tsv 的行（query_id, supergroup, nn1..nn5 及对应距离）
    - taxonomy: load_taxonomy() 的结果
    - catalog: load_pathogen_catalog() 的结果
    - original_ids: query_id → 原始 ORF ID
    - contig_map: query_id → contig_id
    """
    by_family = _catalog_by_family(catalog)
    original_ids = original_ids or {}
    contig_map = contig_map or {}
    rows: list[RiskRow] = []

    for placement in placements:
        query_id = placement.get("query_id", "")
        nn1 = str(placement.get("nn1", "") or "").strip()
        if not nn1:
            continue
        tax = taxonomy.get(nn1)
        family = tax.family if tax else ""
        genus = tax.genus if tax else ""


        entries = by_family.get(family, []) if family else []
        cat, cat_ex, stat, stat_ex = _pick_category(entries)
        neighbors = _collect_neighbors(placement, taxonomy)
        hint = _hint(family, cat, cat_ex, stat, stat_ex)

        rows.append(RiskRow(
            query_id=query_id,
            original_id=original_ids.get(query_id, ""),
            contig_id=contig_map.get(query_id, ""),
            supergroup=placement.get("supergroup", ""),
            nn1=nn1,
            nn1_dist=_nn_distance(placement, 1),
            family=family,
            genus=genus,
            host_human=tax.human if tax else None,
            host_vertebrate=tax.vertebrate if tax else None,
            host_plant=tax.plant if tax else None,
            pathogen_category=cat,
            pathogen_examples=cat_ex,
            statutory_category=stat,
            statutory_examples=stat_ex,
            hint=hint,
            neighbors=neighbors,
        ))
    return rows


def summarize_risk(rows: list[RiskRow]) -> dict[str, Any]:
    """汇总统计：各级别数量、科级分布、法定传染病命中、互斥宿主证据计数。"""
    level_count = defaultdict(int)
    family_count = defaultdict(int)
    statutory_count = defaultdict(int)
    human_nn1 = 0
    human_nn2_5 = 0
    vertebrate_only = 0
    no_placed_label = 0
    for row in rows:
        level = "statutory" if row.statutory_category else ("pathogen" if row.pathogen_category else "other")
        level_count[level] += 1
        family_count[row.family] += 1
        if row.statutory_category:
            statutory_count[row.statutory_category] += 1





        if not row.neighbors:

            if row.host_human is True:
                human_nn1 += 1
            elif row.host_vertebrate is True:
                vertebrate_only += 1
            else:
                no_placed_label += 1
            continue
        nn1_evidence = next((n for n in row.neighbors if n.rank == 1), None)
        if nn1_evidence is not None and nn1_evidence.human is True:
            human_nn1 += 1
        elif any(n.human is True for n in row.neighbors if 2 <= n.rank <= 5):
            human_nn2_5 += 1
        elif any(n.vertebrate is True for n in row.neighbors if 1 <= n.rank <= 5):
            vertebrate_only += 1
        else:
            no_placed_label += 1
    return {
        "total": len(rows),
        "by_level": dict(level_count),
        "by_family": dict(sorted(family_count.items(), key=lambda x: -x[1])),
        "by_statutory": dict(statutory_count),
        "human_nn1": human_nn1,
        "human_nn2_5": human_nn2_5,
        "vertebrate_only": vertebrate_only,
        "no_placed_label": no_placed_label,
    }
