"""LucaPhylo 距离矩阵到 top-k 放置表和可视化记录。"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.stages.common import write_tsv


@dataclass(frozen=True)
class Neighbor:
    tip_id: str
    distance: float


@dataclass
class Placement:
    query_id: str
    supergroup: str
    neighbors: list[Neighbor]
    method: str = "lucaphylo_stage4"


def summarize_distances(distances: dict[str, float], backbone_ids: set[str], top_k: int) -> list[Neighbor]:
    """返回按距离升序排列的 top-k backbone 近邻。只保留 LucaPhylo 实际输出的距离。"""
    candidates = sorted((distance, tip_id) for tip_id, distance in distances.items() if tip_id in backbone_ids)
    selected = candidates[:top_k]
    if not selected:
        raise ValueError("query 没有可用 backbone 近邻")
    return [Neighbor(tip_id, distance) for distance, tip_id in selected]


def placements_from_matrix(
    matrix_path: str | Path,
    backbone_ids: set[str],
    top_k: int = 5,
    supergroup: str = "",
) -> list[Placement]:
    placements = []
    with Path(matrix_path).open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        column_ids = header[1:]
        backbone_columns = [(index, tip_id) for index, tip_id in enumerate(column_ids) if tip_id in backbone_ids]
        for row in reader:
            if not row or row[0] in backbone_ids:
                continue
            distances = {tip_id: float(row[index + 1]) for index, tip_id in backbone_columns}
            neighbors = summarize_distances(distances, backbone_ids, top_k)
            placements.append(Placement(row[0], supergroup, neighbors))
    return placements


def write_placements_tsv(path: str | Path, placements: list[Placement], top_k: int = 5) -> Path:
    """只输出 LucaPhylo 实际提供的指标：query_id、超群、top-k 近邻及其距离。"""
    fields = ["query_id", "supergroup"]
    for index in range(1, top_k + 1):
        fields += [f"nn{index}", f"nn{index}_dist"]
    fields += ["method"]
    rows = []
    for placement in placements:
        row = {"query_id": placement.query_id, "supergroup": placement.supergroup, "method": placement.method}
        for index, neighbor in enumerate(placement.neighbors, 1):
            row[f"nn{index}"] = neighbor.tip_id
            row[f"nn{index}_dist"] = neighbor.distance
        rows.append(row)
    return write_tsv(path, fields, rows)
