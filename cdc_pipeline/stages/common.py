"""阶段共享的产物和表格工具。"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable


def artifact_path(context, path: str | Path) -> str:
    path = Path(path)
    try:
        return str(path.relative_to(context.run_dir))
    except ValueError:
        return str(path)


def write_tsv(path: str | Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
