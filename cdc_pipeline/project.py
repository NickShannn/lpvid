"""项目展示名、发行名、命令名和版本的唯一读取入口。"""

from __future__ import annotations

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
METADATA_PATH = PROJECT_ROOT / "config" / "project.json"


def project_metadata() -> dict[str, str]:
    with METADATA_PATH.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    required = {"display_name", "distribution", "command", "version", "description"}
    missing = required - set(metadata)
    if missing:
        raise ValueError(f"项目元数据缺少字段: {sorted(missing)}")
    return metadata
