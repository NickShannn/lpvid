"""配置加载：默认科学参数与本机路径分离。"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

from cdc_pipeline.project import project_metadata

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "default.yaml"
LOCAL_CONFIG = PROJECT_ROOT / "config" / "local.yaml"


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """递归合并配置，不修改调用方对象。"""
    result = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"配置文件不存在: {path}")
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"配置根节点必须是映射: {path}")
    return data


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """加载 default.yaml，再合并显式配置或默认 local.yaml。"""
    config = _read_yaml(DEFAULT_CONFIG)
    override_path = Path(path) if path else (LOCAL_CONFIG if LOCAL_CONFIG.exists() else None)
    if override_path:
        if not override_path.is_absolute():
            override_path = (Path.cwd() / override_path).resolve()
        config = deep_merge(config, _read_yaml(override_path))
    config["paths"] = {
        key: str((PROJECT_ROOT / value).resolve()) if isinstance(value, str) and not Path(value).is_absolute() else value
        for key, value in config.get("paths", {}).items()
    }
    portable_prefixes = ("resources/", "runtime/", "third_party/", "assets/")
    for section in ("tools", "assets"):
        for key, value in config.get(section, {}).items():
            if isinstance(value, str) and not Path(value).is_absolute() and value.startswith(portable_prefixes):
                config[section][key] = str(PROJECT_ROOT / value)
    config["config_path"] = str(override_path) if override_path else str(DEFAULT_CONFIG)
    config["project"] = project_metadata()
    return config


def env_python(config: dict[str, Any], name: str) -> Path:
    try:
        return Path(config["environments"][name])
    except KeyError as exc:
        raise KeyError(f"未配置运行环境: {name}") from exc
