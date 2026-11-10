#!/usr/bin/env python3
"""Prepare in-project runtime compatibility libraries without modifying shared environments."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

from cdc_pipeline.config import load_config

REQUIRED = ["libtiff.so.6", "libopenjp2.so.7", "libLerc.so.4", "libjpeg.so.8"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="本机 YAML 配置")
    parser.add_argument("--source-dir", help="包含兼容库的目录；默认取 lucaphylo 环境的 lib")
    args = parser.parse_args()
    config = load_config(args.config)
    source_dir = Path(args.source_dir) if args.source_dir else Path(config["environments"]["lucaphylo"]).resolve().parent.parent / "lib"
    missing = [name for name in REQUIRED if not (source_dir / name).is_file()]
    if missing:
        raise SystemExit(f"缺少兼容库 {missing}；请用 --source-dir 指定包含它们的目录: {source_dir}")
    target_dir = ROOT / "runtime" / "libs"
    target_dir.mkdir(parents=True, exist_ok=True)
    actions = {}
    for name in REQUIRED:
        source = source_dir / name
        target = target_dir / name
        if target.exists() or target.is_symlink():
            if target.is_symlink() and target.resolve() == source.resolve():
                actions[name] = "unchanged"
                continue
            raise SystemExit(f"目标已存在且与来源不同，拒绝覆盖: {target}")
        target.symlink_to(source)
        actions[name] = "linked"
    print(json.dumps({"source_dir": str(source_dir), "target_dir": str(target_dir), "actions": actions}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
