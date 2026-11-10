#!/usr/bin/env python3
"""生成外部资源包清单（路径 / 大小 / SHA256），供发布前核验。

资源包不属于源码仓库，由项目方单独分发；本清单用于：
  - 记录每个资产的来源与校验值
  - 让外部用户核对下载完整性
  - 在每个 run_manifest.json 中记录资源版本

用法：
    python scripts/build_resource_manifest.py                 # 扫描默认资源目录
    python scripts/build_resource_manifest.py --out FILE      # 指定输出
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESOURCE_DIRS = ["assets/phylo", "assets/RdRP_blast_ref", "assets/supergroup24"]
DEFAULT_OUT = ROOT / "assets" / "RESOURCE_PACKAGE_MANIFEST.json"
CHUNK = 8 * 1024 * 1024


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect(directories: list[str]) -> list[dict]:
    entries: list[dict] = []
    for relative in directories:
        base = ROOT / relative
        if not base.is_dir():
            print(f"跳过缺失目录: {relative}", file=sys.stderr)
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            stat = path.stat()
            entries.append({
                "path": str(path.relative_to(ROOT)),
                "size": stat.st_size,
                "sha256": sha256_of(path),
            })
    return entries


def main() -> int:
    parser = argparse.ArgumentParser(description="生成资源包发布清单")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dirs", nargs="*", default=RESOURCE_DIRS)
    args = parser.parse_args()

    print(f"扫描: {', '.join(args.dirs)}", flush=True)
    entries = collect(args.dirs)
    total = sum(entry["size"] for entry in entries)
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "directories": args.dirs,
        "file_count": len(entries),
        "total_bytes": total,
        "total_human": f"{total / 1024 ** 3:.2f} GB" if total >= 1024 ** 3 else f"{total / 1024 ** 2:.1f} MB",
        "files": entries,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"完成: {len(entries)} 个文件, {payload['total_human']} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
