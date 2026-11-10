#!/usr/bin/env python3
"""同步报告前端资产（Chart.js + Tabulator）。

报告 HTML 需要把前端库内嵌进文件才能离线自包含。仓库默认已携带这些资产；
如果资产缺失（例如从 clean clone 部署、或误删），运行本脚本重新下载：

    python scripts/sync_report_assets.py          # 缺失时下载
    python scripts/sync_report_assets.py --force  # 强制重新下载

版本固定，避免上游更新导致报告渲染行为变化。
"""

from __future__ import annotations

import argparse
import shutil
import sys
import urllib.request
from pathlib import Path

ASSET_DIR = Path(__file__).resolve().parent.parent / "third_party" / "report_assets"

SOURCES = {
    "chart.umd.min.js": "https://cdn.jsdelivr.net/npm/chart.js@4.4.7/dist/chart.umd.min.js",
    "tabulator.min.js": "https://cdn.jsdelivr.net/npm/tabulator-tables@6.7.1/dist/js/tabulator.min.js",
    "tabulator.min.css": "https://cdn.jsdelivr.net/npm/tabulator-tables@6.7.1/dist/css/tabulator.min.css",
}


MIN_BYTES = {
    "chart.umd.min.js": 100_000,
    "tabulator.min.js": 200_000,
    "tabulator.min.css": 10_000,
}


def _is_valid(path: Path, name: str) -> bool:
    return path.is_file() and path.stat().st_size >= MIN_BYTES[name]


def _download(name: str, url: str, target: Path) -> None:
    print(f"  下载 {name} <- {url}", flush=True)
    temporary = target.with_name(target.name + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    if temporary.stat().st_size < MIN_BYTES[name]:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"{name} 下载内容过小，疑似失败")
    temporary.replace(target)


def main() -> int:
    parser = argparse.ArgumentParser(description="同步 Chart.js/Tabulator 报告资产")
    parser.add_argument("--force", action="store_true", help="强制重新下载")
    args = parser.parse_args()

    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    print(f"资产目录: {ASSET_DIR}", flush=True)
    failed = []
    for name, url in SOURCES.items():
        target = ASSET_DIR / name
        if not args.force and _is_valid(target, name):
            print(f"  已存在 {name} ({target.stat().st_size // 1024} KB)，跳过", flush=True)
            continue
        try:
            _download(name, url, target)
            print(f"  完成 {name} ({target.stat().st_size // 1024} KB)", flush=True)
        except Exception as exc:
            failed.append(name)
            print(f"  失败 {name}: {exc}", flush=True)

    if failed:
        print(f"\n有 {len(failed)} 个资产下载失败: {', '.join(failed)}", file=sys.stderr)
        print("可在有网络的环境运行本脚本后，把 third_party/report_assets/ 一并部署。", file=sys.stderr)
        return 1
    print("\n全部资产就绪。", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
