#!/usr/bin/env python3
"""从显式 bundle JSON 列表构建单文件离线树，不扫描旧演示数据。"""

import argparse
import json
from pathlib import Path

from cdc_pipeline.visualization import build_tree_html


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundles", required=True, help="包含 bundle 列表的 JSON 文件")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    bundles = json.loads(Path(args.bundles).read_text())
    build_tree_html(bundles, args.output)
    print(args.output)


if __name__ == "__main__":
    main()
