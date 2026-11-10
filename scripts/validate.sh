#!/usr/bin/env bash
# 静态快速验证：不启动大规模模型推理，不修改共享环境。
# 资源包（resources/）未配置时 doctor/prepare 仅提示，Python/JS 测试仍然执行。
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${LPVID_PYTHON:-python3}"
cd "$ROOT"

if ./lpvid doctor; then
  ./lpvid prepare
else
  echo "doctor 报告资源缺失；跳过 prepare，继续执行静态测试。" >&2
fi

PYTHONNOUSERSITE=1 "$PY" -m pytest tests -q || exit 1
node third_party/placement_view/test_core.js || exit 1
printf '\nlpvid static validation: PASS\n'
