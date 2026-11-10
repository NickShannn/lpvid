#!/usr/bin/env python3
"""只重跑 visualization + report 两个阶段（跳过上游签名检查）。

用于"树浏览器/报告前端改动，但不想重跑数据"的场景。

用法：
    python scripts/rerun_visualization_and_report.py runs/real_data/<样本> [...]
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from cdc_pipeline.pipeline import ROUTES, PipelineContext, stage_signature
from cdc_pipeline.project import project_metadata
from cdc_pipeline.state import RunState


def _context(state: RunState, run_dir: Path) -> PipelineContext:
    manifest = state.manifest
    config = dict(manifest["config"])
    config.setdefault("project", project_metadata())
    return PipelineContext(
        route=manifest["route"],
        sample=manifest.get("sample", "sample"),
        inputs=manifest.get("inputs", {}),
        config=config,
        run_dir=run_dir,
        state=state,
        gpu=manifest.get("gpu", "auto"),
    )


def rerun(run_dir: Path) -> str:
    state = RunState(run_dir)
    manifest = state.manifest
    route = manifest["route"]
    stages = ROUTES[route]
    context = _context(state, run_dir)

    for name in stages:
        if name in {"visualization", "report"}:
            continue
        context.data.update(state.snapshot["stages"].get(name, {}).get("context_data", {}))

    results = []
    for name in ("visualization", "report"):
        if name not in stages:
            continue
        index = stages.index(name)
        stage_dir = run_dir / f"step{index + 1:02d}_{name}"
        stage_dir.mkdir(parents=True, exist_ok=True)
        context.stage_dir = stage_dir
        module = importlib.import_module(f"cdc_pipeline.stages.{name}")
        signature = stage_signature(context, name, module.__file__)
        state.begin_stage(name, signature)
        result = module.run(context)
        state.complete_stage(
            name, signature, result.artifacts,
            status=result.status, summary=result.summary, context_data=result.data,
        )
        context.data.update(result.data)
        results.append(f"{name}={result.status}")
    state.finish("done")
    return ", ".join(results)


def main(argv: list[str]) -> int:
    if not argv[1:]:
        print("用法: python scripts/rerun_visualization_and_report.py <run_dir> ...", file=sys.stderr)
        return 2
    for item in argv[1:]:
        run_dir = Path(item).resolve()
        if not (run_dir / "state.json").is_file():
            print(f"跳过（无 state.json）: {run_dir}")
            continue
        try:
            print(f"{run_dir.name}: {rerun(run_dir)}")
        except Exception as exc:
            print(f"{run_dir.name}: FAILED {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
