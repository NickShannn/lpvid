#!/usr/bin/env python3
"""只重跑 report 阶段，跳过所有上游阶段的签名/产物检查。

用于"只想更新报告 UI、不想重跑数据"的场景。上游 context_data 直接从
state.json 汇总，results/ 产物已存在即可。

用法：
    python scripts/rerun_report_only.py runs/real_data/<样本> [...]
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


def _base_context(state: RunState, run_dir: Path, manifest: dict) -> PipelineContext:
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


def _visualization_is_current(state: RunState, context: PipelineContext) -> bool:
    """用 visualization 之前阶段的 context 计算当前签名，禁止旧树被静默复用。"""
    route = context.route
    stages = ROUTES[route]
    if "visualization" not in stages:
        return False
    visualization_index = stages.index("visualization")
    for name in stages[:visualization_index]:
        context.data.update(state.snapshot["stages"].get(name, {}).get("context_data", {}))
    module = importlib.import_module("cdc_pipeline.stages.visualization")
    signature = stage_signature(context, "visualization", module.__file__)
    return state.should_skip("visualization", signature)


def rerun(run_dir: Path) -> str:
    state = RunState(run_dir)
    manifest = state.manifest
    route = manifest["route"]
    context = _base_context(state, run_dir, manifest)

    if not _visualization_is_current(state, _base_context(state, run_dir, manifest)):
        raise RuntimeError(
            "visualization 阶段签名或产物不是当前状态，拒绝只重跑 report；"
            "请先运行 scripts/rerun_visualization_and_report.py 重建树与报告。"
        )

    for name in ROUTES[route]:
        if name == "report":
            continue
        context.data.update(state.snapshot["stages"].get(name, {}).get("context_data", {}))

    stage_dir = run_dir / f"step{len(ROUTES[route]):02d}_report"
    stage_dir.mkdir(parents=True, exist_ok=True)
    context.stage_dir = stage_dir
    module = importlib.import_module("cdc_pipeline.stages.report")
    signature = stage_signature(context, "report", module.__file__)
    state.begin_stage("report", signature)
    result = module.run(context)
    state.complete_stage(
        "report",
        signature,
        result.artifacts,
        status=result.status,
        summary=result.summary,
        context_data=result.data,
    )
    state.finish("done")
    return result.status


def main(argv: list[str]) -> int:
    if not argv[1:]:
        print("用法: python scripts/rerun_report_only.py <run_dir> ...", file=sys.stderr)
        return 2
    for item in argv[1:]:
        run_dir = Path(item).resolve()
        if not (run_dir / "state.json").is_file():
            print(f"跳过（无 state.json）: {run_dir}")
            continue
        try:
            status = rerun(run_dir)
            print(f"{run_dir.name}: report {status}")
        except Exception as exc:
            print(f"{run_dir.name}: FAILED {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
