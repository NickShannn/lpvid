#!/usr/bin/env python3
"""在不重跑数据的前提下刷新已完成阶段的签名。

背景：阶段签名包含模块文件的 mtime；代码编辑后即使内容回滚，签名也会变化，
导致已有运行被误判为需要重跑。本脚本用当前代码 + 当前产物重算签名，
只对状态为 done/skipped/not_applicable 且产物完整的阶段生效。

用法：
    python scripts/refresh_stage_signatures.py runs/real_data/<样本> [...]
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from cdc_pipeline.pipeline import ROUTES, PipelineContext, stage_signature
from cdc_pipeline.project import project_metadata
from cdc_pipeline.state import RunState, TERMINAL_SKIP_STATUSES


def _artifacts_complete(state: RunState, artifacts: dict) -> bool:
    for value in (artifacts or {}).values():
        if not isinstance(value, str):
            continue
        path = Path(value)
        if not path.is_absolute():
            path = state.run_dir / path
        if value and not path.exists():
            return False
    return True


def refresh(run_dir: Path, restore: bool = False) -> tuple[list[str], list[str]]:
    state = RunState(run_dir)
    manifest = state.manifest
    config = dict(manifest["config"])
    config.setdefault("project", project_metadata())
    context = PipelineContext(
        route=manifest["route"],
        sample=manifest.get("sample", "sample"),
        inputs=manifest.get("inputs", {}),
        config=config,
        run_dir=run_dir,
        state=state,
        gpu=manifest.get("gpu", "auto"),
    )
    checkpoint = {}
    if restore and state.checkpoint_path.exists():
        checkpoint = json.loads(state.checkpoint_path.read_text(encoding="utf-8")).get("stages", {})

    changed: list[str] = []
    restored: list[str] = []
    for name in ROUTES[manifest["route"]]:
        stage = state.snapshot["stages"].get(name)
        if not stage:
            continue
        if stage.get("status") not in TERMINAL_SKIP_STATUSES:
            saved = checkpoint.get(name) if restore else None
            if saved and saved.get("status") in TERMINAL_SKIP_STATUSES and _artifacts_complete(state, saved.get("artifacts", {})):
                state._data["stages"][name] = saved
                stage = saved
                restored.append(name)
            else:
                context.data.update(stage.get("context_data", {}))
                continue
        context.data.update(stage.get("context_data", {}))
        module = importlib.import_module(f"cdc_pipeline.stages.{name}")
        signature = stage_signature(context, name, module.__file__)
        if signature != stage.get("signature"):
            state._data["stages"][name]["signature"] = signature
            changed.append(name)
    if changed or restored:
        state._flush()
    return changed, restored


def force_align(run_dir: Path) -> list[str]:
    """把有完整产物的阶段强制对齐到当前签名（状态 done），report/visualization 留给 resume 重建。"""
    state = RunState(run_dir)
    manifest = state.manifest
    config = dict(manifest["config"])
    config.setdefault("project", project_metadata())
    context = PipelineContext(
        route=manifest["route"],
        sample=manifest.get("sample", "sample"),
        inputs=manifest.get("inputs", {}),
        config=config,
        run_dir=run_dir,
        state=state,
        gpu=manifest.get("gpu", "auto"),
    )
    aligned: list[str] = []
    for name in ROUTES[manifest["route"]]:

        if name in {"report", "visualization"}:
            state._data["stages"].pop(name, None)
            continue
        stage = state._data["stages"].get(name)
        if not stage:
            continue
        if not _artifacts_complete(state, stage.get("artifacts", {})):
            context.data.update(stage.get("context_data", {}))
            continue
        context.data.update(stage.get("context_data", {}))
        module = importlib.import_module(f"cdc_pipeline.stages.{name}")
        signature = stage_signature(context, name, module.__file__)
        if stage.get("status") not in TERMINAL_SKIP_STATUSES:
            stage["status"] = "done"
            stage["error"] = None
        if stage.get("signature") != signature:
            stage["signature"] = signature
            aligned.append(name)
    state._flush()
    return aligned


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    restore = "--restore" in argv[1:]
    align = "--align" in argv[1:]
    if not args:
        print("用法: python scripts/refresh_stage_signatures.py [--restore|--align] <run_dir> ...", file=sys.stderr)
        return 2
    for item in args:
        run_dir = Path(item).resolve()
        if not (run_dir / "run_manifest.json").is_file():
            print(f"跳过（无 run_manifest.json）: {run_dir}")
            continue
        if align:
            aligned = force_align(run_dir)
            print(f"{run_dir.name}: 对齐 {aligned or '无变化'}（visualization/report 已清除，将重建）")
        else:
            changed, restored = refresh(run_dir, restore=restore)
            print(f"{run_dir.name}: 恢复 {restored or '无'} | 更新签名 {changed or '无'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
