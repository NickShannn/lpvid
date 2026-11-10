"""路线定义和通用阶段执行器。"""

from __future__ import annotations

import hashlib
import importlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cdc_pipeline.project import project_metadata
from cdc_pipeline.state import RunState, compute_signature, file_fingerprint


def code_fingerprint(path: str | Path) -> dict[str, Any]:
    """阶段代码指纹：只依赖文件内容（sha256），不含 mtime。

    mtime 参与签名会带来误判：编辑代码后即使回滚到完全相同的内容，
    也会让所有既有运行重跑。数据文件仍用 file_fingerprint（含 mtime），
    因为数据变化需要用廉价方式立刻感知。
    """
    target = Path(path)
    digest = hashlib.sha256()
    with target.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(target.resolve()), "sha256": digest.hexdigest()}

ROUTES = {
    "proteins": ["prefilter", "detect", "representatives", "supergroup", "classifier", "phylo", "visualization", "report"],
    "contigs": ["orf", "prefilter", "detect", "representatives", "supergroup", "classifier", "phylo", "visualization", "report"],
    "reads": ["preprocess", "orf", "prefilter", "detect", "representatives", "supergroup", "classifier", "phylo", "abundance", "visualization", "report"],
}


@dataclass
class StageResult:
    artifacts: dict[str, Any] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)
    data: dict[str, Any] = field(default_factory=dict)
    status: str = "done"


@dataclass
class PipelineContext:
    route: str
    sample: str
    inputs: dict[str, Any]
    config: dict[str, Any]
    run_dir: Path
    state: RunState
    gpu: str = "auto"
    data: dict[str, Any] = field(default_factory=dict)
    stage_dir: Path | None = None


@dataclass
class RunSummary:
    run_dir: Path
    status: str
    stages: list[str]


def stages_for(route: str) -> list[str]:
    if route not in ROUTES:
        raise ValueError(f"未知输入路线: {route}")
    return list(ROUTES[route])


def _fingerprint_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _fingerprint_value(child) for key, child in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_fingerprint_value(child) for child in value]
    if isinstance(value, Path) or (isinstance(value, str) and len(value) < 4096 and Path(value).exists()):
        path = Path(value)
        if path.is_file():
            return file_fingerprint(path)
        return str(path)
    return value


def _taxonomy_dependencies() -> list[Path]:
    root = Path(__file__).parent.parent
    files: list[Path] = []
    for taxonomy in (root / "resources" / "phylo" / "taxonomy", root / "assets" / "phylo" / "taxonomy"):
        if taxonomy.is_dir():
            files.extend(sorted(taxonomy.glob("*.csv")))
    return files


def stage_dependencies(name: str) -> list[Path]:
    if name == "report":
        return [
            Path(__file__).parent / "html_report.py",
            Path(__file__).parent / "html_template_v2.py",
            Path(__file__).parent / "risk_assessment.py",

            Path(__file__).parent.parent / "scripts" / "verify_tree_explorer.py",
            *_taxonomy_dependencies(),
            Path(__file__).parent.parent / "assets" / "risk" / "china_pathogen_catalog.csv",
        ]
    if name == "visualization":
        root = Path(__file__).parent.parent
        view = root / "third_party" / "placement_view"
        return [
            Path(__file__).parent / "visualization.py",
            Path(__file__).parent / "risk_assessment.py",
            *_taxonomy_dependencies(),
            view / "tree_explorer_template.html",
            view / "phylo_core.js",
            view / "vendor" / "jspdf.umd.min.js",
            view / "vendor" / "svg2pdf.umd.min.js",
        ]
    return []


def stage_signature(context: PipelineContext, name: str, module_path: str) -> str:
    payload = {
        "stage": name,
        "route": context.route,
        "sample": context.sample,
        "inputs": _fingerprint_value(context.inputs),
        "data": _fingerprint_value(context.data),
        "config": {key: value for key, value in context.config.items() if key != "project"},
        "code": code_fingerprint(module_path),
    }
    dependencies = stage_dependencies(name)
    if dependencies:
        payload["dependencies"] = [code_fingerprint(path) for path in dependencies]
    return compute_signature(payload)


def run_pipeline(
    route: str,
    inputs: dict[str, Any],
    output_dir: str | Path,
    config: dict[str, Any],
    *,
    sample: str = "sample",
    gpu: str = "auto",
    resume: bool = False,
) -> RunSummary:
    stages = stages_for(route)


    run_dir = Path(output_dir).resolve()
    inputs = {
        key: str(Path(value).resolve()) if isinstance(value, str) and Path(value).exists() else value
        for key, value in inputs.items()
    }
    if not resume and run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"输出目录非空，请换新目录或显式 resume: {run_dir}")
    state = RunState(run_dir)
    metadata = {"route": route, "sample": sample, "inputs": inputs, "gpu": gpu}
    if resume:
        original = state.manifest
        if original["route"] != route:
            raise ValueError(f"恢复路线不一致: {original['route']} != {route}")
    else:
        state.start(metadata, config)
    context = PipelineContext(route, sample, inputs, config, run_dir, state, gpu)
    for index, name in enumerate(stages):
        module = importlib.import_module(f"cdc_pipeline.stages.{name}")
        signature = stage_signature(context, name, module.__file__)
        if state.should_skip(name, signature):
            print(f"[{index + 1}/{len(stages)}] {name}: SKIP (signature and artifacts unchanged)", flush=True)
            context.data.update(state.snapshot["stages"][name].get("context_data", {}))
            continue
        print(f"[{index + 1}/{len(stages)}] {name}: START", flush=True)
        state.invalidate_stages(stages[index + 1:], reason=f"upstream_changed:{name}")
        stage_dir = run_dir / f"step{index + 1:02d}_{name}"
        stage_dir.mkdir(parents=True, exist_ok=True)
        context.stage_dir = stage_dir
        state.begin_stage(name, signature)
        try:
            result = module.run(context)
            state.complete_stage(
                name,
                signature,
                result.artifacts,
                status=result.status,
                summary=result.summary,
                context_data=result.data,
            )
            context.data.update(result.data)
            print(f"[{index + 1}/{len(stages)}] {name}: {result.status.upper()}", flush=True)
        except Exception as exc:
            state.fail_stage(name, signature, str(exc))
            print(f"[{index + 1}/{len(stages)}] {name}: FAILED: {exc}", flush=True)
            raise
    state.finish("done")
    return RunSummary(run_dir, "done", stages)


def resume_pipeline(run_dir: str | Path, config: dict[str, Any] | None = None, gpu: str | None = None) -> RunSummary:
    state = RunState(run_dir)
    manifest = state.manifest
    runtime_config = config or manifest["config"]
    runtime_config.setdefault("project", project_metadata())
    return run_pipeline(
        manifest["route"],
        manifest["inputs"],
        run_dir,
        runtime_config,
        sample=manifest.get("sample", "sample"),
        gpu=gpu or manifest.get("gpu", "auto"),
        resume=True,
    )
