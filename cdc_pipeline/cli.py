"""lpvid 的 lpvid 命令行入口。"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from cdc_pipeline.config import load_config
from cdc_pipeline.pipeline import resume_pipeline, run_pipeline
from cdc_pipeline.project import project_metadata
from cdc_pipeline.state import RunState


def _config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", help="本机 YAML 配置；默认合并 config/local.yaml")


def build_parser() -> argparse.ArgumentParser:
    metadata = project_metadata()
    parser = argparse.ArgumentParser(
        prog=metadata["command"],
        description=f"{metadata['display_name']}: {metadata['description']}",
    )
    parser.add_argument("--version", action="version", version=metadata["version"])
    sub = parser.add_subparsers(dest="command", required=True)

    doctor = sub.add_parser("doctor", help="检查环境、工具和资产")
    _config_arg(doctor)

    prepare = sub.add_parser("prepare", help="同步并构建 24 超群参考资产")
    _config_arg(prepare)
    prepare.add_argument("--prewarm", action="store_true", help="预留骨架树特征预热")

    for route in ("reads", "contigs", "proteins"):
        run = sub.add_parser(f"run-{route}", help=f"运行 {route} 输入路线")
        _config_arg(run)
        if route == "reads":
            run.add_argument("--r1", required=True)
            run.add_argument("--r2")
        else:
            run.add_argument("--input", required=True)
        run.add_argument("--out", required=True)
        run.add_argument("--sample", default="sample")
        run.add_argument("--gpu-id", default="auto")

    resume = sub.add_parser("resume", help="恢复已有运行")
    resume.add_argument("run_dir")
    _config_arg(resume)
    resume.add_argument("--gpu-id")

    status = sub.add_parser("status", help="查看运行状态")
    status.add_argument("run_dir")
    status.add_argument("--json", action="store_true")
    return parser


def doctor(config: dict) -> tuple[bool, list[dict[str, str]]]:
    checks: list[dict[str, str]] = []
    for name, path in config["environments"].items():
        checks.append({"name": f"env:{name}", "path": path, "status": "ok" if Path(path).is_file() else "missing"})
    required = {
        "fastorf": config["tools"]["fastorf"],
        "diamond": config["tools"]["diamond"],
        "cd_hit": config["tools"]["cd_hit"],
        "fastp": str(Path(config["environments"]["qc"]).parent / "fastp"),
        "bowtie2": str(Path(config["environments"]["qc"]).parent / "bowtie2"),
        "megahit": str(Path(config["environments"]["qc"]).parent / "megahit"),
        "samtools": str(Path(config["environments"]["qc"]).parent / "samtools"),
        "refseq_db": config["assets"]["refseq_db"],
        "viral_taxids": config["assets"]["viral_taxids"],
        "cellular_taxids": config["assets"]["cellular_taxids"],
        "phylo_base": config["assets"]["phylo_base_model"] + "/config.json",
        "rdrp_model": config["assets"]["rdrp_model_dir"] + "/models/RdRP/protein/binary_class_v2/luca_base/matrix/protv2full_20261009100409/checkpoint-95936/model.safetensors",
        "rdrp_label": config["assets"]["rdrp_model_dir"] + "/models/RdRP/protein/binary_class_v2/luca_base/matrix/protv2full_20261009100409/best/label.txt",
        "rdrp_libtiff": config["assets"]["runtime_lib_dir"] + "/libtiff.so.6",
        "rdrp_libjpeg": config["assets"]["runtime_lib_dir"] + "/libjpeg.so.8",
        "supergroup24_db": config["paths"]["assets_dir"] + "/supergroup24/reference.dmnd",
        "lucaphylo_inference": config["paths"]["third_party_dir"] + "/lucaphylo/infer_stage_4.py",
        "tree_template": config["paths"]["third_party_dir"] + "/placement_view/tree_explorer_template.html",
    }
    if config["qc"]["rrna_deplete"]:
        required["rrna_ssu"] = config["assets"]["rrna_ssu"] + ".1.bt2"
        required["rrna_lsu"] = config["assets"]["rrna_lsu"] + ".1.bt2"
    for name, path in required.items():
        checks.append({"name": name, "path": path, "status": "ok" if Path(path).exists() else "missing"})
    from cdc_pipeline.assets import load_supergroup_map
    asset_count = len(load_supergroup_map(config))
    checks.append({"name": "phylo_supergroups", "path": config["assets"]["phylo_root"], "status": "ok" if asset_count >= 24 else "missing"})
    return all(row["status"] == "ok" for row in checks), checks


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "status":
            snapshot = RunState(args.run_dir).snapshot
            if args.json:
                print(json.dumps(snapshot, ensure_ascii=False, indent=2))
            else:
                print(f"status: {snapshot.get('status')}")
                for name, stage in snapshot.get("stages", {}).items():
                    print(f"{name}\t{stage.get('status')}\t{stage.get('error') or ''}")
            return 0
        if args.command == "doctor":
            ok, checks = doctor(load_config(args.config))
            for row in checks:
                print(f"[{row['status'].upper()}] {row['name']}\t{row['path']}")
            return 0 if ok else 1
        if args.command == "prepare":
            from cdc_pipeline.assets import prepare_assets
            result = prepare_assets(load_config(args.config), prewarm=args.prewarm)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0
        if args.command.startswith("run-"):
            route = args.command.removeprefix("run-")
            inputs = {"r1": args.r1, "r2": args.r2} if route == "reads" else {"input": args.input}
            summary = run_pipeline(route, inputs, args.out, load_config(args.config), sample=args.sample, gpu=args.gpu_id)
            print(f"{summary.status}\t{summary.run_dir}")
            return 0
        if args.command == "resume":
            config = load_config(args.config) if args.config else None
            summary = resume_pipeline(args.run_dir, config=config, gpu=args.gpu_id)
            print(f"{summary.status}\t{summary.run_dir}")
            return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
