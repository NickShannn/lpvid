"""RdRP 输入筛查、模型桥接和推理适配器。"""

from __future__ import annotations

import csv
import io
import os
import shutil
from pathlib import Path

from cdc_pipeline.fasta import FastaRecord, read_fasta
from cdc_pipeline.prediction import ModelPrediction
from cdc_pipeline.process import run_command

TOXIC_CHARS = set("ATCGUN")
TIME_STR = "protv2full_20261009100409"
STEP = "95936"


def is_toxic_sequence(sequence: str) -> bool:
    sequence = sequence.upper()
    total = sum(char.isalpha() for char in sequence)
    atcgu = sum(char in TOXIC_CHARS for char in sequence)
    return total > 0 and atcgu >= 0.9 * total and not sequence.startswith("M")


def screen_rdrp_input(records: list[FastaRecord]) -> tuple[list[FastaRecord], list[FastaRecord]]:
    clean, excluded = [], []
    for record in records:
        (excluded if is_toxic_sequence(record.sequence) else clean).append(record)
    return clean, excluded


def parse_rdrp_text(text: str) -> dict[str, ModelPrediction]:
    result = {}
    for row in csv.DictReader(io.StringIO(text)):
        if not row or not row.get("seq_id"):
            continue
        query_id = row["seq_id"].strip().lstrip(">").split()[0]
        result[query_id] = ModelPrediction(query_id, float(row["prob"]), str(row.get("label", "")))
    return result


def parse_rdrp(path: str | Path) -> dict[str, ModelPrediction]:
    return parse_rdrp_text(Path(path).read_text(encoding="utf-8"))


def prepare_model_bridge(repo: Path, destination: Path) -> Path:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    for source in repo.iterdir():
        if source.name in {"models", "logs"}:
            continue
        target = destination / source.name
        target.symlink_to(source.resolve(), target_is_directory=source.is_dir())
    protein = Path("RdRP/protein")
    checkpoint_source = repo / "models" / protein / "binary_class_v2/luca_base/matrix" / TIME_STR / f"checkpoint-{STEP}"
    checkpoint_target = destination / "models" / protein / "binary_class/luca_base/matrix" / TIME_STR / f"checkpoint-{STEP}"
    checkpoint_target.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_target.symlink_to(checkpoint_source.resolve(), target_is_directory=True)
    v2_model_target = destination / "models" / protein / "binary_class_v2"
    v2_model_target.parent.mkdir(parents=True, exist_ok=True)
    v2_model_target.symlink_to((repo / "models" / protein / "binary_class_v2").resolve(), target_is_directory=True)
    label_source = checkpoint_source.parent / "best" / "label.txt"
    for task_type in ["binary_class", "binary_class_v2"]:
        label_target = destination / "logs" / protein / task_type / "luca_base/matrix" / TIME_STR / "label.txt"
        label_target.parent.mkdir(parents=True, exist_ok=True)
        label_target.symlink_to(label_source.resolve())
    return destination


def run_rdrp(input_faa: str | Path, config: dict, workdir: str | Path, gpu: str | int) -> Path:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    expected = len(list(read_fasta(input_faa)))
    output = workdir / "rdrp_predictions.csv"
    if expected == 0:
        output.write_text("seq_id,seq,prob,label_index,label\n")
        return output
    repo = Path(config["assets"]["rdrp_model_dir"])
    bridge = prepare_model_bridge(repo, workdir / "model_bridge")
    raw = workdir / "preds.csv"
    if raw.exists():
        raw.unlink()
    compat = str(config["assets"]["runtime_lib_dir"])
    ld_library = compat + (os.pathsep + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")
    run_command(
        [
            config["environments"]["rdrp"], "predict.py",
            "--input_file", Path(input_faa).resolve(), "--seq_type", "prot",
            "--model_path", bridge, "--dataset_name", "RdRP", "--dataset_type", "protein",
            "--task_type", "binary_class", "--task_level_type", "seq_level",
            "--model_type", "luca_base", "--input_type", "matrix", "--input_mode", "single",
            "--time_str", TIME_STR, "--step", STEP, "--save_path", raw,
            "--llm_truncation_seq_length", config["detect"]["rdrp_truncation"],
            "--print_per_num", "5000", "--gpu_id", "0", "--threshold", config["model_threshold"],
        ],
        env_name="rdrp", cwd=repo / "src", log_path=workdir / "rdrp.log",
        env={
            "CUDA_VISIBLE_DEVICES": str(gpu),
            "LD_LIBRARY_PATH": ld_library,
            "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD": "1",
            "OMP_NUM_THREADS": "8",
        },
    )
    if not raw.is_file():
        raise FileNotFoundError(f"RdRP 识别模型未产出 {raw}")
    shutil.copyfile(raw, output)
    rows = parse_rdrp(output)
    if len(rows) != expected:
        raise ValueError(f"RdRP 输入/结果行数不一致: {expected}/{len(rows)}")
    return output
