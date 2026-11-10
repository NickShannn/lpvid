"""reads 质控、去 rRNA 和 MEGAHIT 组装。"""

from __future__ import annotations

import gzip
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from cdc_pipeline.fasta import read_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path
from cdc_pipeline.state import atomic_write_json


@dataclass
class PreprocessResult:
    contigs_fna: Path
    clean_r1: Path
    clean_r2: Path | None
    clean_read_count: int
    qc_summary: dict[str, Any]
    assembly_stats: dict[str, Any]


def fastp_args(r1: str, r2: str | None, config: dict, workdir: Path) -> list[str]:
    args = [
        "fastp", "-i", str(r1), "-o", str(workdir / "trimmed_r1.fastq.gz"),
        "--json", str(workdir / "qc.json"), "--html", str(workdir / "qc.html"),
        "--thread", str(config["qc"]["threads"]), "--detect_adapter_for_pe",
    ]
    if r2:
        args += ["-I", str(r2), "-O", str(workdir / "trimmed_r2.fastq.gz")]
    return args


def rrna_args(index: str, r1: str, r2: str | None, workdir: Path, tag: str, threads: int = 16) -> list[str]:
    args = ["bowtie2", "-x", str(index), "-p", str(threads), "-S", "/dev/null", "--quiet"]
    if r2:
        args += [
            "--un-conc-gz", str(workdir / f"{tag}_keep_%.fastq.gz"),
            "-1", str(r1), "-2", str(r2),
        ]
    else:
        args += ["--un-gz", str(workdir / f"{tag}_keep_1.fastq.gz"), "-U", str(r1)]
    return args


def megahit_args(reads: list[str], config: dict, output_dir: Path) -> list[str]:
    assembly = config["assembly"]
    args = [
        "megahit", "-t", str(assembly["threads"]),
        "--min-contig-len", str(assembly["min_contig_len"]),
        "-m", str(int(assembly["memory_gb"]) * 1024 ** 3),
        "--force", "-o", str(output_dir),
    ]
    if len(reads) == 2:
        args += ["-1", reads[0], "-2", reads[1]]
    elif len(reads) == 1:
        args += ["-r", reads[0]]
    else:
        raise ValueError("MEGAHIT 需要一个或两个 FASTQ 输入")
    return args


def count_fastq_reads(path: str | Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"FASTQ 行数不是 4 的倍数: {path}")
    return lines // 4


def _assembly_stats(path: Path) -> dict[str, Any]:
    lengths = sorted((record.length for record in read_fasta(path)), reverse=True)
    total = sum(lengths)
    accumulated = 0
    n50 = 0
    for length in lengths:
        accumulated += length
        if accumulated >= total / 2:
            n50 = length
            break
    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "n50": n50,
        "max_len": lengths[0] if lengths else 0,
    }


def run_preprocess(inputs: dict[str, str], config: dict, workdir: Path) -> PreprocessResult:
    workdir.mkdir(parents=True, exist_ok=True)
    r1 = inputs["r1"]
    r2 = inputs.get("r2") or None
    run_command(fastp_args(r1, r2, config, workdir), env_name="qc", cwd=workdir, log_path=workdir / "fastp.log")
    current1 = workdir / "trimmed_r1.fastq.gz"
    current2 = workdir / "trimmed_r2.fastq.gz" if r2 else None
    rrna_summary = []
    if config["qc"]["rrna_deplete"]:
        for tag, key in [("ssu", "rrna_ssu"), ("lsu", "rrna_lsu")]:
            before = count_fastq_reads(current1)
            run_command(
                rrna_args(config["assets"][key], current1, current2, workdir, tag, config["qc"]["threads"]),
                env_name="qc", cwd=workdir, log_path=workdir / f"rrna_{tag}.log",
            )
            next1 = workdir / f"{tag}_keep_1.fastq.gz"
            next2 = workdir / f"{tag}_keep_2.fastq.gz" if r2 else None
            if not next1.exists() or (r2 and not next2.exists()):
                raise FileNotFoundError(f"bowtie2 {tag} 未产生预期 clean reads")
            after = count_fastq_reads(next1)
            rrna_summary.append({"reference": tag, "reads_or_pairs_before": before, "reads_or_pairs_after": after, "removed": before - after})
            current1, current2 = next1, next2
    clean1 = workdir / "clean_r1.fastq.gz"
    clean2 = workdir / "clean_r2.fastq.gz" if r2 else None
    os.replace(current1, clean1)
    if r2:
        os.replace(current2, clean2)
    clean_count = count_fastq_reads(clean1) + (count_fastq_reads(clean2) if clean2 else 0)
    qc_json = json.loads((workdir / "qc.json").read_text())
    qc_summary = {
        "before_filtering": qc_json.get("summary", {}).get("before_filtering", {}),
        "after_filtering": qc_json.get("summary", {}).get("after_filtering", {}),
        "rrna": rrna_summary,
        "clean_read_count_for_rpm": clean_count,
        "count_unit": "reads; PE counts R1+R2",
    }
    megahit_out = workdir / "megahit_out"
    reads = [str(clean1)] + ([str(clean2)] if clean2 else [])
    run_command(megahit_args(reads, config, megahit_out), env_name="qc", cwd=workdir, log_path=workdir / "megahit.log")
    assembled = megahit_out / "final.contigs.fa"
    if not assembled.is_file():
        raise FileNotFoundError(f"MEGAHIT 未产出 {assembled}")
    contigs = workdir / "contigs.fna"
    shutil.copyfile(assembled, contigs)
    stats = _assembly_stats(contigs)
    atomic_write_json(workdir / "assembly_stats.json", stats)
    atomic_write_json(workdir / "qc_summary.json", qc_summary)
    return PreprocessResult(contigs, clean1, clean2, clean_count, qc_summary, stats)


def run(context: PipelineContext) -> StageResult:
    result = run_preprocess(context.inputs, context.config, context.stage_dir)
    artifacts = {
        "contigs_fna": artifact_path(context, result.contigs_fna),
        "clean_r1": artifact_path(context, result.clean_r1),
        "clean_r2": artifact_path(context, result.clean_r2) if result.clean_r2 else None,
        "qc_html": artifact_path(context, context.stage_dir / "qc.html"),
        "qc_summary": artifact_path(context, context.stage_dir / "qc_summary.json"),
        "assembly_stats": artifact_path(context, context.stage_dir / "assembly_stats.json"),
    }
    return StageResult(
        artifacts=artifacts,
        summary={"qc": result.qc_summary, "assembly": result.assembly_stats},
        data={
            "contigs_fna": str(result.contigs_fna),
            "clean_r1": str(result.clean_r1),
            "clean_r2": str(result.clean_r2) if result.clean_r2 else None,
            "clean_read_count": result.clean_read_count,
            "qc_summary": result.qc_summary,
            "assembly_stats": result.assembly_stats,
        },
    )
