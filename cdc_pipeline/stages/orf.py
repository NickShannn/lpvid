"""fastorf 调用和 ORF→contig→稳定 query ID 映射。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.fasta import assign_stable_ids, read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path, write_tsv


@dataclass(frozen=True)
class OrfMapRow:
    orf_id: str
    contig_id: str
    description: str
    strand: str
    frame: int
    start: int
    end: int
    nt_len: int
    aa_len: int
    start_codon: str
    stop_codon: str
    status: str


@dataclass
class OrfResult:
    proteins_faa: Path
    raw_faa: Path
    cds_fnt: Path
    gff3: Path
    raw_map: Path
    orf_contig_map: Path
    n_orfs: int


def parse_orf_map_row(parts: list[str]) -> OrfMapRow:
    if len(parts) < 12:
        raise ValueError(f"fastorf TSV 需要至少 12 列，实际 {len(parts)}")
    return OrfMapRow(
        orf_id=parts[0], contig_id=parts[1], description=parts[2], strand=parts[3],
        frame=int(parts[4]), start=int(parts[5]), end=int(parts[6]), nt_len=int(parts[7]),
        aa_len=int(parts[8]), start_codon=parts[9], stop_codon=parts[10], status=parts[11],
    )


def read_orf_map(path: str | Path) -> list[OrfMapRow]:
    rows = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if parts[0] == "orf_id":
                continue
            rows.append(parse_orf_map_row(parts))
    return rows


def normalize_orf_outputs(raw_faa: Path, raw_map: Path, sample: str, workdir: Path) -> OrfResult:
    records = list(read_fasta(raw_faa))
    rows = read_orf_map(raw_map)
    if len(records) != len(rows):
        raise ValueError(f"fastorf FAA({len(records)}) 与 TSV({len(rows)}) 行数不一致")
    for record, row in zip(records, rows):
        if record.record_id != row.orf_id:
            raise ValueError(f"fastorf FAA/TSV ID 错位: {record.record_id} != {row.orf_id}")
    assigned, id_mapping = assign_stable_ids(records, sample)
    proteins = workdir / "proteins.faa"
    write_fasta(proteins, assigned)
    unified = workdir / "orf_contig_map.tsv"
    fields = [
        "sample", "query_id", "orf_id", "contig_id", "original_description", "strand", "frame",
        "start", "end", "nt_len", "aa_len", "start_codon", "stop_codon", "status",
    ]
    output_rows = []
    for id_row, orf_row in zip(id_mapping, rows):
        output_rows.append({
            "sample": sample,
            "query_id": id_row["query_id"],
            "orf_id": orf_row.orf_id,
            "contig_id": orf_row.contig_id,
            "original_description": id_row["original_description"],
            "strand": orf_row.strand,
            "frame": orf_row.frame,
            "start": orf_row.start,
            "end": orf_row.end,
            "nt_len": orf_row.nt_len,
            "aa_len": orf_row.aa_len,
            "start_codon": orf_row.start_codon,
            "stop_codon": orf_row.stop_codon,
            "status": orf_row.status,
        })
    write_tsv(unified, fields, output_rows)
    return OrfResult(
        proteins_faa=proteins,
        raw_faa=raw_faa,
        cds_fnt=workdir / "orfs.cds.fnt",
        gff3=workdir / "orfs.gff3",
        raw_map=raw_map,
        orf_contig_map=unified,
        n_orfs=len(records),
    )


def run_fastorf(contigs: str | Path, sample: str, config: dict, workdir: Path) -> OrfResult:
    workdir.mkdir(parents=True, exist_ok=True)
    raw_faa = workdir / "orfs.raw.faa"
    raw_map = workdir / "fastorf_map.tsv"
    run_command(
        [
            config["tools"]["fastorf"], "--in", contigs,
            "--faa", raw_faa, "--fnt", workdir / "orfs.cds.fnt",
            "--tsv", raw_map, "--gff3", workdir / "orfs.gff3",
            "--threads", config["orf"]["threads"],
        ],
        env_name="lucaphylo", cwd=workdir, log_path=workdir / "fastorf.log",
    )
    return normalize_orf_outputs(raw_faa, raw_map, sample, workdir)


def run(context: PipelineContext) -> StageResult:
    contigs = context.data.get("contigs_fna") or context.inputs["input"]
    result = run_fastorf(contigs, context.sample, context.config, context.stage_dir)
    return StageResult(
        artifacts={
            "proteins_faa": artifact_path(context, result.proteins_faa),
            "raw_faa": artifact_path(context, result.raw_faa),
            "cds_fnt": artifact_path(context, result.cds_fnt),
            "gff3": artifact_path(context, result.gff3),
            "raw_map": artifact_path(context, result.raw_map),
            "orf_contig_map": artifact_path(context, result.orf_contig_map),
        },
        summary={"n_orfs": result.n_orfs},
        data={
            "contigs_fna": str(contigs),
            "proteins_faa": str(result.proteins_faa),
            "orf_contig_map": str(result.orf_contig_map),
            "n_orfs": result.n_orfs,
        },
    )
