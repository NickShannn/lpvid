"""clean reads 回比阳性 contigs，计算 RPM、深度和覆盖度。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.fasta import read_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path, read_tsv, write_tsv
from cdc_pipeline.stages.preprocess import count_fastq_reads
from cdc_pipeline.state import atomic_write_json


@dataclass
class AbundanceResult:
    table: Path
    supergroup_table: Path
    bam: Path | None
    bai: Path | None
    summary: Path
    counts: dict[str, int]
    status: str


def calculate_rpm(mapped_reads: int, total_clean_reads: int) -> float:
    return mapped_reads / total_clean_reads * 1_000_000 if total_clean_reads else 0.0


def clean_read_denominator(r1: int, r2: int = 0) -> int:
    return r1 + r2


def build_abundance_rows(records, mapped: dict[str, int], coverage: dict[str, dict], total_clean_reads: int) -> list[dict]:
    rows = []
    for record in records:
        contig = record.record_id
        stats = coverage.get(contig, {})
        rows.append({
            "contig_id": contig,
            "length": record.length,
            "mapped_reads": mapped.get(contig, 0),
            "RPM": calculate_rpm(mapped.get(contig, 0), total_clean_reads),
            "mean_depth": float(stats.get("meandepth", 0.0)),
            "coverage_breadth_pct": float(stats.get("coverage", 0.0)),
            "total_clean_reads": total_clean_reads,
            "count_unit": "primary_mapped_reads",
        })
    return rows


def aggregate_supergroups(rows: list[dict], contig_supergroups: dict[str, set[str]], total_clean_reads: int) -> list[dict]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        groups = contig_supergroups.get(row["contig_id"], set())
        if len(groups) == 1:
            grouped.setdefault(next(iter(groups)), []).append(row)
    output = []
    for supergroup, members in sorted(grouped.items()):
        length = sum(int(row["length"]) for row in members)
        mapped = sum(int(row["mapped_reads"]) for row in members)
        output.append({
            "supergroup": supergroup,
            "n_contigs": len(members),
            "mapped_reads": mapped,
            "RPM": calculate_rpm(mapped, total_clean_reads),
            "mean_depth_weighted": sum(float(row["mean_depth"]) * int(row["length"]) for row in members) / length if length else 0.0,
            "coverage_breadth_pct": sum(float(row["coverage_breadth_pct"]) * int(row["length"]) for row in members) / length if length else 0.0,
            "total_clean_reads": total_clean_reads,
        })
    return output


def _read_idxstats(path: Path) -> dict[str, int]:
    mapped = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                mapped[parts[0]] = int(parts[2])
    return mapped


def _read_coverage(path: Path) -> dict[str, dict]:
    with path.open(newline="") as handle:
        return {
            row["#rname"]: {"meandepth": float(row["meandepth"]), "coverage": float(row["coverage"])}
            for row in csv.DictReader(handle, delimiter="\t")
        }


def run_abundance(
    clean_reads: dict[str, str],
    positive_contigs: str | Path,
    config: dict,
    workdir: str | Path,
    total_clean_reads: int | None = None,
    contig_supergroups: dict[str, set[str]] | None = None,
) -> AbundanceResult:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    records = list(read_fasta(positive_contigs)) if positive_contigs else []
    total = total_clean_reads if total_clean_reads is not None else clean_read_denominator(
        count_fastq_reads(clean_reads["r1"]), count_fastq_reads(clean_reads["r2"]) if clean_reads.get("r2") else 0
    )
    table = workdir / "abundance.tsv"
    supergroup_table = workdir / "supergroup_abundance.tsv"
    summary = workdir / "abundance_summary.json"
    fields = ["contig_id", "length", "mapped_reads", "RPM", "mean_depth", "coverage_breadth_pct", "total_clean_reads", "count_unit"]
    if not records:
        write_tsv(table, fields, [])
        write_tsv(supergroup_table, ["supergroup", "n_contigs", "mapped_reads", "RPM", "mean_depth_weighted", "coverage_breadth_pct", "total_clean_reads"], [])
        write_tsv(workdir / "abundance_contig_assignments.tsv", ["contig_id", "supergroups", "status"], [])
        counts = {"contigs": 0, "mapped_reads": 0, "total_clean_reads": total, "ambiguous_contigs": 0}
        atomic_write_json(summary, {"status": "skipped_no_positive_contigs", "counts": counts})
        return AbundanceResult(table, supergroup_table, None, None, summary, counts, "skipped")

    index = workdir / "contigs_bt2"
    run_command(["bowtie2-build", "--threads", config["abundance"]["threads"], positive_contigs, index], env_name="qc", cwd=workdir, log_path=workdir / "bowtie2_build.log")
    sam = workdir / "mapped.sam"
    command = [config["abundance"]["bowtie2_mode"], "-x", index, "--threads", config["abundance"]["threads"], "--no-unal", "-S", sam]
    if clean_reads.get("r2"):
        command += ["-1", clean_reads["r1"], "-2", clean_reads["r2"]]
    else:
        command += ["-U", clean_reads["r1"]]
    run_command(["bowtie2", *command], env_name="qc", cwd=workdir, log_path=workdir / "bowtie2.log")
    unsorted_bam = workdir / "mapped.unsorted.bam"
    all_bam = workdir / "mapped.all.bam"
    bam = workdir / "mapped.bam"
    run_command(["samtools", "view", "-b", "-o", unsorted_bam, sam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log")
    run_command(["samtools", "sort", "-@", config["abundance"]["threads"], "-o", all_bam, unsorted_bam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log")
    run_command(["samtools", "view", "-b", "-F", "2308", "-o", bam, all_bam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log")
    run_command(["samtools", "index", bam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log")
    idxstats = workdir / "idxstats.tsv"
    coverage_path = workdir / "coverage.tsv"
    run_command(["samtools", "idxstats", bam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log", stdout_path=idxstats)
    run_command(["samtools", "coverage", bam], env_name="qc", cwd=workdir, log_path=workdir / "samtools.log", stdout_path=coverage_path)
    rows = build_abundance_rows(records, _read_idxstats(idxstats), _read_coverage(coverage_path), total)
    write_tsv(table, fields, rows)
    groups = contig_supergroups or {}
    supergroup_rows = aggregate_supergroups(rows, groups, total)
    write_tsv(supergroup_table, ["supergroup", "n_contigs", "mapped_reads", "RPM", "mean_depth_weighted", "coverage_breadth_pct", "total_clean_reads"], supergroup_rows)
    assignment_rows = [{
        "contig_id": record.record_id,
        "supergroups": ";".join(sorted(groups.get(record.record_id, set()))),
        "status": "unique" if len(groups.get(record.record_id, set())) == 1 else "ambiguous" if len(groups.get(record.record_id, set())) > 1 else "unassigned",
    } for record in records]
    write_tsv(workdir / "abundance_contig_assignments.tsv", ["contig_id", "supergroups", "status"], assignment_rows)
    counts = {
        "contigs": len(records),
        "mapped_reads": sum(row["mapped_reads"] for row in rows),
        "total_clean_reads": total,
        "ambiguous_contigs": sum(row["status"] == "ambiguous" for row in assignment_rows),
    }
    atomic_write_json(summary, {"status": "done", "counts": counts, "formula": "RPM=primary_mapped_reads/total_clean_reads*1e6", "pe_denominator": "R1+R2 reads"})
    for path in [sam, unsorted_bam, all_bam]:
        path.unlink(missing_ok=True)
    return AbundanceResult(table, supergroup_table, bam, Path(str(bam) + ".bai"), summary, counts, "done")


def _classifier_rows(context: PipelineContext) -> list[dict]:
    value = context.data.get("classifier_assignments")
    if isinstance(value, (list, tuple)):
        return list(value)
    path = Path(value) if value else None
    if not path or not path.is_file():
        candidates = sorted(context.run_dir.glob("step*_classifier/classifier_assignments.tsv"))
        path = candidates[0] if candidates else None
    return read_tsv(path) if path and path.is_file() else []


def _effective_contig_supergroups(context: PipelineContext) -> dict[str, set[str]]:
    """member → representative → 有效超群；参考匹配结果优先于分类。"""
    representative_rows = read_tsv(context.data["representative_map"])
    representative_map = {
        row["query_id"]: row["representative_id"] for row in representative_rows
        if row.get("query_id") and row.get("representative_id")
    }
    assignments = {row["query_id"]: row for row in read_tsv(context.data["supergroup_assignments"]) if row.get("query_id")}
    classifier = {row["query_id"]: row for row in _classifier_rows(context) if row.get("query_id")}
    effective = {}
    for representative in set(representative_map.values()):
        assignment = assignments.get(representative, {})
        supergroup = str(assignment.get("supergroup") or "").strip()
        if not supergroup:
            supergroup = str(classifier.get(representative, {}).get("supergroup") or "").strip()
        if supergroup:
            effective[representative] = supergroup

    contig_groups: dict[str, set[str]] = {}
    path = context.data.get("orf_contig_map")
    if path and Path(path).is_file():
        with Path(path).open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                representative = representative_map.get(row.get("query_id", ""))
                supergroup = effective.get(representative)
                contig = row.get("contig_id", "")
                if supergroup and contig:
                    contig_groups.setdefault(contig, set()).add(supergroup)
    return contig_groups


def run(context: PipelineContext) -> StageResult:
    contig_groups = _effective_contig_supergroups(context)
    result = run_abundance(
        {"r1": context.data["clean_r1"], "r2": context.data.get("clean_r2")},
        context.data.get("positive_contigs_fna"),
        context.config,
        context.stage_dir,
        context.data.get("clean_read_count"),
        contig_groups,
    )
    return StageResult(
        artifacts={
            "abundance": artifact_path(context, result.table),
            "supergroup_abundance": artifact_path(context, result.supergroup_table),
            "assignments": artifact_path(context, context.stage_dir / "abundance_contig_assignments.tsv"),
            "bam": artifact_path(context, result.bam) if result.bam else None,
            "bai": artifact_path(context, result.bai) if result.bai else None,
            "summary": artifact_path(context, result.summary),
        },
        summary={**result.counts, "reason": "" if result.status == "done" else "no_positive_contigs"},
        data={"abundance_tsv": str(result.table), "supergroup_abundance_tsv": str(result.supergroup_table), "abundance_counts": result.counts},
        status=result.status,
    )
