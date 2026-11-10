"""长度分层、RefSeq DIAMOND 粗筛、三路分诊和 DISCARD QA。"""

from __future__ import annotations

import random
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

from cdc_pipeline.assets import normalize_diamond_sequence
from cdc_pipeline.fasta import FastaRecord, assign_stable_ids, read_fasta, write_fasta, write_id_map
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path, write_tsv


@dataclass(frozen=True)
class PrefilterHit:
    qseqid: str
    staxids: str
    sseqid: str
    pident: float
    qcovhsp: float
    qlen: int
    bitscore: float

    def taxids(self) -> set[int]:
        result = set()
        for value in self.staxids.split(";"):
            try:
                result.add(int(value))
            except ValueError:
                continue
        return result


@dataclass
class LengthPartition:
    too_short: list[FastaRecord]
    holdout: list[FastaRecord]
    main: list[FastaRecord]

    @property
    def all_ge_input(self) -> list[FastaRecord]:
        return self.holdout + self.main


@dataclass
class PrefilterResult:
    model_input_faa: Path
    calls_tsv: Path
    hits_tsv: Path
    qa_faa: Path
    counts: dict[str, int]
    status: str
    reason: str


def partition_lengths(records: Iterable[FastaRecord], input_min: int = 100, analysis_min: int = 150) -> LengthPartition:
    too_short, holdout, main = [], [], []
    for record in records:
        if record.length < input_min:
            too_short.append(record)
        elif record.length < analysis_min:
            holdout.append(record)
        else:
            main.append(record)
    return LengthPartition(too_short, holdout, main)


def _has_taxid(hit: PrefilterHit, allowed: set[int]) -> bool:
    return bool(hit.taxids() & allowed)


def triage_hits(
    hits: Iterable[PrefilterHit],
    viral_taxids: set[int],
    cellular_taxids: set[int],
    pident: float = 70,
    qcov: float = 70,
) -> dict[str, str]:
    grouped: dict[str, list[PrefilterHit]] = {}
    for hit in hits:
        grouped.setdefault(hit.qseqid, []).append(hit)
    calls = {}
    for query_id, query_hits in grouped.items():
        if any(_has_taxid(hit, viral_taxids) for hit in query_hits):
            calls[query_id] = "KEEP_viral"
        elif any(_has_taxid(hit, cellular_taxids) and hit.pident >= pident and hit.qcovhsp >= qcov for hit in query_hits):
            calls[query_id] = "DISCARD"
        else:
            calls[query_id] = "KEEP_unknown"
    return calls


def resolve_triage(
    query_ids: Iterable[str],
    hits: Iterable[PrefilterHit],
    viral_taxids: set[int],
    cellular_taxids: set[int],
    pident: float = 70,
    qcov: float = 70,
) -> dict[str, str]:
    calls = {query_id: "KEEP_unknown" for query_id in query_ids}
    calls.update(triage_hits(hits, viral_taxids, cellular_taxids, pident, qcov))
    return calls


def load_taxids(path: str | Path) -> set[int]:
    result = set()
    with Path(path).open() as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                try:
                    result.add(int(value))
                except ValueError:
                    continue
    return result


def read_hits(path: str | Path) -> list[PrefilterHit]:
    hits = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            hits.append(PrefilterHit(
                qseqid=parts[0], staxids=parts[1], sseqid=parts[2],
                pident=float(parts[3]), qcovhsp=float(parts[4]), qlen=int(parts[5]), bitscore=float(parts[6]),
            ))
    return hits


def choose_qa(discards: list[FastaRecord], fraction: float, seed: int) -> list[FastaRecord]:
    if not discards or fraction <= 0:
        return []
    count = max(1, int(len(discards) * fraction))
    return random.Random(seed).sample(discards, min(count, len(discards)))


def run_prefilter(input_faa: str | Path, sample: str, config: dict, workdir: Path) -> PrefilterResult:
    workdir.mkdir(parents=True, exist_ok=True)
    records, id_mapping = assign_stable_ids(read_fasta(input_faa), sample)
    stable_faa = workdir / "input.stable.faa"
    write_fasta(stable_faa, records)
    write_id_map(workdir / "input_id_map.tsv", id_mapping)
    partition = partition_lengths(records, config["min_input_aa"], config["min_analysis_aa"])
    write_fasta(workdir / "below_min_input.faa", partition.too_short)
    write_fasta(workdir / "short_holdout.faa", partition.holdout)
    write_fasta(workdir / "ge_min_input.faa", partition.all_ge_input)
    write_fasta(workdir / "main_analysis.faa", partition.main)

    hits_path = workdir / "refseq_hits.tsv"
    hits_path.write_text("")
    calls: dict[str, str] = {record.header_id: "KEEP_unknown" for record in partition.main}
    hits: list[PrefilterHit] = []
    enabled = bool(config["prefilter"]["enabled"])
    cellular: set[int] = set()
    if enabled and partition.main:
        diamond_records = []
        for record in partition.main:
            normalized, _ = normalize_diamond_sequence(record.sequence)
            diamond_records.append(replace(record, sequence=normalized))
        diamond_query = workdir / "diamond_query.faa"
        write_fasta(diamond_query, diamond_records)
        tmpdir = workdir / "diamond_tmp"
        tmpdir.mkdir(exist_ok=True)
        run_command(
            [
                config["tools"]["diamond"], "blastp", "-d", config["assets"]["refseq_db"],
                "-q", diamond_query, "--outfmt", "6", "qseqid", "staxids", "sseqid", "pident",
                "qcovhsp", "qlen", "bitscore", "--max-target-seqs", config["prefilter"]["max_target_seqs"],
                "--threads", config["prefilter"]["threads"], "--evalue", config["prefilter"]["evalue"],
                "--tmpdir", tmpdir, "-o", hits_path,
            ],
            env_name="lucaphylo", cwd=workdir, log_path=workdir / "diamond.log",
        )
        hits = read_hits(hits_path)
        viral = load_taxids(config["assets"]["viral_taxids"])
        cellular = load_taxids(config["assets"]["cellular_taxids"])
        calls = resolve_triage(
            [record.header_id for record in partition.main], hits, viral, cellular,
            config["prefilter"]["pident"], config["prefilter"]["qcov"],
        )

    keep_viral = [record for record in partition.main if calls[record.header_id] == "KEEP_viral"]
    keep_unknown = [record for record in partition.main if calls[record.header_id] == "KEEP_unknown"]
    discards = [record for record in partition.main if calls[record.header_id] == "DISCARD"]
    write_fasta(workdir / "KEEP_viral.faa", keep_viral)
    write_fasta(workdir / "KEEP_unknown.faa", keep_unknown)
    qa_records = choose_qa(discards, config["prefilter"]["qa_fraction"], config["prefilter"]["qa_seed"])
    qa_path = workdir / "qa_discard_sample.faa"
    write_fasta(qa_path, qa_records)
    qa_ids = {record.header_id for record in qa_records}
    model_records = keep_viral + keep_unknown + qa_records
    model_input = workdir / "model_input.faa"
    write_fasta(model_input, model_records)

    best_cell = {}
    for hit in hits:
        if hit.qseqid in calls and calls[hit.qseqid] == "DISCARD" and hit.qseqid not in best_cell \
                and _has_taxid(hit, cellular) \
                and hit.pident >= config["prefilter"]["pident"] and hit.qcovhsp >= config["prefilter"]["qcov"]:
            best_cell[hit.qseqid] = hit
    write_tsv(
        workdir / "DISCARD.tsv",
        ["qseqid", "best_cellular_sseqid", "staxids", "pident", "qcovhsp"],
        [{
            "qseqid": record.header_id,
            "best_cellular_sseqid": getattr(best_cell.get(record.header_id), "sseqid", ""),
            "staxids": getattr(best_cell.get(record.header_id), "staxids", ""),
            "pident": getattr(best_cell.get(record.header_id), "pident", ""),
            "qcovhsp": getattr(best_cell.get(record.header_id), "qcovhsp", ""),
        } for record in discards],
    )
    call_rows = [{
        "sample": sample,
        "query_id": record.header_id,
        "aa_len": record.length,
        "prefilter_call": calls[record.header_id],
        "analysis_set": "qa_discard" if record.header_id in qa_ids else "main",
    } for record in partition.main]
    calls_path = workdir / "prefilter_calls.tsv"
    write_tsv(calls_path, ["sample", "query_id", "aa_len", "prefilter_call", "analysis_set"], call_rows)
    counts = {
        "input": len(records),
        "below_min_input": len(partition.too_short),
        "short_holdout": len(partition.holdout),
        "main": len(partition.main),
        "KEEP_viral": len(keep_viral),
        "KEEP_unknown": len(keep_unknown),
        "DISCARD": len(discards),
        "qa_discard": len(qa_records),
        "model_input": len(model_records),
        "diamond_hits": len(hits),
        "prefilter_enabled": int(enabled),
    }
    status = "done" if model_records else "skipped"
    reason = "" if model_records else "no_prefilter_candidates"
    return PrefilterResult(model_input, calls_path, hits_path, qa_path, counts, status, reason)


def run(context: PipelineContext) -> StageResult:
    input_faa = context.data.get("proteins_faa") or context.inputs["input"]
    result = run_prefilter(input_faa, context.sample, context.config, context.stage_dir)
    artifacts = {
        name: artifact_path(context, context.stage_dir / name)
        for name in [
            "input.stable.faa", "input_id_map.tsv", "below_min_input.faa", "short_holdout.faa",
            "ge_min_input.faa", "main_analysis.faa", "refseq_hits.tsv", "KEEP_viral.faa",
            "KEEP_unknown.faa", "DISCARD.tsv", "qa_discard_sample.faa", "model_input.faa",
            "prefilter_calls.tsv",
        ]
    }
    return StageResult(
        artifacts=artifacts,
        summary={**result.counts, "reason": result.reason},
        data={
            "model_input_faa": str(result.model_input_faa),
            "prefilter_calls": str(result.calls_tsv),
            "prefilter_counts": result.counts,
            "prefilter_status": result.reason,
        },
        status=result.status,
    )
