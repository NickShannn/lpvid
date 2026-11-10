"""单一路径 RdRP 识别、ORF/contig 判定和阳性序列输出。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.fasta import read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.prediction import ModelPrediction
from cdc_pipeline.process import pick_free_gpu
from cdc_pipeline.stages.common import artifact_path, read_tsv, write_tsv
from cdc_pipeline.stages.rdrp_runner import parse_rdrp, run_rdrp, screen_rdrp_input


@dataclass
class DetectionRow:
    query_id: str
    rdrp_prob: float | None
    rdrp_positive: bool
    rdrp_status: str
    call: str
    analysis_set: str = "main"


def classify_predictions(
    predictions: dict[str, ModelPrediction],
    threshold: float = 0.5,
    excluded: set[str] | None = None,
) -> list[DetectionRow]:
    excluded = excluded or set()
    rows = []
    for query_id in dict.fromkeys([*predictions, *sorted(excluded)]):
        prediction = predictions.get(query_id)
        positive = bool(prediction and prediction.prob >= threshold)
        if query_id in excluded:
            status = "excluded_input_validation"
        elif prediction:
            status = "ok"
        else:
            status = "missing"
        rows.append(DetectionRow(
            query_id=query_id,
            rdrp_prob=prediction.prob if prediction else None,
            rdrp_positive=positive,
            rdrp_status=status,
            call="positive" if positive else "negative",
        ))
    return rows


def summarize_contigs(rows: list[DetectionRow], query_to_contig: dict[str, str]) -> list[dict]:
    grouped: dict[str, list[DetectionRow]] = {}
    for row in rows:
        if row.analysis_set == "main" and row.query_id in query_to_contig:
            grouped.setdefault(query_to_contig[row.query_id], []).append(row)
    output = []
    for contig_id, members in grouped.items():
        probs = [row.rdrp_prob for row in members if row.rdrp_prob is not None]
        positive = any(row.rdrp_positive for row in members)
        output.append({
            "contig_id": contig_id,
            "n_orf_model_input": len(members),
            "rdrp_max_prob": max(probs) if probs else "",
            "call": "positive" if positive else "negative",
        })
    return output


def run(context: PipelineContext) -> StageResult:
    input_faa = Path(context.data["model_input_faa"])
    records = list(read_fasta(input_faa))
    prefilter = {row["query_id"]: row for row in read_tsv(context.data["prefilter_calls"])}
    clean, excluded = screen_rdrp_input(records)
    clean_faa = context.stage_dir / "rdrp_input.faa"
    write_fasta(clean_faa, clean)
    excluded_faa = context.stage_dir / "rdrp_excluded.faa"
    write_fasta(excluded_faa, excluded)
    if clean:
        gpu = pick_free_gpu(context.gpu)
        prediction_path = run_rdrp(clean_faa, context.config, context.stage_dir / "rdrp", gpu)
        predictions = parse_rdrp(prediction_path)
        missing = [record.header_id for record in clean if record.header_id not in predictions]
        if missing:
            raise ValueError(f"RdRP 结果缺失: {missing[:5]}")
    else:
        prediction_path = context.stage_dir / "rdrp" / "rdrp_predictions.csv"
        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        prediction_path.write_text("seq_id,seq,prob,label_index,label\n", encoding="utf-8")
        predictions = {}
    rows = classify_predictions(
        predictions, context.config["model_threshold"], {record.header_id for record in excluded},
    )
    for row in rows:
        row.analysis_set = prefilter.get(row.query_id, {}).get("analysis_set", "main")
    record_by_id = {record.header_id: record for record in records}
    positive = [record_by_id[row.query_id] for row in rows if row.analysis_set == "main" and row.rdrp_positive]
    positive_faa = context.stage_dir / "rdrp_positive.faa"
    write_fasta(positive_faa, positive)
    fields = ["query_id", "rdrp_prob", "rdrp_positive", "rdrp_status", "call", "analysis_set"]
    predictions_tsv = context.stage_dir / "orf_predictions.tsv"
    write_tsv(predictions_tsv, fields, [row.__dict__ for row in rows])

    query_to_contig = {}
    positive_contigs_faa = None
    if context.data.get("orf_contig_map"):
        query_to_contig = {row["query_id"]: row["contig_id"] for row in read_tsv(context.data["orf_contig_map"])}
    contig_rows = summarize_contigs(rows, query_to_contig)
    contig_tsv = context.stage_dir / "contig_level_calls.tsv"
    write_tsv(contig_tsv, ["contig_id", "n_orf_model_input", "rdrp_max_prob", "call"], contig_rows)
    if context.data.get("contigs_fna"):
        positive_ids = {row["contig_id"] for row in contig_rows if row["call"] == "positive"}
        positive_contigs_faa = context.stage_dir / "rdrp_positive_contigs.fna"
        write_fasta(positive_contigs_faa, [record for record in read_fasta(context.data["contigs_fna"]) if record.record_id in positive_ids])
    counts = {
        "model_input": len(records),
        "rdrp_excluded": len(excluded),
        "main_rdrp_positive_orfs": len(positive),
        "qa_rdrp_positive_orfs": sum(row.analysis_set == "qa_discard" and row.rdrp_positive for row in rows),
        "positive_contigs": sum(row["call"] == "positive" for row in contig_rows),
        "reason": "" if records else "no_prefilter_candidates",
    }
    return StageResult(
        artifacts={
            "orf_predictions": artifact_path(context, predictions_tsv),
            "contig_calls": artifact_path(context, contig_tsv),
            "positive_faa": artifact_path(context, positive_faa),
            "positive_contigs_fna": artifact_path(context, positive_contigs_faa) if positive_contigs_faa else None,
            "rdrp_input": artifact_path(context, clean_faa),
            "rdrp_predictions": artifact_path(context, prediction_path),
            "rdrp_excluded": artifact_path(context, excluded_faa),
        },
        summary=counts,
        data={
            "positive_faa": str(positive_faa),
            "orf_predictions": str(predictions_tsv),
            "contig_calls": str(contig_tsv),
            "positive_contigs_fna": str(positive_contigs_faa) if positive_contigs_faa else None,
            "detect_counts": counts,
            "rdrp_predictions": str(prediction_path),
            "rdrp_excluded": str(excluded_faa),
        },
        status="done" if records else "skipped",
    )
