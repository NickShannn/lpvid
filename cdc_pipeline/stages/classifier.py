"""对参考序列比对未命中的序列执行通用分类。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.fasta import read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import pick_free_gpu, run_command
from cdc_pipeline.stages.common import artifact_path, write_tsv
from cdc_pipeline.supergroup_mapping import get_tree_availability


@dataclass
class ClassifierAssignment:
    query_id: str
    supergroup: str
    confidence: float
    top1_label: str
    top1_prob: float
    top2_label: str
    top2_prob: float
    top3_label: str
    top3_prob: float
    tree_status: str


def parse_predictions(csv_path: Path) -> list[dict]:
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _load_label_names(classifier_dir: Path) -> dict[str, str]:
    labels_path = classifier_dir / "labels_supergroups.tsv"
    mapping = {}
    if labels_path.is_file():
        with labels_path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                label = row.get("label", "").strip()
                name = row.get("supergroup_name", "").strip()
                if label and name and name != "unnamed_novel":
                    mapping[label] = name
    return mapping


def run_classifier(
    unassigned_faa: str | Path,
    config: dict,
    workdir: str | Path,
    gpu: str | int,
) -> list[ClassifierAssignment]:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    records = list(read_fasta(unassigned_faa))
    if not records:
        return []

    classifier_dir = Path(config["assets"]["classifier_pipeline"])
    output_dir = workdir / "classifier_output"
    output_dir.mkdir(exist_ok=True)

    selected_gpu = pick_free_gpu(gpu)
    run_command(
        [
            config["environments"]["lucaphylo"],
            classifier_dir / "classify.py",
            "--input", Path(unassigned_faa).resolve(),
            "--output_dir", output_dir,
            "--cuda_device", "0",
        ],
        env_name="lucaphylo",
        cwd=classifier_dir,
        log_path=workdir / "classifier.log",
        env={
            "CUDA_VISIBLE_DEVICES": str(selected_gpu),
            "PYTHONNOUSERSITE": "1",
        },
    )

    predictions_path = output_dir / "superclade_predictions.csv"
    if not predictions_path.is_file():
        raise FileNotFoundError(f"分类器未产出: {predictions_path}")

    raw_predictions = parse_predictions(predictions_path)
    tree_availability = get_tree_availability()

    label_names = _load_label_names(classifier_dir)
    results = []
    for pred in raw_predictions:
        raw_sg = pred.get("predicted_supergroup", "").strip()

        sg = label_names.get(raw_sg, raw_sg)
        tree_status = tree_availability.get(sg, "no_tree")
        results.append(ClassifierAssignment(
            query_id=pred["seq_id"],
            supergroup=sg,
            confidence=float(pred.get("confidence", 0)),
            top1_label=label_names.get(pred.get("top1_label", ""), pred.get("top1_label", sg)),
            top1_prob=float(pred.get("top1_prob", pred.get("confidence", 0))),
            top2_label=label_names.get(pred.get("top2_label", ""), pred.get("top2_label", "")),
            top2_prob=float(pred.get("top2_prob", 0)),
            top3_label=label_names.get(pred.get("top3_label", ""), pred.get("top3_label", "")),
            top3_prob=float(pred.get("top3_prob", 0)),
            tree_status=tree_status,
        ))
    return results


def run(context: PipelineContext) -> StageResult:
    unassigned_faa = Path(context.data.get("unassigned_faa", ""))
    if not unassigned_faa.is_file():
        return StageResult(
            artifacts={},
            summary={"status": "skipped", "reason": "no_unassigned"},
            data={"classifier_assignments": [], "classifier_status": "skipped_no_unassigned"},
            status="skipped",
        )

    assignments = run_classifier(
        unassigned_faa, context.config, context.stage_dir, context.gpu,
    )

    tsv_path = context.stage_dir / "classifier_assignments.tsv"
    fields = [
        "query_id", "supergroup", "confidence",
        "top1_label", "top1_prob", "top2_label", "top2_prob", "top3_label", "top3_prob",
        "tree_status",
    ]
    write_tsv(tsv_path, fields, [row.__dict__ for row in assignments])

    tree_availability = get_tree_availability()
    groups_with_tree: dict[str, Path] = {}
    groups_no_tree: dict[str, Path] = {}
    group_dir = context.stage_dir / "groups"
    no_tree_dir = context.stage_dir / "no_tree_groups"

    records = {r.header_id: r for r in read_fasta(unassigned_faa)}
    for sg in sorted({row.supergroup for row in assignments if row.supergroup}):
        query_ids = {row.query_id for row in assignments if row.supergroup == sg}
        seqs = [records[qid] for qid in query_ids if qid in records]
        safe_name = sg.replace("/", "_")
        status = tree_availability.get(sg, "no_tree")
        if status == "available":
            path = group_dir / f"{safe_name}.faa"
            write_fasta(path, seqs)
            groups_with_tree[sg] = path
        else:
            path = no_tree_dir / f"{safe_name}.faa"
            write_fasta(path, seqs)
            groups_no_tree[sg] = path

    remaining_faa = context.stage_dir / "remaining_unclassified.faa"
    classified_ids = {row.query_id for row in assignments if row.supergroup}
    remaining = [r for r in records.values() if r.header_id not in classified_ids]
    write_fasta(remaining_faa, remaining)

    counts = {
        "classified": len([row for row in assignments if row.supergroup]),
        "classified_with_tree": len([row for row in assignments if row.tree_status == "available"]),
        "classified_no_tree": len([row for row in assignments if row.supergroup and row.tree_status != "available"]),
        "remaining_unclassified": len(remaining),
    }

    return StageResult(
        artifacts={
            "assignments": artifact_path(context, tsv_path),
            "remaining": artifact_path(context, remaining_faa),
            "groups": {sg: artifact_path(context, path) for sg, path in groups_with_tree.items()},
            "no_tree_groups": {sg: artifact_path(context, path) for sg, path in groups_no_tree.items()},
        },
        summary={**counts, "reason": "" if assignments else "no_predictions"},
        data={
            "classifier_assignments": [row.__dict__ for row in assignments],
            "classifier_groups_with_tree": {sg: str(path) for sg, path in groups_with_tree.items()},
            "classifier_groups_no_tree": {sg: str(path) for sg, path in groups_no_tree.items()},
            "classifier_counts": counts,
            "classifier_status": "done" if assignments else "skipped_no_predictions",
        },
        status="done" if assignments else "skipped",
    )
