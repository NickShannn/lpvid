"""基于参考序列比对结果的 180 超群分配。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.assets import normalize_diamond_sequence
from cdc_pipeline.fasta import read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path, write_tsv
from cdc_pipeline.supergroup_mapping import (
    get_tree_availability,
    normalize_supergroup_name,
)


@dataclass(frozen=True)
class Hit:
    query_id: str
    reference_id: str
    supergroup: str
    bitscore: float
    evalue: float
    pident: float
    qtitle: str = ""
    qlen: int = 0
    slen: int = 0
    qstart: int = 0
    qend: int = 0
    sstart: int = 0
    send: int = 0
    qcovhsp: float = 0.0
    scovhsp: float = 0.0
    length: int = 0


@dataclass
class SupergroupResult:
    assignments_tsv: Path
    hits_tsv: Path
    groups: dict[str, Path]
    no_tree_groups: dict[str, Path]
    unassigned_faa: Path
    counts: dict[str, int]
    status: str


def select_best_hit(hits: list[Hit]) -> Hit:
    if not hits:
        raise ValueError("空命中列表没有参考序列")
    return min(hits, key=lambda hit: (-hit.bitscore, hit.evalue, -hit.pident, hit.reference_id))


def summarize_assignments(query_ids: list[str], best_hits: dict[str, Hit]) -> list[dict]:
    rows = []
    for query_id in query_ids:
        hit = best_hits.get(query_id)
        rows.append({
            "query_id": query_id,
            "supergroup": hit.supergroup if hit else "",
            "status": "assigned" if hit else "unassigned",
            "reference_id": hit.reference_id if hit else "",
            "alignment_e_value": hit.evalue if hit else "",
            "alignment_score": hit.bitscore if hit else "",
            "reference_identity": hit.pident if hit else "",
            "qcovhsp": hit.qcovhsp if hit else "",
            "scovhsp": hit.scovhsp if hit else "",
        })
    return rows


def _reference_map(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {row["reference_id"]: row["supergroup"] for row in csv.DictReader(handle, delimiter="\t")}


def _read_hits(path: Path, reference_map: dict[str, str]) -> list[Hit]:
    hits = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 15:
                continue
            hits.append(Hit(
                query_id=p[0], qtitle=p[1], qlen=int(p[2]), reference_id=p[3], slen=int(p[4]),
                pident=float(p[5]), qstart=int(p[6]), qend=int(p[7]), sstart=int(p[8]), send=int(p[9]),
                qcovhsp=float(p[10]), scovhsp=float(p[11]), length=int(p[12]), evalue=float(p[13]),
                bitscore=float(p[14]), supergroup=reference_map[p[3]],
            ))
    return hits


def run_supergroup(representatives_faa: str | Path, config: dict, workdir: str | Path) -> SupergroupResult:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    records = list(read_fasta(representatives_faa))
    hits_path = workdir / "supergroup_hits.tsv"
    hits_path.write_text("")
    hits = []
    if records:
        from dataclasses import replace
        normalized_records = []
        for record in records:
            normalized, _ = normalize_diamond_sequence(record.sequence)
            normalized_records.append(replace(record, sequence=normalized))
        query = workdir / "reference_query.faa"
        write_fasta(query, normalized_records)


        db_path = Path(config["assets"]["rdrp_180_db"])
        ref_csv = Path(config["assets"]["rdrp_180_csv"])

        run_command(
            [
                config["tools"]["diamond"], "blastp", "--ultra-sensitive", "-q", query,
                "-d", db_path, "-o", hits_path, "-e", config["supergroup"]["evalue"],
                "--threads", config["supergroup"]["threads"],
                "--max-target-seqs", config["supergroup"]["max_target_seqs"],
                "--outfmt", "6", "qseqid", "qtitle", "qlen", "sseqid", "slen", "pident",
                "qstart", "qend", "sstart", "send", "qcovhsp", "scovhsp", "length", "evalue", "bitscore",
            ],
            env_name="lucaphylo", cwd=workdir, log_path=workdir / "reference_alignment.log",
        )


        raw_map = {}
        with ref_csv.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                sseqid = row.get("RdRp_sseqid", row.get("RdRP_sseqid", "")).strip()
                sg = row.get("super_group", row.get("Super_group", "")).strip()
                if sseqid and sg:
                    raw_map[sseqid] = sg


        raw_hits = _read_hits(hits_path, raw_map)
        from dataclasses import replace as _replace
        hits = [
            _replace(hit, supergroup=normalize_supergroup_name(hit.supergroup))
            for hit in raw_hits
        ]

    grouped_hits: dict[str, list[Hit]] = {}
    for hit in hits:
        grouped_hits.setdefault(hit.query_id, []).append(hit)
    best = {query_id: select_best_hit(group) for query_id, group in grouped_hits.items()}
    record_by_id = {record.header_id: record for record in records}


    tree_availability = get_tree_availability()


    assignment_rows = []
    for query_id in record_by_id:
        hit = best.get(query_id)
        if hit and hit.supergroup:
            tree_status = tree_availability.get(hit.supergroup, "no_tree")
            status = "assigned" if tree_status == "available" else f"classified_{tree_status}"
        else:
            status = "unassigned"
        assignment_rows.append({
            "query_id": query_id,
            "supergroup": hit.supergroup if hit else "",
            "status": status,
            "classification_method": "reference match",
            "reference_id": hit.reference_id if hit else "",
            "reference_identity": hit.pident if hit else "",
            "alignment_e_value": hit.evalue if hit else "",
            "alignment_score": hit.bitscore if hit else "",
            "qcovhsp": hit.qcovhsp if hit else "",
            "scovhsp": hit.scovhsp if hit else "",
        })

    assignments_path = workdir / "supergroup_assignments.tsv"
    write_tsv(assignments_path, [
        "query_id", "supergroup", "status", "classification_method",
        "reference_id", "reference_identity", "alignment_e_value", "alignment_score", "qcovhsp", "scovhsp",
    ], assignment_rows)
    hit_fields = ["query_id", "reference_id", "supergroup", "bitscore", "evalue", "pident", "qtitle", "qlen", "slen", "qstart", "qend", "sstart", "send", "qcovhsp", "scovhsp", "length"]
    write_tsv(workdir / "supergroup_hits_normalized.tsv", hit_fields, [hit.__dict__ for hit in hits])


    groups: dict[str, Path] = {}
    no_tree_groups: dict[str, Path] = {}
    group_dir = workdir / "groups"
    no_tree_dir = workdir / "no_tree_groups"
    for sg in sorted({hit.supergroup for hit in best.values() if hit.supergroup}):
        query_ids = {query_id for query_id, hit in best.items() if hit.supergroup == sg}
        seqs = [record for record in records if record.header_id in query_ids]
        safe_name = sg.replace("/", "_")
        tree_status = tree_availability.get(sg, "no_tree")
        if tree_status == "available":
            group_path = group_dir / f"{safe_name}.faa"
            write_fasta(group_path, seqs)
            groups[sg] = group_path
        else:
            no_tree_path = no_tree_dir / f"{safe_name}.faa"
            write_fasta(no_tree_path, seqs)
            no_tree_groups[sg] = no_tree_path

    unassigned = workdir / "unassigned.faa"
    unassigned_ids = {row["query_id"] for row in assignment_rows if row["status"] == "unassigned"}
    write_fasta(unassigned, [record for record in records if record.header_id in unassigned_ids])

    counts = {sg: sum(hit.supergroup == sg for hit in best.values()) for sg in set(groups) | set(no_tree_groups)}
    counts["unassigned"] = len(unassigned_ids)
    counts["classified_no_tree"] = sum(1 for r in assignment_rows if r["status"].startswith("classified_") and r["status"] != "assigned")
    counts["representatives"] = len(records)
    return SupergroupResult(assignments_path, hits_path, groups, no_tree_groups, unassigned, counts, "done" if records else "skipped")


def run(context: PipelineContext) -> StageResult:
    result = run_supergroup(context.data["representatives_faa"], context.config, context.stage_dir)
    return StageResult(
        artifacts={
            "assignments": artifact_path(context, result.assignments_tsv),
            "hits": artifact_path(context, result.hits_tsv),
            "normalized_hits": artifact_path(context, context.stage_dir / "supergroup_hits_normalized.tsv"),
            "unassigned": artifact_path(context, result.unassigned_faa),
            "groups": {sg: artifact_path(context, path) for sg, path in result.groups.items()},
            "no_tree_groups": {sg: artifact_path(context, path) for sg, path in result.no_tree_groups.items()},
        },
        summary={**result.counts, "reason": "" if result.counts["representatives"] else "no_representatives"},
        data={
            "supergroup_assignments": str(result.assignments_tsv),
            "supergroup_groups": {sg: str(path) for sg, path in result.groups.items()},
            "supergroup_no_tree_groups": {sg: str(path) for sg, path in result.no_tree_groups.items()},
            "supergroup_counts": result.counts,
            "unassigned_faa": str(result.unassigned_faa),
        },
        status=result.status,
    )
