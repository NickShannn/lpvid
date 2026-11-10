"""RdRP 阳性序列精确去重、CD-HIT 聚类和确定性代表选择。"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from cdc_pipeline.fasta import FastaRecord, read_fasta, write_fasta
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.process import run_command
from cdc_pipeline.stages.common import artifact_path, write_tsv

MOTIF_PATTERNS = {
    "F": re.compile(r"G.K.{3,5}G"),
    "A": re.compile(r"D.{3,5}[DE]"),
    "B": re.compile(r"[ST]G.{2,4}[ST].{2,4}[NQH]"),
    "C": re.compile(r"GD[DN]"),
}


@dataclass(frozen=True)
class Member:
    query_id: str
    aa_len: int
    motif_classes: set[str]

    @property
    def motif_class_count(self) -> int:
        return len(self.motif_classes)


@dataclass
class RepresentativeResult:
    representatives_faa: Path
    mapping_tsv: Path
    motifs_tsv: Path
    counts: dict[str, int]
    status: str


def scan_motifs(sequence: str) -> dict[str, list[int]]:
    return {name: [match.start() for match in pattern.finditer(sequence)] for name, pattern in MOTIF_PATTERNS.items()}


def choose_representative(cluster_members: list[Member]) -> Member:
    if not cluster_members:
        raise ValueError("空簇没有代表序列")
    return min(cluster_members, key=lambda member: (-member.motif_class_count, -member.aa_len, member.query_id))


def parse_cdhit_clusters(path: str | Path) -> list[list[str]]:
    clusters = []
    current = None
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">Cluster"):
                current = []
                clusters.append(current)
            elif line and current is not None and ">" in line:
                current.append(line.split(">", 1)[1].split("...", 1)[0])
    return clusters


def _member(record: FastaRecord) -> Member:
    motifs = scan_motifs(record.sequence)
    return Member(record.header_id, record.length, {name for name, positions in motifs.items() if positions})


def run_representatives(detection_faa: str | Path, config: dict, workdir: str | Path) -> RepresentativeResult:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    records = list(read_fasta(detection_faa))
    record_by_id = {record.header_id: record for record in records}
    exact_groups: dict[str, list[FastaRecord]] = {}
    for record in records:
        exact_groups.setdefault(record.sequence, []).append(record)
    canonical_records = [min(group, key=lambda record: record.header_id) for group in exact_groups.values()]
    canonical_records.sort(key=lambda record: record.header_id)
    canonical_by_id = {record.header_id: record.header_id for record in canonical_records}
    exact_members = {
        canonical.header_id: sorted(group, key=lambda record: record.header_id)
        for canonical in canonical_records
        for group in [exact_groups[canonical.sequence]]
    }
    unique_faa = workdir / "exact_unique.faa"
    write_fasta(unique_faa, canonical_records)

    cluster_ids = []
    if len(canonical_records) == 1:
        cluster_ids = [[canonical_records[0].header_id]]
    elif canonical_records:
        clustered_faa = workdir / "cdhit.faa"
        run_command(
            [
                config["tools"]["cd_hit"], "-i", unique_faa, "-o", clustered_faa,
                "-c", config["cluster"]["identity"], "-aS", config["cluster"]["coverage"],
                "-n", "5", "-T", config["cluster"]["threads"], "-M", "0",
            ],
            env_name="lucaphylo", cwd=workdir, log_path=workdir / "cdhit.log",
        )
        cluster_ids = parse_cdhit_clusters(workdir / "cdhit.faa.clstr")

    output_rows = []
    representatives = []
    for index, cluster in enumerate(cluster_ids, 1):
        cluster_id = f"C{index:06d}"
        expanded = [record for canonical_id in cluster for record in exact_members[canonical_id]]
        representative = choose_representative([_member(record) for record in expanded])
        representatives.append(record_by_id[representative.query_id])
        canonical_for_member = {
            member.header_id: canonical_id
            for canonical_id in cluster
            for member in exact_members[canonical_id]
        }
        for record in expanded:
            member = _member(record)
            output_rows.append({
                "cluster_id": cluster_id,
                "query_id": record.header_id,
                "representative_id": representative.query_id,
                "aa_len": record.length,
                "motif_class_count": member.motif_class_count,
                "motif_classes": ";".join(sorted(member.motif_classes)),
                "is_representative": int(record.header_id == representative.query_id),
                "exact_duplicate_of": canonical_for_member[record.header_id] if record.header_id != canonical_for_member[record.header_id] else "",
            })
    mapped_ids = {row["query_id"] for row in output_rows}
    if mapped_ids != set(record_by_id):
        raise ValueError(f"CD-HIT 簇成员映射不完整: {len(mapped_ids)}/{len(record_by_id)}")

    representatives_faa = workdir / "representatives.faa"
    mapping_tsv = workdir / "representative_map.tsv"
    motifs_tsv = workdir / "motifs.tsv"
    write_fasta(representatives_faa, representatives)
    write_tsv(mapping_tsv, [
        "cluster_id", "query_id", "representative_id", "aa_len", "motif_class_count",
        "motif_classes", "is_representative", "exact_duplicate_of",
    ], output_rows)
    motif_rows = []
    for record in records:
        motifs = scan_motifs(record.sequence)
        motif_rows.append({
            "query_id": record.header_id,
            "aa_len": record.length,
            **{name: ";".join(map(str, motifs[name])) for name in ["F", "A", "B", "C"]},
            "motif_class_count": sum(bool(motifs[name]) for name in ["F", "A", "B", "C"]),
            "motif_classes": ";".join(name for name in ["F", "A", "B", "C"] if motifs[name]),
        })
    write_tsv(motifs_tsv, ["query_id", "aa_len", "F", "A", "B", "C", "motif_class_count", "motif_classes"], motif_rows)
    counts = {
        "positive_input": len(records),
        "exact_unique": len(canonical_records),
        "clusters": len(cluster_ids),
        "representatives": len(representatives),
        "mapped_members": len(output_rows),
    }
    return RepresentativeResult(representatives_faa, mapping_tsv, motifs_tsv, counts, "done" if records else "skipped")


def run(context: PipelineContext) -> StageResult:
    result = run_representatives(context.data["positive_faa"], context.config, context.stage_dir)
    return StageResult(
        artifacts={
            "representatives_faa": artifact_path(context, result.representatives_faa),
            "representative_map": artifact_path(context, result.mapping_tsv),
            "motifs": artifact_path(context, result.motifs_tsv),
        },
        summary={**result.counts, "reason": "" if result.counts["positive_input"] else "no_rdrp_positive"},
        data={
            "representatives_faa": str(result.representatives_faa),
            "representative_map": str(result.mapping_tsv),
            "motifs": str(result.motifs_tsv),
            "representative_counts": result.counts,
        },
        status=result.status,
    )
