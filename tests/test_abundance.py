from types import SimpleNamespace

from cdc_pipeline.fasta import FastaRecord
from cdc_pipeline.stages import abundance as abundance_module
from cdc_pipeline.stages.abundance import (
    AbundanceResult,
    aggregate_supergroups,
    build_abundance_rows,
    calculate_rpm,
    clean_read_denominator,
)
from cdc_pipeline.stages.common import write_tsv


def test_rpm_and_paired_denominator():
    assert calculate_rpm(100, 1_000_000) == 100.0
    assert calculate_rpm(100, 0) == 0.0
    assert clean_read_denominator(1000, 1000) == 2000


def test_zero_mapping_contigs_are_preserved():
    records = [FastaRecord("c1", "A" * 100), FastaRecord("c2", "C" * 200)]
    rows = build_abundance_rows(records, {"c1": 10}, {"c1": {"meandepth": 2.0, "coverage": 50.0}}, 1000)
    assert rows[0]["mapped_reads"] == 10
    assert rows[0]["RPM"] == 10000.0
    assert rows[1]["mapped_reads"] == 0
    assert rows[1]["coverage_breadth_pct"] == 0.0


def test_ambiguous_contigs_are_not_double_counted():
    rows = [
        {"contig_id": "c1", "length": 100, "mapped_reads": 10, "RPM": 10, "mean_depth": 2, "coverage_breadth_pct": 50},
        {"contig_id": "c2", "length": 200, "mapped_reads": 20, "RPM": 20, "mean_depth": 4, "coverage_breadth_pct": 100},
    ]
    aggregated = aggregate_supergroups(rows, {"c1": {"SG1", "SG2"}, "c2": {"SG1"}}, 1000)
    assert len(aggregated) == 1
    assert aggregated[0]["supergroup"] == "SG1"
    assert aggregated[0]["mapped_reads"] == 20


def test_run_maps_all_members_through_every_effective_supergroup(tmp_path, monkeypatch):
    representative_map = write_tsv(
        tmp_path / "representative_map.tsv",
        ["query_id", "representative_id", "is_representative"],
        [
            {"query_id": "R1", "representative_id": "R1", "is_representative": 1},
            {"query_id": "M1", "representative_id": "R1", "is_representative": 0},
            {"query_id": "R2", "representative_id": "R2", "is_representative": 1},
            {"query_id": "R3", "representative_id": "R3", "is_representative": 1},
        ],
    )
    assignments = write_tsv(
        tmp_path / "assignments.tsv",
        ["query_id", "supergroup", "status"],
        [
            {"query_id": "R1", "supergroup": "SG_NO_TREE", "status": "classified_no_tree"},
            {"query_id": "R2", "supergroup": "", "status": "unassigned"},
            {"query_id": "R3", "supergroup": "SG_PENDING", "status": "classified_pending_merge"},
        ],
    )
    orf_map = write_tsv(
        tmp_path / "orf_contig_map.tsv",
        ["query_id", "contig_id"],
        [
            {"query_id": "R1", "contig_id": "c1"},
            {"query_id": "M1", "contig_id": "c2"},
            {"query_id": "R2", "contig_id": "c1"},
            {"query_id": "R3", "contig_id": "c3"},
            {"query_id": "UNRELATED", "contig_id": "c4"},
        ],
    )
    context = SimpleNamespace(
        data={
            "representative_map": str(representative_map),
            "supergroup_assignments": str(assignments),
            "classifier_assignments": [
                {"query_id": "R1", "supergroup": "SG_MUST_NOT_OVERRIDE", "tree_status": "available"},
                {"query_id": "R2", "supergroup": "SG_EMB", "tree_status": "no_tree"},
            ],
            "orf_contig_map": str(orf_map),
            "clean_r1": "unused-r1.fastq",
            "positive_contigs_fna": "",
            "clean_read_count": 100,
        },
        config={},
        run_dir=tmp_path,
        stage_dir=tmp_path,
    )
    captured = {}

    def fake_run_abundance(clean_reads, positive_contigs, config, workdir, total_clean_reads=None, contig_supergroups=None):
        captured.update(contig_supergroups)
        return AbundanceResult(
            table=tmp_path / "abundance.tsv",
            supergroup_table=tmp_path / "supergroup_abundance.tsv",
            bam=None,
            bai=None,
            summary=tmp_path / "abundance_summary.json",
            counts={"contigs": 0, "mapped_reads": 0, "total_clean_reads": 100, "ambiguous_contigs": 1},
            status="done",
        )

    monkeypatch.setattr(abundance_module, "run_abundance", fake_run_abundance)
    result = abundance_module.run(context)
    assert result.status == "done"
    assert captured == {
        "c1": {"SG_NO_TREE", "SG_EMB"},
        "c2": {"SG_NO_TREE"},
        "c3": {"SG_PENDING"},
    }
