from pathlib import Path

import pytest

from cdc_pipeline.config import load_config
from cdc_pipeline.stages.supergroup import Hit, select_best_hit, summarize_assignments
from cdc_pipeline.supergroup_mapping import load_reference_supergroup_map


def hit(bits, evalue, sg, reference_id):
    return Hit("q1", reference_id, sg, bits, evalue, 90.0)


def test_best_hit_uses_bitscore_then_evalue():
    hits = [hit(100, 1e-5, "A", "R1"), hit(100, 1e-20, "B", "2"), hit(90, 0.0, "C", "R3")]
    assert select_best_hit(hits).supergroup == "B"


def test_unassigned_queries_are_not_dropped():
    rows = summarize_assignments(["q1", "q2"], {})
    assert rows[0]["status"] == "unassigned"
    assert rows[1]["status"] == "unassigned"
    assert [row["query_id"] for row in rows] == ["q1", "q2"]


def test_assignment_rows_use_generic_evidence_field_names():
    row = summarize_assignments(["q1"], {"q1": hit(100, 1e-5, "A", "R1")})[0]
    assert row["reference_identity"] == 90.0
    assert row["alignment_e_value"] == 1e-5
    assert row["alignment_score"] == 100
    assert "pident" not in row
    assert "evalue" not in row
    assert "bitscore" not in row


@pytest.mark.skipif(
    not Path(load_config()["assets"]["rdrp_180_db"]).with_suffix(".fas").is_file(),
    reason="外部 RdRP 参考库资源包未配置",
)
def test_rdrp_180_fasta_ids_match_reference_map():
    config = load_config()
    reference_map = load_reference_supergroup_map()
    fasta_path = Path(config["assets"]["rdrp_180_db"]).with_suffix(".fas")
    fasta_ids = []
    with fasta_path.open() as handle:
        fasta_ids = [line[1:].split()[0] for line in handle if line.startswith(">")]

    missing = sorted(set(fasta_ids) - set(reference_map))
    assert len(fasta_ids) == len(set(fasta_ids)), "Duplicate sequence IDs in RdRP reference FASTA"
    assert not missing, f"{len(missing)} FASTA IDs missing from reference CSV: {missing}"
