import csv

from cdc_pipeline.fasta import FastaRecord, read_fasta, write_fasta
from cdc_pipeline.stages.representatives import (
    Member,
    choose_representative,
    parse_cdhit_clusters,
    run_representatives,
    scan_motifs,
)


def member(query_id, length, motifs):
    return Member(query_id, length, motifs)


def test_motif_scan_uses_existing_fabc_definitions():
    motifs = scan_motifs("ACDGDD")
    assert motifs["C"] == [3]
    assert motifs["F"] == []
    assert sum(bool(value) for value in motifs.values()) == 1


def test_representative_prioritizes_motif_classes_then_length_then_id():
    members = [member("q1", 500, {"A", "B"}), member("q2", 300, {"A", "B", "C"})]
    assert choose_representative(members).query_id == "q2"
    assert choose_representative([member("q2", 300, set()), member("q1", 300, set())]).query_id == "q1"


def test_parse_cdhit_clusters(tmp_path):
    path = tmp_path / "x.clstr"
    path.write_text(">Cluster 0\n0\t100aa, >Q1... *\n1\t90aa, >Q2... at 90%\n>Cluster 1\n0\t80aa, >Q3... *\n")
    assert parse_cdhit_clusters(path) == [["Q1", "Q2"], ["Q3"]]


def test_exact_duplicates_expand_back_to_cluster_mapping(tmp_path):
    source = tmp_path / "positive.faa"
    write_fasta(source, [FastaRecord("q1", "ACDGDD"), FastaRecord("q2", "ACDGDD")])
    config = {"cluster": {"identity": 0.9, "coverage": 0.8, "threads": 2}, "tools": {"cd_hit": "cd-hit"}}
    result = run_representatives(source, config, tmp_path / "work")
    representatives = list(read_fasta(result.representatives_faa))
    assert len(representatives) == 1
    with result.mapping_tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["representative_id"] for row in rows} == {"q1"}
