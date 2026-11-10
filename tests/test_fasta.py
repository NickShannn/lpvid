import gzip

from cdc_pipeline.fasta import assign_stable_ids, read_fasta, stable_id, write_fasta


def test_fasta_preserves_original_id_and_assigns_stable_id(tmp_path):
    path = tmp_path / "x.faa"
    path.write_text(">orig id with spaces\nACD\n")
    records = list(read_fasta(path))
    assert records[0].record_id == "orig"
    assert records[0].description == "orig id with spaces"
    assert stable_id("sample", 1) == "Q000001"
    assigned, mapping = assign_stable_ids(records, "sample")
    assert assigned[0].query_id == "Q000001"
    assert mapping[0]["original_id"] == "orig"
    assert mapping[0]["sample"] == "sample"


def test_fasta_roundtrip_gzip_and_empty_file(tmp_path):
    source = tmp_path / "x.faa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">a\nACD\n>b description\nEFG\n")
    records = list(read_fasta(source))
    assigned, _ = assign_stable_ids(records, "s")
    output = tmp_path / "out.faa"
    write_fasta(output, assigned)
    reread = list(read_fasta(output))
    assert [r.record_id for r in reread] == ["Q000001", "Q000002"]
    assert [r.sequence for r in reread] == ["ACD", "EFG"]
    empty = tmp_path / "empty.faa"
    write_fasta(empty, [])
    assert list(read_fasta(empty)) == []
