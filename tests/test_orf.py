import csv

from cdc_pipeline.fasta import read_fasta
from cdc_pipeline.stages.orf import normalize_orf_outputs, parse_orf_map_row


def test_parse_orf_map_row():
    row = parse_orf_map_row(["SEQ000001_ORF1", "ctg1", "", "+", "1", "1", "300", "300", "99", "ATG", "TAA", "complete"])
    assert (row.orf_id, row.contig_id, row.aa_len) == ("SEQ000001_ORF1", "ctg1", 99)
    assert row.start == 1
    assert row.end == 300


def test_stable_orf_output_preserves_contig_mapping(tmp_path):
    raw_faa = tmp_path / "raw.faa"
    raw_map = tmp_path / "map.tsv"
    with raw_faa.open("w") as handle:
        handle.write(">SEQ000001_ORF1 [ctg1] description\nACD\n")
    raw_map.write_text(
        "orf_id\torig_id\tdesc\tstrand\tframe\tstart\tend\tnt_len\taa_len\tstart_codon\tstop_codon\tstatus\n"
        "SEQ000001_ORF1\tctg1\t\t+\t1\t1\t12\t12\t3\tATG\tTAA\tcomplete\n"
    )
    result = normalize_orf_outputs(raw_faa, raw_map, "sample", tmp_path)
    stable = list(read_fasta(result.proteins_faa))
    assert stable[0].record_id == "Q000001"
    with result.orf_contig_map.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["query_id"] == "Q000001"
    assert rows[0]["orf_id"] == "SEQ000001_ORF1"
    assert rows[0]["contig_id"] == "ctg1"
