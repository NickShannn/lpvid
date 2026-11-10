from cdc_pipeline.fasta import FastaRecord
from cdc_pipeline.stages.prefilter import (
    PrefilterHit,
    choose_qa,
    partition_lengths,
    resolve_triage,
    triage_hits,
)


def hit(qseqid, taxid, pident, qcov):
    return PrefilterHit(qseqid, str(taxid), "subject", pident, qcov, 500, 100.0)


def record(record_id, length):
    return FastaRecord(record_id=record_id, sequence="A" * length, description=record_id)


def test_viral_hit_wins_over_cellular_hit():
    hits = [hit("q1", 2, 100, 100), hit("q1", 10239, 10, 10)]
    assert triage_hits(hits, {10239}, {2})["q1"] == "KEEP_viral"


def test_cellular_discard_requires_both_thresholds():
    hits = [hit("q1", 2, 70, 70), hit("q2", 2, 70, 69.9), hit("q3", 2, 69.9, 100)]
    calls = triage_hits(hits, set(), {2})
    assert calls["q1"] == "DISCARD"
    assert calls["q2"] == "KEEP_unknown"
    assert calls["q3"] == "KEEP_unknown"


def test_no_hit_is_explicitly_kept_unknown():
    calls = resolve_triage(["q1", "q2"], [hit("q1", 2, 100, 100)], set(), {2})
    assert calls == {"q1": "DISCARD", "q2": "KEEP_unknown"}


def test_short_sequences_are_holdout_not_deleted():
    part = partition_lengths([record("a", 99), record("b", 100), record("c", 149), record("d", 150)])
    assert [r.record_id for r in part.too_short] == ["a"]
    assert [r.record_id for r in part.holdout] == ["b", "c"]
    assert [r.record_id for r in part.main] == ["d"]
    assert [r.record_id for r in part.all_ge_input] == ["b", "c", "d"]


def test_unified_150aa_threshold_has_no_holdout():
    part = partition_lengths([record("a", 149), record("b", 150)], input_min=150, analysis_min=150)
    assert [r.record_id for r in part.too_short] == ["a"]
    assert [r.record_id for r in part.main] == ["b"]
    assert part.holdout == []


def test_qa_sample_is_nonempty_and_deterministic():
    discards = [record(str(i), 200) for i in range(10)]
    first = choose_qa(discards, 0.01, 42)
    second = choose_qa(discards, 0.01, 42)
    assert len(first) == 1
    assert [r.record_id for r in first] == [r.record_id for r in second]
    assert choose_qa([], 0.01, 42) == []
