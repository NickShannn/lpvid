from cdc_pipeline.prediction import ModelPrediction
from cdc_pipeline.stages.detect import classify_predictions, summarize_contigs


def pred(value):
    return ModelPrediction("", value, "test")


def rows_by_id(rows):
    return {row.query_id: row for row in rows}


def test_single_rdrp_threshold_and_probability_provenance():
    rows = classify_predictions({"q1": pred(0.1), "q2": pred(0.9)}, threshold=0.5)
    assert rows_by_id(rows)["q1"].call == "negative"
    assert rows_by_id(rows)["q2"].call == "positive"
    assert rows_by_id(rows)["q1"].rdrp_prob == 0.1
    assert rows_by_id(rows)["q2"].rdrp_prob == 0.9
    assert rows_by_id(rows)["q1"].rdrp_positive is False
    assert rows_by_id(rows)["q2"].rdrp_positive is True


def test_excluded_input_is_not_a_negative_probability():
    rows = classify_predictions({}, excluded={"q1"})
    row = rows_by_id(rows)["q1"]
    assert row.call == "negative"
    assert row.rdrp_prob is None
    assert row.rdrp_status == "excluded_input_validation"


def test_contig_any_orf_positive_rule():
    rows = classify_predictions({"q1": pred(0.9), "q2": pred(0.1)}, threshold=0.5)
    mapping = {"q1": "ctg1", "q2": "ctg1"}
    calls = summarize_contigs(rows, mapping)
    assert len(calls) == 1
    assert calls[0]["contig_id"] == "ctg1"
    assert calls[0]["call"] == "positive"
    assert calls[0]["rdrp_max_prob"] == 0.9
