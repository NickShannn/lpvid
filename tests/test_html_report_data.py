import json
import math
import re
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import pytest

import cdc_pipeline.html_report as html_report
from cdc_pipeline.html_report import join_sequence_rows, length_statistics, supergroup_shares
from cdc_pipeline.stages.common import read_tsv, write_tsv


def test_supergroup_shares_include_unassigned_and_sum_100():
    rows = [
        {"query_id": "Q1", "supergroup": "SG1"},
        {"query_id": "Q2", "supergroup": "SG1"},
        {"query_id": "Q3", "supergroup": ""},
    ]
    shares = supergroup_shares(rows)
    assert shares[0]["count"] == 2
    assert math.isclose(sum(row["percent"] for row in shares), 100.0)
    assert shares[-1]["supergroup"] == "Unassigned"
    assert shares[-1]["percent"] == 33.33


def test_join_sequence_rows_connects_all_result_tables(tmp_path):
    rep = write_tsv(tmp_path / "rep.tsv", ["query_id", "representative_id", "aa_len", "motif_class_count", "motif_classes", "is_representative"], [
        {"query_id": "Q1", "representative_id": "Q1", "aa_len": 200, "motif_class_count": 2, "motif_classes": "A;C", "is_representative": 1},
        {"query_id": "Q2", "representative_id": "Q1", "aa_len": 190, "motif_class_count": 1, "motif_classes": "C", "is_representative": 0},
    ])
    detection = write_tsv(tmp_path / "detect.tsv", ["query_id", "call", "rdrp_prob"], [{"query_id": "Q1", "call": "positive", "rdrp_prob": 0.8}])
    assignments = write_tsv(tmp_path / "sg.tsv", ["query_id", "supergroup", "reference_id", "alignment_e_value", "alignment_score", "reference_identity"], [{"query_id": "Q1", "supergroup": "SG1", "reference_id": "R1", "alignment_e_value": 1e-5, "alignment_score": 100, "reference_identity": 80}])
    placements = write_tsv(tmp_path / "place.tsv", ["query_id", "nn1", "nn1_dist"], [{"query_id": "Q1", "nn1": "tip1", "nn1_dist": 0.2}])
    orf_map = write_tsv(tmp_path / "map.tsv", ["query_id", "contig_id"], [{"query_id": "Q1", "contig_id": "c1"}])
    abundance = write_tsv(tmp_path / "ab.tsv", ["contig_id", "RPM", "mean_depth", "coverage_breadth_pct"], [{"contig_id": "c1", "RPM": 10, "mean_depth": 2, "coverage_breadth_pct": 100}])
    context = SimpleNamespace(data={
        "representative_map": str(rep), "orf_predictions": str(detection), "supergroup_assignments": str(assignments),
        "placements_tsv": str(placements), "orf_contig_map": str(orf_map), "abundance_tsv": str(abundance),
    })
    rows = join_sequence_rows(context)
    assert len(rows) == 1
    row = rows[0]
    assert row["query_id"] == "Q1"
    assert row["aa_len"] == 200
    assert row["supergroup"] == "SG1"
    assert row["nn1"] == "tip1"
    assert row["contig_id"] == "c1"
    assert row["original_orf_id"] == ""
    assert row["rpm"] == 10.0
    for legacy in ("host_human", "host_mammal", "host_vertebrate", "host_suppression_reason"):
        assert legacy not in row


def _reconciled_context(tmp_path):
    representative_rows = [
        {"query_id": "R1", "representative_id": "R1", "aa_len": 300, "is_representative": 1},
        {"query_id": "M1", "representative_id": "R1", "aa_len": 290, "is_representative": 0},
        {"query_id": "M2", "representative_id": "R1", "aa_len": 280, "is_representative": 0},
        {"query_id": "R2", "representative_id": "R2", "aa_len": 200, "is_representative": 1},
        {"query_id": "R3", "representative_id": "R3", "aa_len": 220, "is_representative": 1},
        {"query_id": "R4", "representative_id": "R4", "aa_len": 400, "is_representative": 1},
        {"query_id": "M4", "representative_id": "R4", "aa_len": 390, "is_representative": 0},
        {"query_id": "R5_Picobirnaviridae", "representative_id": "R5_Picobirnaviridae", "aa_len": 240, "is_representative": 1},
        {"query_id": "R6", "representative_id": "R6", "aa_len": "", "is_representative": 1},
    ]
    rep = write_tsv(
        tmp_path / "rep.tsv",
        ["query_id", "representative_id", "aa_len", "is_representative"],
        representative_rows,
    )
    detection = write_tsv(
        tmp_path / "detect.tsv",
        ["query_id", "call", "rdrp_prob"],
        [
            {"query_id": "R1", "call": "positive", "rdrp_prob": "NaN"},
            {"query_id": "UNRELATED", "call": "positive", "rdrp_prob": "not-a-number"},
        ],
    )
    assignments = write_tsv(
        tmp_path / "sg.tsv",
        ["query_id", "supergroup", "status", "classification_method", "reference_id", "alignment_e_value", "alignment_score", "reference_identity"],
        [
            {"query_id": "R1", "supergroup": "SG1", "status": "assigned", "classification_method": "reference match", "alignment_e_value": "1e-5"},
            {"query_id": "R2", "supergroup": "SG1", "status": "assigned", "classification_method": "reference match", "alignment_e_value": "1e9999"},
            {"query_id": "R3", "supergroup": "", "status": "unassigned", "classification_method": "reference match"},
            {"query_id": "R4", "supergroup": "SG_NO_TREE", "status": "classified_no_tree", "classification_method": "reference match"},
            {"query_id": "R5_Picobirnaviridae", "supergroup": "SG_PBV", "status": "assigned", "classification_method": "reference match"},
            {"query_id": "R6", "supergroup": "", "status": "unassigned", "classification_method": "reference match"},
        ],
    )
    classifier = [
        {"query_id": "R3", "supergroup": "SG_EMB", "tree_status": "no_tree"},
        {"query_id": "R1", "supergroup": "SG_MUST_NOT_OVERRIDE_REFERENCE", "tree_status": "available"},
    ]
    orf_map = write_tsv(
        tmp_path / "map.tsv",
        ["query_id", "orf_id", "contig_id"],
        [
            {"query_id": "R1", "orf_id": "o1", "contig_id": "c1"},
            {"query_id": "M1", "orf_id": "o2", "contig_id": "c1"},
            {"query_id": "M2", "orf_id": "o3", "contig_id": "c2"},
            {"query_id": "R2", "orf_id": "o4", "contig_id": "c3"},
            {"query_id": "R3", "orf_id": "o5", "contig_id": "c3"},
            {"query_id": "R4", "orf_id": "o6", "contig_id": "c4"},
            {"query_id": "M4", "orf_id": "o7", "contig_id": "c5"},
            {"query_id": "R5_Picobirnaviridae", "orf_id": "o8", "contig_id": "c6"},
            {"query_id": "UNRELATED", "orf_id": "o9", "contig_id": "c999"},
        ],
    )
    abundance = write_tsv(
        tmp_path / "ab.tsv",
        ["contig_id", "length", "mapped_reads", "RPM", "mean_depth", "coverage_breadth_pct", "total_clean_reads"],
        [
            {"contig_id": "c1", "length": 100, "mapped_reads": 10, "RPM": 100, "mean_depth": 2, "coverage_breadth_pct": 50, "total_clean_reads": 100000},
            {"contig_id": "c2", "length": 100, "mapped_reads": 20, "RPM": 200, "mean_depth": 4, "coverage_breadth_pct": 100, "total_clean_reads": 100000},
            {"contig_id": "c3", "length": 200, "mapped_reads": 30, "RPM": "1e9999", "mean_depth": 6, "coverage_breadth_pct": 75, "total_clean_reads": 100000},
            {"contig_id": "c4", "length": 100, "mapped_reads": 40, "RPM": 400, "mean_depth": 8, "coverage_breadth_pct": 80, "total_clean_reads": 100000},
            {"contig_id": "c6", "length": 100, "mapped_reads": 60, "RPM": 600, "mean_depth": 12, "coverage_breadth_pct": 100, "total_clean_reads": 100000},
            {"contig_id": "c999", "length": 100, "mapped_reads": 999, "RPM": 9990, "mean_depth": 99, "coverage_breadth_pct": 100, "total_clean_reads": 100000},
        ],
    )
    context = SimpleNamespace(
        data={
            "representative_map": str(rep),
            "orf_predictions": str(detection),
            "supergroup_assignments": str(assignments),
            "classifier_assignments": classifier,
            "orf_contig_map": str(orf_map),
            "abundance_tsv": str(abundance),
        },
        run_dir=tmp_path,
        route="reads",
        sample="fixture",
        config={"project": {"display_name": "test"}},
    )
    summary = {
        "risk_rows": [
            {"query_id": "R1", "family": "NotPicobirnaviridae"},
            SimpleNamespace(query_id="R5_Picobirnaviridae", family="Picobirnaviridae"),
        ]
    }
    return context, summary


def test_enrichment_reconciles_all_positive_members_and_effective_supergroups(tmp_path, monkeypatch):
    context, summary = _reconciled_context(tmp_path)
    real_read_tsv = html_report.read_tsv

    def reject_bulk_reads(path):
        assert Path(path).name not in {"detect.tsv", "map.tsv"}, "large ORF tables must be streamed"
        return real_read_tsv(path)

    monkeypatch.setattr(html_report, "read_tsv", reject_bulk_reads)
    rows = join_sequence_rows(context)
    enrich = getattr(html_report, "enrich_report_summary", None)
    assert callable(enrich)
    assert enrich(context, rows, summary) is summary

    stats = {row["supergroup"]: row for row in summary["supergroup_stats"]}
    assert stats["SG1"]["representatives"] == 2
    assert stats["SG1"]["positive_orfs"] == 4
    assert stats["SG1"]["positive_contigs"] == 3
    assert stats["SG1"]["positive_orf_ids"] == ["M1", "M2", "R1", "R2"]
    assert stats["SG1"]["contig_ids"] == ["c1", "c2", "c3"]
    assert stats["SG1"]["length_statistics"]["count"] == 2
    assert stats["SG_EMB"]["classification_methods"] == {"classifier": 1}
    assert stats["SG_EMB"]["tree_statuses"] == {"classified_no_tree": 1}
    assert stats["SG_NO_TREE"]["positive_orfs"] == 2
    assert stats["SG_NO_TREE"]["contig_ids"] == ["c4", "c5"]
    assert stats["SG_NO_TREE"]["tree_statuses"] == {"classified_no_tree": 1}
    assert stats["Unassigned"]["representatives"] == 1
    assert stats["Unassigned"]["length_statistics"]["count"] == 0

    assert summary["multi_supergroup_contig_count"] == 1
    assert summary["multi_supergroup_contigs"] == [
        {"contig_id": "c3", "supergroups": ["SG1", "SG_EMB"]}
    ]
    assert stats["SG1"]["abundance"]["n_contigs"] == 2
    assert stats["SG1"]["abundance"]["mapped_reads"] == 30
    assert stats["SG1"]["abundance"]["RPM"] == 300.0
    assert stats["SG1"]["abundance"]["mean_depth_weighted"] == 3.0
    assert stats["SG1"]["abundance"]["coverage_breadth_pct"] == 75.0
    assert stats["SG_EMB"]["abundance"]["n_contigs"] == 0
    assert stats["SG_EMB"]["abundance"]["RPM"] is None
    assert stats["SG_EMB"]["abundance"]["mapped_reads"] is None
    assert stats["SG_NO_TREE"]["abundance"]["mapped_reads"] == 40
    for stat in stats.values():
        assert "host" not in stat

    expected_available_contigs = {"c1", "c2", "c3", "c4", "c6"}
    assert {row["contig_id"] for row in summary["abundance_rows"]} == expected_available_contigs
    assert "c5" not in {row["contig_id"] for row in summary["abundance_rows"]}
    assert "c999" not in {row["contig_id"] for row in summary["abundance_rows"]}
    assert summary["_abundance_rows"] == summary["abundance_rows"]
    for legacy in ("host_rows", "host_coverage", "host_suppressed_contigs", "host_pbv_unpredicted_contigs", "_host_rows"):
        assert legacy not in summary

    by_query = {row["query_id"]: row for row in rows}
    assert by_query["R1"]["rdrp_prob"] == ""
    assert by_query["R2"]["evalue"] == ""
    assert stats["SG1"]["classification_methods"] == {"reference match": 2}
    c3_abundance = next(row for row in summary["abundance_rows"] if row["contig_id"] == "c3")
    assert c3_abundance["rpm"] == ""
    json.dumps(
        {
            "rows": rows,
            "supergroup_stats": summary["supergroup_stats"],
            "abundance_rows": summary["abundance_rows"],
            "multi_supergroup_contigs": summary["multi_supergroup_contigs"],
        },
        allow_nan=False,
    )


def test_bilingual_build_reuses_cached_filtered_tables(tmp_path, monkeypatch):
    import cdc_pipeline.html_template_v2 as template

    context, summary = _reconciled_context(tmp_path)
    summary["risk_rows"] = []
    reads = Counter()
    real_open = Path.open

    def counting_open(self, mode="r", *args, **kwargs):
        if "r" in mode and self.name in {"rep.tsv", "detect.tsv", "map.tsv", "ab.tsv"}:
            reads[self.name] += 1
        return real_open(self, mode, *args, **kwargs)

    rendered = []

    def fake_render(context, rows, summary, lang="zh"):
        rendered.append((lang, len(rows), len(summary["supergroup_stats"])))
        return lang

    monkeypatch.setattr(Path, "open", counting_open)
    monkeypatch.setattr(template, "render_html", fake_render)
    html_report.build_html_report(context, tmp_path / "zh.html", summary, lang="zh")
    html_report.build_html_report(context, tmp_path / "en.html", summary, lang="en")
    assert rendered == [("zh", 6, 5), ("en", 6, 5)]
    assert reads == {"rep.tsv": 1, "detect.tsv": 1, "map.tsv": 1, "ab.tsv": 1}


def _missing_abundance_context(tmp_path):
    rep = write_tsv(
        tmp_path / "rep.tsv",
        ["query_id", "representative_id", "aa_len", "is_representative"],
        [
            {"query_id": "Q1", "representative_id": "Q1", "aa_len": 200, "is_representative": 1},
            {"query_id": "Q2", "representative_id": "Q2", "aa_len": 210, "is_representative": 1},
        ],
    )
    assignments = write_tsv(
        tmp_path / "sg.tsv", ["query_id", "supergroup", "status"],
        [
            {"query_id": "Q1", "supergroup": "SG_MISSING", "status": "assigned"},
            {"query_id": "Q2", "supergroup": "SG_ZERO", "status": "assigned"},
        ],
    )
    orf_map = write_tsv(
        tmp_path / "map.tsv",
        ["query_id", "orf_id", "contig_id"],
        [
            {"query_id": "Q1", "orf_id": "o1", "contig_id": "c1"},
            {"query_id": "Q2", "orf_id": "o2", "contig_id": "c2"},
        ],
    )
    abundance = write_tsv(
        tmp_path / "ab.tsv",
        ["contig_id", "length", "mapped_reads", "RPM", "mean_depth", "coverage_breadth_pct", "total_clean_reads"],
        [
            {"contig_id": "c1", "length": 100, "mapped_reads": "", "RPM": "NaN",
             "mean_depth": "NaN", "coverage_breadth_pct": "", "total_clean_reads": ""},
            {"contig_id": "c2", "length": 100, "mapped_reads": 0, "RPM": 0,
             "mean_depth": 0, "coverage_breadth_pct": 0, "total_clean_reads": 0},
        ],
    )
    context = SimpleNamespace(
        data={
            "representative_map": str(rep),
            "supergroup_assignments": str(assignments),
            "orf_contig_map": str(orf_map),
            "abundance_tsv": str(abundance),
        },
        run_dir=tmp_path,
        route="reads",
        sample="fixture",
        config={"project": {"display_name": "test"}},
    )
    summary = {"run_id": "fixture", "route": "reads", "sample": "fixture", "counts": {"representatives": 2}, "risk_rows": [], "links": {}}
    return context, summary


def _card(text, name):
    body = re.sub(r"<script.*?</script>", "", text, flags=re.S)
    return re.search(rf'<article class="supergroup-card"><h3>{name}</h3>.*?</article>', body, re.S).group(0)


def test_weighted_mean_returns_none_without_valid_pairs():
    assert html_report._weighted_mean([{"mean_depth": "", "length": 100}], "mean_depth") is None
    assert html_report._weighted_mean([{"mean_depth": "NaN", "length": ""}], "mean_depth") is None
    assert html_report._weighted_mean([{"mean_depth": 0.0, "length": 100}], "mean_depth") == 0.0
    assert html_report._weighted_mean([{"mean_depth": 2.0, "length": 100}], "mean_depth") == 2.0


def test_abundance_summary_computes_only_metrics_with_valid_values():
    rows = {"c1": {"rpm": "", "mapped_reads": 12, "mean_depth": 2.0, "coverage_breadth_pct": "", "length": 100}}
    summary = html_report._abundance_summary({"c1"}, set(), rows)
    assert summary["mapped_reads"] == 12
    assert summary["RPM"] is None
    assert summary["mean_depth_weighted"] == 2.0
    assert summary["coverage_breadth_pct"] is None


def test_enrichment_missing_abundance_metrics_are_unavailable_but_real_zero_is_kept(tmp_path):
    from cdc_pipeline.html_template_v2 import render_html

    context, summary = _missing_abundance_context(tmp_path)
    rows = join_sequence_rows(context)
    html_report.enrich_report_summary(context, rows, summary)
    stats = {row["supergroup"]: row for row in summary["supergroup_stats"]}

    missing = stats["SG_MISSING"]["abundance"]
    assert missing["n_contigs"] == 1
    assert missing["mapped_reads"] is None
    assert missing["RPM"] is None
    assert missing["mean_depth_weighted"] is None
    assert missing["coverage_breadth_pct"] is None

    zero = stats["SG_ZERO"]["abundance"]
    assert zero["n_contigs"] == 1
    assert zero["mapped_reads"] == 0
    assert zero["RPM"] == 0.0
    assert zero["mean_depth_weighted"] == 0.0
    assert zero["coverage_breadth_pct"] == 0.0

    text = render_html(context, rows, summary, lang="en")
    missing_card = _card(text, "SG_MISSING")
    assert "RPM unavailable" in missing_card
    assert "mapped reads unavailable" in missing_card
    assert "weighted depth unavailable" in missing_card
    assert "coverage unavailable" in missing_card
    for pattern in (r"RPM 0", r"mapped reads 0", r"weighted depth 0", r"coverage 0"):
        assert not re.search(pattern, missing_card, re.I)

    zero_card = _card(text, "SG_ZERO")
    assert "RPM 0" in zero_card
    assert "mapped reads 0" in zero_card
    assert "weighted depth 0" in zero_card
    assert "coverage 0%" in zero_card


def test_length_statistics_and_empty_inputs():
    stats = length_statistics([{"aa_len": 100}, {"aa_len": 200}, {"aa_len": 300}])
    assert stats["min"] == 100
    assert stats["median"] == 200
    assert stats["max"] == 300
    assert sum(bin_["count"] for bin_ in stats["bins"]) == 3
    assert length_statistics([])["bins"] == []
    assert supergroup_shares([]) == []


def test_enrichment_without_upstream_probability_table_never_fabricates_rows(tmp_path):
    context, summary = _reconciled_context(tmp_path)
    context.data.pop("orf_contig_map", None)
    rows = join_sequence_rows(context)
    html_report.enrich_report_summary(context, rows, summary)

    assert summary["abundance_rows"] == []
    for key in summary:
        assert "host" not in key.lower(), key
