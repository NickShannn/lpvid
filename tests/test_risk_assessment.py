import json

import pytest

from cdc_pipeline.risk_assessment import (
    PathogenEntry,
    RiskRow,
    TaxonomyEntry,
    build_risk_rows,
    summarize_risk,
)

MISSING = object()


def _tax(
    tip_id: str,
    *,
    family: str = "",
    genus: str = "",
    human: bool = False,
    vertebrate: bool = False,
    plant: bool = False,
) -> TaxonomyEntry:
    return TaxonomyEntry(
        tip_id=tip_id,
        family=family,
        genus=genus,
        species="",
        cluster_taxon="",
        human=human,
        vertebrate=vertebrate,
        plant=plant,
    )


def test_human_nn1_is_counted_per_row_not_by_global_minimum_distance():
    taxonomy = {
        "human_a": _tax("human_a", family="Flaviviridae", genus="Flavivirus", human=True, vertebrate=True),
        "human_b": _tax("human_b", family="Flaviviridae", genus="Flavivirus", human=True, vertebrate=True),
        "other_a": _tax("other_a", family="Flaviviridae", genus="Flavivirus"),
        "other_b": _tax("other_b", family="Flaviviridae", genus="Flavivirus"),
    }
    placements = [
        {"query_id": "Q1", "supergroup": "SG1", "nn1": "human_a", "nn1_dist": "0.10",
         "nn2": "other_a", "nn2_dist": "2.00"},
        {"query_id": "Q2", "supergroup": "SG1", "nn1": "human_b", "nn1_dist": "0.95",
         "nn2": "other_b", "nn2_dist": "2.50"},
    ]

    rows = build_risk_rows(placements, taxonomy, [])
    summary = summarize_risk(rows)

    # Both rows have a human nn1 even though only one is the global minimum distance.
    assert summary.get("human_nn1", 0) == 2
    assert summary.get("human_nn2_5", 0) == 0
    assert summary.get("vertebrate_only", 0) == 0
    assert summary.get("no_placed_label", 0) == 0
    assert summary["total"] == 2

    for row, expected_dist in zip(rows, (0.10, 0.95)):
        assert [neighbor.rank for neighbor in row.neighbors] == [1, 2]
        assert row.neighbors[0].reference_id == row.nn1
        assert row.neighbors[0].distance == pytest.approx(expected_dist)
        assert row.neighbors[0].human is True
        assert row.neighbors[0].vertebrate is True
        assert row.neighbors[1].reference_id.startswith("other_")
        assert row.neighbors[1].distance == pytest.approx(2.0 if row.query_id == "Q1" else 2.5)


def test_human_nn2_5_is_counted_separately_and_keeps_true_rank():
    taxonomy = {
        "nonhuman_nn1": _tax("nonhuman_nn1", family="Picornaviridae", genus="Enterovirus"),
        "human_nn2": _tax("human_nn2", family="Coronaviridae", genus="Betacoronavirus", human=True, vertebrate=True),
    }
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "nonhuman_nn1", "nn1_dist": "0.20",
        "nn2": "human_nn2", "nn2_dist": "0.70",
        "nn3": "tip_without_taxonomy", "nn3_dist": "1.40",
    }]

    rows = build_risk_rows(placements, taxonomy, [])
    summary = summarize_risk(rows)

    assert summary.get("human_nn1", 0) == 0
    assert summary.get("human_nn2_5", 0) == 1
    assert summary.get("vertebrate_only", 0) == 0
    assert summary.get("no_placed_label", 0) == 0

    row = rows[0]
    assert [(neighbor.rank, neighbor.reference_id) for neighbor in row.neighbors] == [
        (1, "nonhuman_nn1"), (2, "human_nn2"), (3, "tip_without_taxonomy"),
    ]
    assert row.neighbors[1].rank == 2
    assert row.neighbors[1].distance == pytest.approx(0.70)
    assert row.neighbors[1].human is True
    assert row.neighbors[1].vertebrate is True
    assert row.neighbors[2].rank == 3
    assert row.neighbors[2].family == ""
    assert row.neighbors[2].genus == ""
    # Unknown taxonomy => unknown host evidence, never a negative claim.
    assert row.neighbors[2].human is None
    assert row.neighbors[2].vertebrate is None
    assert row.neighbors[2].plant is None


def test_vertebrate_only_and_no_placed_label_are_mutually_exclusive():
    taxonomy = {
        "plain_nn1": _tax("plain_nn1", family="Picornaviridae", genus="Enterovirus"),
        "vertebrate_nn2": _tax("vertebrate_nn2", family="Picornaviridae", genus="Enterovirus", vertebrate=True),
        "plain_other": _tax("plain_other", family="Picornaviridae", genus="Enterovirus"),
        "plant_nn3": _tax("plant_nn3", family="Virgaviridae", genus="Tobamovirus", plant=True),
    }
    placements = [
        {"query_id": "Q1", "supergroup": "SG1", "nn1": "plain_nn1", "nn1_dist": "0.3",
         "nn2": "vertebrate_nn2", "nn2_dist": "0.8", "nn3": "plant_nn3", "nn3_dist": "1.1"},
        {"query_id": "Q2", "supergroup": "SG1", "nn1": "plain_other", "nn2": "plain_other"},
    ]

    rows = build_risk_rows(placements, taxonomy, [])
    summary = summarize_risk(rows)

    assert summary.get("human_nn1", 0) == 0
    assert summary.get("human_nn2_5", 0) == 0
    assert summary.get("vertebrate_only", 0) == 1
    assert summary.get("no_placed_label", 0) == 1
    assert (
        summary["human_nn1"] + summary["human_nn2_5"]
        + summary["vertebrate_only"] + summary["no_placed_label"]
    ) == summary["total"] == 2


def test_catalog_and_statutory_matching_stays_anchored_to_nn1_family():
    taxonomy = {
        "nn1_tip": _tax("nn1_tip", family="Unlistedviridae", genus="Unlistedvirus"),
        "nn2_tip": _tax("nn2_tip", family="Coronaviridae", genus="Betacoronavirus", human=True, vertebrate=True),
    }
    catalog = [
        PathogenEntry(
            seq="1", name_cn="示例冠状病毒", name_en="Example coronavirus",
            family_cn="冠状病毒科", family_en="Coronaviridae",
            category="第二类", statutory="乙类", covered="yes",
        )
    ]
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "nn1_tip", "nn1_dist": "0.2",
        "nn2": "nn2_tip", "nn2_dist": "0.9",
    }]

    row = build_risk_rows(placements, taxonomy, catalog)[0]

    # Catalog/statutory fields come from the nn1 family only.
    assert row.family == "Unlistedviridae"
    assert row.genus == "Unlistedvirus"
    assert row.pathogen_category == ""
    assert row.statutory_category == ""
    assert row.pathogen_examples == ""
    assert row.statutory_examples == ""
    assert "未匹配到上述名录收录成员" in row.hint
    assert "不排除名录外或未识别病毒" in row.hint
    assert "需进一步核查" in row.hint
    assert "不构成新序列风险分级或病原确证" in row.hint
    assert "常规监测" not in row.hint
    # Evidence for nn2 remains available without contaminating the nn1 catalog lookup.
    assert row.neighbors[1].family == "Coronaviridae"
    assert row.neighbors[1].genus == "Betacoronavirus"
    assert row.neighbors[1].human is True


def test_blank_nn_columns_are_skipped_without_renumbering_ranks():
    taxonomy = {
        "first": _tax("first", family="FamA", genus="GenusA"),
        "third": _tax("third", family="FamB", genus="GenusB", vertebrate=True),
    }
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "first", "nn1_dist": "0.1",
        "nn2": "", "nn2_dist": "",
        "nn3": "third", "nn3_dist": "1.3",
        "nn4": "", "nn4_dist": "",
        "nn5": "", "nn5_dist": "",
    }]

    row = build_risk_rows(placements, taxonomy, [])[0]

    assert [(neighbor.rank, neighbor.reference_id) for neighbor in row.neighbors] == [
        (1, "first"), (3, "third"),
    ]
    assert row.neighbors[1].rank == 3
    assert row.neighbors[1].distance == pytest.approx(1.3)
    assert row.neighbors[1].vertebrate is True


def test_risk_row_and_neighbor_serialize_to_structured_dicts():
    taxonomy = {
        "human_tip": _tax("human_tip", family="Flaviviridae", genus="Flavivirus", human=True, vertebrate=True),
        "other_tip": _tax("other_tip", family="Flaviviridae", genus="Flavivirus"),
    }
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "human_tip", "nn1_dist": "0.15",
        "nn2": "other_tip", "nn2_dist": "1.5",
    }]

    row = build_risk_rows(placements, taxonomy, [])[0]

    neighbor_payload = row.neighbors[0].to_dict()
    assert isinstance(neighbor_payload, dict)
    assert neighbor_payload["rank"] == 1
    assert neighbor_payload["reference_id"] == "human_tip"
    assert neighbor_payload["distance"] == pytest.approx(0.15)
    assert neighbor_payload["family"] == "Flaviviridae"
    assert neighbor_payload["genus"] == "Flavivirus"
    assert neighbor_payload["human"] is True
    assert neighbor_payload["vertebrate"] is True
    assert neighbor_payload["plant"] is False
    assert not isinstance(row.neighbors[0].to_dict(), str)

    payload = row.to_dict()
    assert isinstance(payload, dict)
    assert isinstance(payload["neighbors"], list)
    assert all(isinstance(item, dict) for item in payload["neighbors"])
    assert payload["neighbors"][1]["rank"] == 2
    assert payload["neighbors"][1]["reference_id"] == "other_tip"

    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False)
    decoded = json.loads(encoded)
    assert isinstance(decoded, dict)
    assert decoded["query_id"] == "Q1"
    assert isinstance(decoded["neighbors"][0], dict)
    assert decoded["neighbors"][0]["rank"] == 1
    assert decoded["neighbors"][1]["distance"] == pytest.approx(1.5)


def test_nn1_and_nn2_both_human_counts_only_human_nn1():
    taxonomy = {
        "nn1_human": _tax("nn1_human", family="Flaviviridae", genus="Flavivirus", human=True),
        "nn2_human": _tax("nn2_human", family="Flaviviridae", genus="Flavivirus", human=True),
    }
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "nn1_human", "nn1_dist": "0.2",
        "nn2": "nn2_human", "nn2_dist": "0.8",
    }]

    rows = build_risk_rows(placements, taxonomy, [])
    summary = summarize_risk(rows)

    assert summary["human_nn1"] == 1
    assert summary["human_nn2_5"] == 0
    assert summary["vertebrate_only"] == 0
    assert summary["no_placed_label"] == 0
    assert [neighbor.human for neighbor in rows[0].neighbors] == [True, True]


def test_nn4_human_and_nn5_vertebrate_evidence_is_collected_and_summarized():
    taxonomy = {
        "plain_nn1": _tax("plain_nn1", family="Picornaviridae", genus="Enterovirus"),
        "human_nn4": _tax("human_nn4", family="Coronaviridae", genus="Betacoronavirus", human=True),
        "vertebrate_nn5": _tax("vertebrate_nn5", family="Astroviridae", genus="Avastrovirus", vertebrate=True),
    }
    placements = [
        {"query_id": "Q1", "supergroup": "SG1", "nn1": "plain_nn1", "nn1_dist": "0.2",
         "nn2": "", "nn3": "", "nn4": "human_nn4", "nn4_dist": "1.4"},
        {"query_id": "Q2", "supergroup": "SG1", "nn1": "plain_nn1", "nn1_dist": "0.3",
         "nn2": "", "nn3": "", "nn4": "", "nn5": "vertebrate_nn5", "nn5_dist": "1.8"},
    ]

    rows = build_risk_rows(placements, taxonomy, [])
    summary = summarize_risk(rows)

    assert [(n.rank, n.reference_id) for n in rows[0].neighbors] == [(1, "plain_nn1"), (4, "human_nn4")]
    assert rows[0].neighbors[1].human is True
    assert rows[0].neighbors[1].distance == pytest.approx(1.4)
    assert [(n.rank, n.reference_id) for n in rows[1].neighbors] == [(1, "plain_nn1"), (5, "vertebrate_nn5")]
    assert rows[1].neighbors[1].vertebrate is True
    assert rows[1].neighbors[1].distance == pytest.approx(1.8)

    assert summary["human_nn1"] == 0
    assert summary["human_nn2_5"] == 1
    assert summary["vertebrate_only"] == 1
    assert summary["no_placed_label"] == 0
    assert (
        summary["human_nn1"] + summary["human_nn2_5"]
        + summary["vertebrate_only"] + summary["no_placed_label"]
    ) == summary["total"] == 2


def test_nn1_catalog_hits_describe_member_categories_in_both_dimensions():
    taxonomy = {
        "corona_tip": _tax("corona_tip", family="Coronaviridae", genus="Betacoronavirus"),
    }
    catalog = [
        PathogenEntry(
            seq="1", name_cn="示例冠状病毒", name_en="Example coronavirus",
            family_cn="冠状病毒科", family_en="Coronaviridae",
            category="第二类", statutory="乙类", covered="yes",
        )
    ]
    placements = [{
        "query_id": "Q1", "supergroup": "SG1",
        "nn1": "corona_tip", "nn1_dist": "0.1",
    }]

    row = build_risk_rows(placements, taxonomy, catalog)[0]

    assert row.pathogen_category == "第二类"
    assert row.statutory_category == "乙类"
    assert "Coronaviridae（nn1 所属科）的已知成员中" in row.hint
    assert "有名录收录成员涉及" in row.hint
    assert "《病原微生物目录》第二类（示例冠状病毒）" in row.hint
    assert "《传染病防治法》乙类（示例冠状病毒）" in row.hint
    assert "和《传染病防治法》乙类" in row.hint
    assert "不构成新序列风险分级或病原确证" in row.hint
    # Never claim the whole family is classified as a category.
    assert "Coronaviridae 在《" not in row.hint
    assert "属第二类" not in row.hint


def test_catalog_category_only_and_statutory_only_hints_keep_available_dimensions():
    taxonomy = {
        "cat_tip": _tax("cat_tip", family="CatOnlyviridae", genus="Catvirus"),
        "stat_tip": _tax("stat_tip", family="StatOnlyviridae", genus="Statvirus"),
    }
    catalog = [
        PathogenEntry(
            seq="1", name_cn="目录病毒", name_en="Catalog virus",
            family_cn="目录科", family_en="CatOnlyviridae",
            category="第三类", statutory="", covered="yes",
        ),
        PathogenEntry(
            seq="2", name_cn="法定病毒", name_en="Statutory virus",
            family_cn="法定科", family_en="StatOnlyviridae",
            category="", statutory="丙类", covered="yes",
        ),
    ]
    placements = [
        {"query_id": "Q1", "supergroup": "SG1", "nn1": "cat_tip", "nn1_dist": "0.1"},
        {"query_id": "Q2", "supergroup": "SG1", "nn1": "stat_tip", "nn1_dist": "0.2"},
    ]

    rows = build_risk_rows(placements, taxonomy, catalog)
    cat_hint = rows[0].hint
    stat_hint = rows[1].hint

    assert "《病原微生物目录》第三类" in cat_hint
    assert "《传染病防治法》" not in cat_hint
    assert "《传染病防治法》丙类" in stat_hint
    assert "《病原微生物目录》" not in stat_hint
    for hint in (cat_hint, stat_hint):
        assert "有名录收录成员涉及" in hint
        assert "不构成新序列风险分级或病原确证" in hint


def test_missing_nn1_taxonomy_keeps_row_evidence_and_skips_catalog_matching():
    taxonomy = {
        "known_human": _tax("known_human", family="Coronaviridae", genus="Betacoronavirus",
                            human=True, vertebrate=True),
    }
    # An empty-family catalog row must not match an unidentified nn1 family.
    catalog = [
        PathogenEntry(
            seq="1", name_cn="空科病毒", name_en="Empty family virus",
            family_cn="", family_en="",
            category="第一类", statutory="甲类", covered="no",
        )
    ]
    placements = [
        {"query_id": "Q1", "supergroup": "SG1", "nn1": "unknown_primary", "nn1_dist": "0.42",
         "nn2": "known_human", "nn2_dist": "0.90"},
        {"query_id": "Q2", "supergroup": "SG1", "nn1": "unknown_only", "nn1_dist": "0.50"},
    ]

    rows = build_risk_rows(placements, taxonomy, catalog)
    assert len(rows) == 2

    first, second = rows
    assert first.nn1 == "unknown_primary"
    assert first.nn1_dist == pytest.approx(0.42)
    assert first.family == ""
    assert first.genus == ""
    assert first.host_human is None
    assert first.host_vertebrate is None
    assert first.host_plant is None
    assert first.pathogen_category == ""
    assert first.statutory_category == ""
    assert first.neighbors[0].family == ""
    assert first.neighbors[0].human is None
    assert first.neighbors[1].reference_id == "known_human"
    assert first.neighbors[1].human is True
    assert "nn1 所属科未识别" in first.hint
    assert "未匹配到上述名录收录成员" in first.hint
    assert "不排除名录外或未识别病毒" in first.hint
    assert "常规监测" not in first.hint
    assert second.nn1 == "unknown_only"
    assert second.nn1_dist == pytest.approx(0.50)
    assert second.pathogen_category == ""
    assert len(second.neighbors) == 1
    assert second.neighbors[0].human is None
    assert second.neighbors[0].vertebrate is None

    summary = summarize_risk(rows)
    assert summary["human_nn1"] == 0
    assert summary["human_nn2_5"] == 1
    assert summary["vertebrate_only"] == 0
    assert summary["no_placed_label"] == 1
    assert summary["total"] == 2


def _legacy_row(
    query_id: str,
    *,
    host_human: bool | None = False,
    host_vertebrate: bool | None = False,
    host_plant: bool | None = False,
) -> RiskRow:
    return RiskRow(
        query_id=query_id,
        original_id="",
        contig_id="",
        supergroup="SG1",
        nn1="legacy_tip",
        nn1_dist=0.1,
        family="Legacyviridae",
        genus="Legacyvirus",
        host_human=host_human,
        host_vertebrate=host_vertebrate,
        host_plant=host_plant,
        pathogen_category="",
        pathogen_examples="",
        statutory_category="",
        statutory_examples="",
        hint="",
    )


def test_legacy_rows_without_neighbors_fall_back_to_legacy_host_flags():
    rows = [
        _legacy_row("Q1", host_human=True, host_vertebrate=True),
        _legacy_row("Q2", host_human=None, host_vertebrate=True),
        _legacy_row("Q3", host_human=None, host_vertebrate=None),
    ]

    summary = summarize_risk(rows)

    assert summary["human_nn1"] == 1
    assert summary["human_nn2_5"] == 0
    assert summary["vertebrate_only"] == 1
    assert summary["no_placed_label"] == 1
    assert (
        summary["human_nn1"] + summary["human_nn2_5"]
        + summary["vertebrate_only"] + summary["no_placed_label"]
    ) == summary["total"] == 3


@pytest.mark.parametrize(
    "bad_distance",
    [
        MISSING, None, "", "   ", "NA", "nan", "NaN", "inf",
        "-Infinity", "Infinity", float("nan"), float("inf"), float("-inf"),
    ],
    ids=[
        "missing", "none", "empty", "spaces", "text", "nan_lower", "nan_upper",
        "inf_lower", "negative_infinity", "infinity", "nan_float", "inf_float",
        "negative_inf_float",
    ],
)
def test_missing_or_nonfinite_distances_become_none_and_serialize_cleanly(bad_distance):
    taxonomy = {
        "nn1": _tax("nn1", family="FamA", genus="GenusA"),
        "nn2": _tax("nn2", family="FamB", genus="GenusB"),
    }
    placement = {"query_id": "Q1", "supergroup": "SG1", "nn1": "nn1", "nn2": "nn2"}
    if bad_distance is not MISSING:
        placement["nn1_dist"] = bad_distance
        placement["nn2_dist"] = bad_distance

    row = build_risk_rows([placement], taxonomy, [])[0]

    assert row.nn1_dist is None
    assert row.neighbors[0].distance is None
    assert row.neighbors[1].distance is None

    encoded = json.dumps(row.to_dict(), ensure_ascii=False, allow_nan=False)
    assert "NaN" not in encoded
    assert "Infinity" not in encoded
    decoded = json.loads(encoded)
    assert decoded["nn1_dist"] is None
    assert [neighbor["distance"] for neighbor in decoded["neighbors"]] == [None, None]
