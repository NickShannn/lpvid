import csv

from cdc_pipeline.placement import placements_from_matrix, summarize_distances


def test_summarize_distances_excludes_query_and_returns_topk():
    neighbors = summarize_distances(
        {"ref1": 2.0, "ref2": 4.0, "query": 0.0}, backbone_ids={"ref1", "ref2"}, top_k=2
    )
    assert [n.tip_id for n in neighbors] == ["ref1", "ref2"]
    assert [n.distance for n in neighbors] == [2.0, 4.0]


def test_placements_from_full_matrix(tmp_path):
    path = tmp_path / "matrix.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "ref1", "ref2", "Q1"])
        writer.writerow(["ref1", 0, 5, 2])
        writer.writerow(["ref2", 5, 0, 4])
        writer.writerow(["Q1", 2, 4, 0])
    placements = placements_from_matrix(path, {"ref1", "ref2"}, top_k=2, supergroup="SG")
    assert len(placements) == 1
    assert placements[0].query_id == "Q1"
    assert placements[0].supergroup == "SG"
    assert [n.tip_id for n in placements[0].neighbors] == ["ref1", "ref2"]
    assert [n.distance for n in placements[0].neighbors] == [2.0, 4.0]
