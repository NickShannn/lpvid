import json
import re
from pathlib import Path
from types import SimpleNamespace

from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.placement import Neighbor, Placement
from cdc_pipeline.visualization import build_bundle, build_tree_html, run
from cdc_pipeline.stages.common import write_tsv


def make_bundle(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:0.1,B:0.2):0.3,C:0.4);\n")
    taxonomy = tmp_path / "taxonomy.csv"
    taxonomy.write_text(
        "Cluster_sseqid,Realm,Class,Family,Human,Vertebrates,Plant\n"
        "A,Riboviria,ClassA,Fam,1,1,0\n"
        "B,Riboviria,ClassA,Fam,0,1,0\n"
        "C,Riboviria,ClassB,Other,0,0,1\n"
    )
    placement = Placement("Q1", "TEST_SG", [Neighbor("A", 0.2), Neighbor("B", 0.4)])
    return build_bundle("TEST_SG", tree, taxonomy, [placement], {"run_id": "test"})


def test_bundle_has_real_query_schema_and_taxonomy(tmp_path):
    bundle = make_bundle(tmp_path)
    assert bundle["kind"] == "virid"
    assert bundle["n_backbone"] == 3
    assert bundle["n_query"] == len(bundle["placements"]) == 1
    assert bundle["truth"] is None
    assert bundle["host_marks"]["A"] == [1, 1, 0]
    assert bundle["value_counts"]["Family"]["Fam"] == 2


def test_bundle_placements_carry_query_and_neighbors_only(tmp_path):
    bundle = make_bundle(tmp_path)
    row = bundle["placements"][0]
    assert set(row) == {"query_id", "original_id", "supergroup", "neighbors"}
    assert row["neighbors"] == [["A", 0.2], ["B", 0.4]]


def test_html_contains_only_current_run_bundle(tmp_path):
    output = tmp_path / "tree.html"
    build_tree_html([make_bundle(tmp_path)], output)
    text = output.read_text()
    assert "TEST_SG" in text
    assert "/*__PHYLO_CORE__*/" not in text
    assert "/*__VIRID_BUNDLES__*/" not in text
    assert "bh pol genes" not in text.lower()


def test_visualization_run_delivers_annotations_under_tables(tmp_path, monkeypatch):
    stage_dir = tmp_path / "step08_visualization"
    stage_dir.mkdir()
    placements = write_tsv(stage_dir / "placements.tsv", ["query_id", "supergroup"], [])
    context = SimpleNamespace(
        run_dir=tmp_path,
        stage_dir=stage_dir,
        config={"phylo": {"top_k": 5}},
        data={"placements_tsv": str(placements), "phylo_matrices": {}},
    )
    monkeypatch.setattr("cdc_pipeline.visualization.load_supergroup_map", lambda config: {})

    result = run(context)

    annotation = tmp_path / "results" / "tables" / "placement_annotations.tsv"
    assert annotation.is_file()
    assert not (tmp_path / "results" / "placement_annotations.tsv").exists()
    assert result.data["placement_annotations"] == str(annotation)
    assert (tmp_path / "results" / "trees").is_dir()


def _bundle_json(tree_html):
    match = re.search(
        r'<script id="virid-bundles" type="application/json">\s*(.*?)\s*</script>',
        tree_html.read_text(encoding="utf-8"),
        re.S,
    )
    assert match
    return json.loads(match.group(1))


def _run_context(tmp_path):
    stage_dir = tmp_path / "step08_visualization"
    stage_dir.mkdir()
    tree = tmp_path / "TEST_SG.nwk"
    tree.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")
    taxonomy = tmp_path / "TEST_SG.csv"
    taxonomy.write_text(
        "Cluster_sseqid,Family,Human,Vertebrates,Plant\n"
        "A,Picobirnaviridae,0,0,0\n"
        "B,Otherviridae,1,1,0\n",
        encoding="utf-8",
    )
    placements = write_tsv(
        stage_dir / "placements.tsv",
        ["query_id", "supergroup", "nn1", "nn1_dist"],
        [
            {"query_id": "Q_A", "supergroup": "TEST_SG", "nn1": "A", "nn1_dist": 0.2},
            {"query_id": "Q_B", "supergroup": "TEST_SG", "nn1": "B", "nn1_dist": 0.3},
        ],
    )
    orf_map = write_tsv(
        tmp_path / "orf_contig_map.tsv",
        ["query_id", "orf_id", "contig_id"],
        [
            {"query_id": "Q_A", "orf_id": "o1", "contig_id": "c1"},
            {"query_id": "Q_B", "orf_id": "o2", "contig_id": "c2"},
        ],
    )
    return SimpleNamespace(
        run_dir=tmp_path,
        stage_dir=stage_dir,
        config={"phylo": {"top_k": 5}, "assets": {"phylo_root": str(tmp_path)}},
        data={
            "placements_tsv": str(placements),
            "phylo_matrices": {"TEST_SG": str(tmp_path / "matrix.tsv")},
            "orf_contig_map": str(orf_map),
            "_fake_asset": SimpleNamespace(tree=tree, taxonomy=taxonomy, adapter=tmp_path / "adapter"),
        },
    )


def _patch_assets(monkeypatch, context):
    monkeypatch.setattr("cdc_pipeline.visualization.load_supergroup_map", lambda config: {"TEST_SG": context.data["_fake_asset"]})
    monkeypatch.setattr("cdc_pipeline.supergroup_mapping.lucaphylo_to_asset_name", lambda name: name)


def test_visualization_run_embeds_original_ids_and_neighbors(tmp_path, monkeypatch):
    context = _run_context(tmp_path)
    _patch_assets(monkeypatch, context)

    result = run(context)

    data = _bundle_json(Path(result.data["tree_html"]))
    by_query = {row["query_id"]: row for row in data["TEST_SG"]["placements"]}
    assert by_query["Q_A"]["original_id"] == "o1"
    assert by_query["Q_A"]["neighbors"] == [["A", 0.2]]
    assert by_query["Q_B"]["original_id"] == "o2"
    for row in by_query.values():
        assert "host" not in row
        assert "host_suppression_reason" not in row


def test_tree_verifier_rejects_legacy_tokens_and_stale_host_payload(tmp_path):
    from scripts.verify_tree_explorer import check_file

    output = tmp_path / "tree.html"
    build_tree_html([make_bundle(tmp_path)], output)
    good_text = output.read_text(encoding="utf-8")

    output.write_text(good_text + "\n<!-- HOST_CLASSES hostClassForTip 宿主大类 -->", encoding="utf-8")
    problems = check_file(output)
    assert any("HOST_CLASSES" in problem for problem in problems)
    assert any("hostClassForTip" in problem for problem in problems)
    assert any("宿主大类" in problem for problem in problems)

    stale = good_text.replace(
        '"query_id": "Q1"', '"query_id": "Q1", "host": {"human": 0.9}', 1)
    output.write_text(stale, encoding="utf-8")
    problems = check_file(output)
    assert any("旧版宿主概率字段" in problem for problem in problems)


def test_tree_verifier_supports_optional_original_id_requirement(tmp_path):
    from scripts.verify_tree_explorer import check_file

    bundle = make_bundle(tmp_path)  # placement original_id 为空
    output = tmp_path / "tree.html"
    build_tree_html([bundle], output)
    assert any("original_id" in problem for problem in check_file(output))
    assert check_file(output, require_original_ids=False) == []

    output.write_text(output.read_text(encoding="utf-8") + "\n<!-- HOST_CLASSES -->", encoding="utf-8")
    assert any("HOST_CLASSES" in problem for problem in check_file(output, require_original_ids=False))


def test_tree_template_keeps_taxonomy_icons_and_drops_probability_payload():
    source = PROJECT_ROOT / "third_party" / "placement_view"
    core = (source / "phylo_core.js").read_text(encoding="utf-8")
    template = (source / "tree_explorer_template.html").read_text(encoding="utf-8")
    assert "HOST_ICONS" in template
    assert "host_marks" in template
    assert "HOST_CLASSES" not in core
    assert "HOST_LABEL" not in core
    assert "hostClassForTip" not in core
    assert "host_class" not in core
    assert "宿主大类" not in template
    assert "host_class" not in template
    assert "hostPct" not in template
    assert "hostRow" not in template
    assert "host_suppression_reason" not in template
