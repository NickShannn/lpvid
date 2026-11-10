import csv
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest

import cdc_pipeline.stages.report as report_stage
from cdc_pipeline.pipeline import StageResult
from cdc_pipeline.risk_assessment import NeighborEvidence, RiskRow
from cdc_pipeline.state import RunState
from cdc_pipeline.stages.common import read_tsv, write_tsv
from cdc_pipeline.stages.report import build_report, plot_bars, summary_markdown
from cdc_pipeline.visualization import build_tree_html


KNOWN_NN1 = "mono-chu-like_YP_009337847.1_Hubei_qinvirus-like_virus_1"
ROOT_WHITELIST = {
    "main_results_report_zh.html",
    "main_results_report_en.html",
    "phylogenetic_explorer.html",
    "README.md",
    "tables",
    "sequences",
    "trees",
}
FILE_SUFFIXES = {".html", ".md", ".tsv", ".json", ".faa", ".fna", ".nwk", ".svg", ".png"}


class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.rows.append(dict(attrs))


def _file_links(document: Path):
    parser = _Links()
    parser.feed(document.read_text(encoding="utf-8"))
    for attrs in parser.rows:
        href = attrs.get("href", "")
        parsed = urlsplit(href)
        if not parsed.scheme and not parsed.netloc and Path(parsed.path).suffix in FILE_SUFFIXES:
            yield attrs, parsed


def _make_context(tmp_path):
    run_dir = tmp_path / "run"
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    stage_dir = run_dir / "step08_report"
    stage_dir.mkdir(parents=True)
    results = run_dir / "results"
    trees = results / "trees"
    trees.mkdir(parents=True)
    tree = trees / "SG_TEST.nwk"
    tree.write_text("(reference:0.1);\n", encoding="utf-8")

    def table(name, fields, rows):
        return write_tsv(upstream / name, fields, rows)

    representative_map = table(
        "representative_map.tsv",
        ["query_id", "representative_id", "aa_len", "motif_class_count", "motif_classes", "is_representative"],
        [{"query_id": "Q1", "representative_id": "Q1", "aa_len": 240, "motif_class_count": 2, "motif_classes": "A;C", "is_representative": 1}],
    )
    orf_map = table(
        "orf_contig_map.tsv",
        ["query_id", "orf_id", "contig_id"],
        [
            {"query_id": "Q1", "orf_id": "original-orf-1", "contig_id": "contig-1"},
            {"query_id": "UNRELATED", "orf_id": "do-not-load", "contig_id": "other-contig"},
        ],
    )
    placements = table(
        "placements.tsv",
        ["query_id", "supergroup", "nn1", "nn1_dist", "nn2", "nn2_dist", "nn3", "nn3_dist", "nn4", "nn4_dist", "nn5", "nn5_dist"],
        [{"query_id": "Q1", "supergroup": "SG_TEST", "nn1": KNOWN_NN1, "nn1_dist": "0.125", "nn3": "UNKNOWN_REFERENCE", "nn3_dist": ""}],
    )
    data = {
        "detect_counts": {"main_rdrp_positive_orfs": 1, "positive_contigs": 1},
        "representative_counts": {"representatives": 1},
        "supergroup_counts": {"SG_TEST": 1, "representatives": 1},
        "phylo_counts": {"placed": 1},
        "prefilter_calls": str(table("prefilter_calls.tsv", ["query_id", "analysis_set"], [{"query_id": "Q1", "analysis_set": "main"}])),
        "orf_predictions": str(table("orf_predictions.tsv", ["query_id", "analysis_set", "call", "rdrp_prob"], [{"query_id": "Q1", "analysis_set": "main", "call": "positive", "rdrp_prob": 0.8}])),
        "contig_calls": str(table("contig_level_calls.tsv", ["contig_id", "call"], [{"contig_id": "contig-1", "call": "positive"}])),
        "positive_faa": str(upstream / "rdrp_positive.faa"),
        "positive_contigs_fna": str(upstream / "rdrp_positive_contigs.fna"),
        "representatives_faa": str(upstream / "representatives.faa"),
        "representative_map": str(representative_map),
        "motifs": str(table("motifs.tsv", ["query_id", "motif_classes"], [{"query_id": "Q1", "motif_classes": "A;C"}])),
        "supergroup_assignments": str(table("supergroup_assignments.tsv", ["query_id", "supergroup", "status", "classification_method"], [{"query_id": "Q1", "supergroup": "SG_TEST", "status": "assigned", "classification_method": "reference match"}])),
        "placements_tsv": str(placements),
        "placement_annotations": str(table("placement_annotations.tsv", ["query_id", "original_id", "supergroup", "nn1", "nn1_dist"], [{"query_id": "Q1", "original_id": "original-orf-1", "supergroup": "SG_TEST", "nn1": KNOWN_NN1, "nn1_dist": 0.125}])),
        "orf_contig_map": str(orf_map),
        "tree_files": {"SG_TEST": str(tree)},
    }
    Path(data["positive_faa"]).write_text(">Q1\nMPEPTIDE\n", encoding="utf-8")
    Path(data["positive_contigs_fna"]).write_text(">contig-1\nACGT\n", encoding="utf-8")
    Path(data["representatives_faa"]).write_text(">Q1\nMPEPTIDE\n", encoding="utf-8")
    visualization_dir = run_dir / "step07_visualization"
    visualization_dir.mkdir()
    tree_html = build_tree_html([], visualization_dir / "tree_explorer.html", notice="fixture tree explorer")
    data["tree_html"] = str(tree_html)

    config = {"report": {"formats": [], "dpi": 80}, "project": {"display_name": "CDC fixture"}}
    state = RunState(run_dir)
    state.start({"route": "proteins", "sample": "fixture", "inputs": {}}, config)
    from types import SimpleNamespace
    return SimpleNamespace(
        route="proteins", sample="fixture", inputs={}, config=config,
        run_dir=run_dir, state=state, data=data, stage_dir=stage_dir,
    )


def test_plot_bars_writes_svg_and_png_for_multiple_categories(tmp_path):
    config = {"report": {"formats": ["svg", "png"], "dpi": 80}}
    outputs = plot_bars(tmp_path / "counts", ["SG1", "SG2"], [2, 1], "Representatives", "Count", config)
    assert len(outputs) == 2
    assert all(path.is_file() and path.stat().st_size > 0 for path in outputs)
    assert "Representatives" in (tmp_path / "counts.svg").read_text()


def test_single_category_uses_table_instead_of_one_bar_chart(tmp_path):
    config = {"report": {"formats": ["svg", "png"], "dpi": 80}}
    assert plot_bars(tmp_path / "single", ["SG1"], [1], "Title", "Count", config) == []


def test_summary_contains_scientific_boundaries_and_links():
    summary = {
        "run_id": "run1",
        "route": "proteins",
        "sample": "sample",
        "counts": {"positive_orfs": 2, "representatives": 1, "placed": 1},
        "warnings": ["1 条代表序列的骨架树正在建设中，暂未纳入系统发育放置"],
        "links": {"placements": "tables/placements.tsv", "tree": "phylogenetic_explorer.html"},
        "stages": {"visualization": {"status": "done"}},
    }
    text = summary_markdown(summary)
    assert "系统发育放置" in text
    assert "placements.tsv" in text
    assert "骨架树正在建设中" in text
    assert "使用边界" in text


def test_build_report_creates_categorized_delivery_and_bilingual_readme(tmp_path):
    context = _make_context(tmp_path)
    result = build_report(context)
    results = context.run_dir / "results"

    assert {path.name for path in results.iterdir()} == ROOT_WHITELIST
    assert (results / "tables" / "orf_predictions.tsv").is_file()
    assert (results / "tables" / "placement_annotations.tsv").is_file()
    assert (results / "tables" / "risk_warning.tsv").is_file()
    assert (results / "tables" / "neighbor_evidence.tsv").is_file()
    assert (results / "tables" / "run_summary.json").is_file()
    assert (results / "tables" / "stage_metrics.tsv").is_file()
    assert (results / "sequences" / "rdrp_positive_orfs.faa").is_file()
    assert (results / "sequences" / "rdrp_positive_contigs.fna").is_file()
    assert (results / "sequences" / "representative_orfs.faa").is_file()
    assert (results / "trees" / "SG_TEST.nwk").is_file()
    assert not (results / "report").exists()
    assert not (results / "summary.md").exists()
    assert not (results / "figures").exists()
    assert (context.stage_dir / "summary.md").is_file()
    assert (context.stage_dir / "summary.json").is_file()

    readme = (results / "README.md").read_text(encoding="utf-8")
    assert "文件索引" in readme and "File index" in readme
    for forbidden in ("Luca" + "Prot", "Luca" + "Virus", "Luca" + "Host"):
        assert forbidden not in readme
    delivered = [path for path in results.rglob("*") if path.is_file()]
    for path in delivered:
        assert path.relative_to(results).as_posix() in readme
    assert result["summary"]["supergroup_stats"]
    assert all(not link.startswith("../") for link in result["summary"]["links"].values())


def test_build_report_delivers_only_current_schema_tables(tmp_path):
    context = _make_context(tmp_path)
    result = build_report(context)
    results = context.run_dir / "results"

    delivered_tables = {path.name for path in (results / "tables").iterdir()}
    assert "abundance.tsv" not in delivered_tables  # proteins 路线没有丰度上游
    run_summary = json.loads((results / "tables" / "run_summary.json").read_text(encoding="utf-8"))
    for key in run_summary:
        assert "host" not in key.lower(), key
    assert "abundance_rows" in run_summary
    assert result["summary"]["abundance_rows"] == []


def test_build_report_readme_and_html_file_links_resolve_safely(tmp_path):
    context = _make_context(tmp_path)
    result = build_report(context)
    results = context.run_dir / "results"
    documents = [
        results / "README.md",
        result["report_zh"],
        result["report_en"],
    ]
    for document in documents:
        links = list(_file_links(document))
        assert links, document
        for attrs, parsed in links:
            assert attrs.get("target") == "_blank", (document, attrs)
            assert "noopener" in attrs.get("rel", "").split(), (document, attrs)
            assert Path(parsed.path).name != "tree_explorer.html"
            target = (document.parent / unquote(parsed.path)).resolve()
            target.relative_to(results.resolve())
            assert target.is_file(), (document, target)


def test_build_report_writes_structured_summary_and_neighbor_long_table(tmp_path):
    context = _make_context(tmp_path)
    result = build_report(context)
    results = context.run_dir / "results"
    summary = json.loads((results / "tables" / "run_summary.json").read_text(encoding="utf-8"))

    assert isinstance(summary["risk_rows"], list)
    row = summary["risk_rows"][0]
    assert isinstance(row, dict)
    assert row["original_id"] == "original-orf-1"
    assert row["contig_id"] == "contig-1"
    assert isinstance(row["neighbors"], list)
    assert [neighbor["rank"] for neighbor in row["neighbors"]] == [1, 3]
    unknown = row["neighbors"][1]
    assert unknown["reference_id"] == "UNKNOWN_REFERENCE"
    assert unknown["distance"] is None
    assert unknown["family"] == unknown["genus"] == ""
    assert unknown["human"] is None
    assert unknown["vertebrate"] is None
    assert unknown["plant"] is None
    assert isinstance(result["summary"]["risk_rows"][0], RiskRow)

    neighbors = read_tsv(results / "tables" / "neighbor_evidence.tsv")
    assert list(neighbors[0]) == ["query_id", "rank", "reference", "distance", "family", "genus", "human", "vertebrate", "plant"]
    assert [(row["query_id"], row["rank"]) for row in neighbors] == [("Q1", "1"), ("Q1", "3")]
    assert neighbors[1]["reference"] == "UNKNOWN_REFERENCE"
    assert neighbors[1]["distance"] == ""
    assert neighbors[1]["human"] == neighbors[1]["vertebrate"] == neighbors[1]["plant"] == ""
    risk = read_tsv(results / "tables" / "risk_warning.tsv")
    assert risk[0]["original_id"] == "original-orf-1"
    assert risk[0]["contig_id"] == "contig-1"


def test_delivery_run_summary_excludes_state_audit_but_step_summary_keeps_it(tmp_path):
    context = _make_context(tmp_path)
    stage_artifacts = {
        "tree_html": str(context.data["tree_html"]),
        "annotations": str(context.data["placement_annotations"]),
    }
    context.state.begin_stage("visualization", "visualization-signature")
    context.state.complete_stage(
        "visualization",
        "visualization-signature",
        stage_artifacts,
        status="done",
        summary={"bundle_count": 1, "query_count": 1},
        context_data={
            "tree_html": str(context.data["tree_html"]),
            "placement_annotations": str(context.data["placement_annotations"]),
        },
    )

    result = build_report(context)
    results = context.run_dir / "results"
    delivered = json.loads((results / "tables" / "run_summary.json").read_text(encoding="utf-8"))
    step_summary = json.loads((context.stage_dir / "summary.json").read_text(encoding="utf-8"))

    delivered_text = json.dumps(delivered, ensure_ascii=False, default=str)
    assert "context_data" not in delivered_text
    assert "artifacts" not in delivered_text
    assert str(context.run_dir.resolve()) not in delivered_text
    for internal in ("stages", "figures", "_run_dir", "_abundance_rows"):
        assert internal not in delivered
    for essential in (
        "counts", "risk_rows", "risk_summary", "supergroup_stats", "abundance_rows", "funnel",
        "stage_metrics", "timings", "total_minutes", "supergroup_outcome",
        "links", "warnings", "warnings_en",
    ):
        assert essential in delivered, essential
    assert delivered["run_id"] == "run"

    assert step_summary["_run_dir"] == str(context.run_dir)
    assert "figures" in step_summary
    visualization_stage = step_summary["stages"]["visualization"]
    assert visualization_stage["summary"]["query_count"] == 1
    assert visualization_stage["artifacts"] == stage_artifacts
    assert result["summary"]["stages"]["visualization"]["summary"]["query_count"] == 1


def test_build_report_scans_orf_contig_map_once_for_risk_annotations(tmp_path, monkeypatch):
    context = _make_context(tmp_path)
    # Prevent the independent HTML enrichment pass from reading the map so this
    # instrument measures the risk-building pass itself.
    context.data["representative_map"] = str(write_tsv(
        tmp_path / "empty_representatives.tsv",
        ["query_id", "representative_id", "aa_len", "is_representative"], [],
    ))
    mapping_path = str(Path(context.data["orf_contig_map"]).resolve())
    real_dict_reader = csv.DictReader
    scans = []

    def counting_dict_reader(handle, *args, **kwargs):
        if str(Path(getattr(handle, "name", "")).resolve()) == mapping_path:
            scans.append(handle.name)
        return real_dict_reader(handle, *args, **kwargs)

    monkeypatch.setattr(csv, "DictReader", counting_dict_reader)
    result = build_report(context)
    assert len(scans) == 1
    row = result["summary"]["risk_rows"][0]
    assert row.original_id == "original-orf-1"
    assert row.contig_id == "contig-1"


def test_build_report_repeat_cleans_only_known_legacy_delivery_entries(tmp_path):
    context = _make_context(tmp_path)
    build_report(context)
    results = context.run_dir / "results"
    legacy_table = results / "orf_predictions.tsv"
    legacy_sequence = results / "rdrp_positive.faa"
    legacy_report = results / "report" / "main_results_report_zh.html"
    legacy_figure = results / "report" / "figures" / "workflow_overview.png"
    unknown_root = results / "user_keep.txt"
    unknown_report = results / "report" / "user_keep.txt"
    legacy_figure.parent.mkdir(parents=True, exist_ok=True)
    for path in (legacy_table, legacy_sequence, legacy_report, legacy_figure, unknown_root, unknown_report):
        path.write_text("keep" if "user_keep" in path.name else "legacy", encoding="utf-8")

    build_report(context)

    assert not legacy_table.exists()
    assert not legacy_sequence.exists()
    assert not legacy_report.exists()
    assert not legacy_figure.exists()
    assert unknown_root.read_text(encoding="utf-8") == "keep"
    assert unknown_report.read_text(encoding="utf-8") == "keep"


def _empty_tsv(path):
    rows = read_tsv(path)
    assert rows
    write_tsv(path, list(rows[0]), [])


def test_build_report_rerun_with_empty_details_removes_stale_detail_deliveries(tmp_path):
    context = _make_context(tmp_path)
    build_report(context)
    results = context.run_dir / "results"
    assert (context.stage_dir / "supergroup_details.tsv").is_file()
    assert (context.stage_dir / "placement_details.tsv").is_file()
    assert (results / "tables" / "supergroup_details.tsv").is_file()
    assert (results / "tables" / "placement_details.tsv").is_file()

    # 上游从“有明细”变为“空明细”：本轮的报告不得再复制或索引上一轮明细。
    _empty_tsv(context.data["supergroup_assignments"])
    _empty_tsv(context.data["placements_tsv"])
    result = build_report(context)

    assert not (context.stage_dir / "supergroup_details.tsv").exists()
    assert not (context.stage_dir / "placement_details.tsv").exists()
    assert not (results / "tables" / "supergroup_details.tsv").exists()
    assert not (results / "tables" / "placement_details.tsv").exists()
    assert "supergroup_details" not in result["result_files"]
    assert "placement_details" not in result["result_files"]
    assert "supergroup_details.tsv" not in result["summary"]["links"]
    assert "placement_details.tsv" not in result["summary"]["links"]
    readme = (results / "README.md").read_text(encoding="utf-8")
    assert "supergroup_details.tsv" not in readme
    assert "placement_details.tsv" not in readme


def test_build_report_writes_complete_delivery_manifest_last(tmp_path, monkeypatch):
    context = _make_context(tmp_path)
    writes = []
    real_atomic_write_json = report_stage.atomic_write_json

    def recording_atomic_write_json(path, payload):
        writes.append(Path(path))
        return real_atomic_write_json(path, payload)

    monkeypatch.setattr(report_stage, "atomic_write_json", recording_atomic_write_json)
    result = build_report(context)
    results = context.run_dir / "results"
    manifest_path = results / "tables" / "delivery_manifest.json"

    assert writes[-1].resolve() == manifest_path.resolve()
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["status"] == "complete"
    assert payload["format_version"] == 1
    assert payload["run_id"] == "run"
    assert payload["files"] == sorted(payload["files"])
    expected = {
        "main_results_report_zh.html",
        "main_results_report_en.html",
        "phylogenetic_explorer.html",
        "README.md",
        "tables/run_summary.json",
        "tables/delivery_manifest.json",
        "tables/placement_annotations.tsv",
    }
    assert expected <= set(payload["files"])
    for relative in payload["files"]:
        assert not Path(relative).is_absolute()
        assert (results / relative).is_file(), relative
    assert result["result_files"]["delivery_manifest"] == manifest_path
    assert result["summary"]["links"]["delivery_manifest.json"] == "tables/delivery_manifest.json"
    assert "tables/delivery_manifest.json" in (results / "README.md").read_text(encoding="utf-8")


def test_build_report_failure_clears_sentinel_without_cleaning_old_delivery(tmp_path, monkeypatch):
    context = _make_context(tmp_path)
    build_report(context)
    results = context.run_dir / "results"
    manifest_path = results / "tables" / "delivery_manifest.json"
    assert manifest_path.is_file()
    old_summary = results / "tables" / "run_summary.json"
    assert old_summary.is_file()
    legacy_root = results / "orf_predictions.tsv"
    legacy_root.write_text("legacy-root", encoding="utf-8")
    legacy_report = results / "report" / "main_results_report_zh.html"
    legacy_report.parent.mkdir(parents=True, exist_ok=True)
    legacy_report.write_text("legacy-report", encoding="utf-8")

    real_build_html_report = report_stage.build_html_report

    def flaky_build_html_report(context, output_path, summary=None, lang="zh"):
        if lang == "en":
            raise RuntimeError("english report failed")
        return real_build_html_report(context, output_path, summary, lang=lang)

    monkeypatch.setattr(report_stage, "build_html_report", flaky_build_html_report)
    with pytest.raises(RuntimeError, match="english report failed"):
        build_report(context)

    assert not manifest_path.exists()
    assert not list((results / "tables").glob("delivery_manifest.json*"))
    assert old_summary.is_file()
    assert legacy_root.read_text(encoding="utf-8") == "legacy-root"
    assert legacy_report.read_text(encoding="utf-8") == "legacy-report"


def _artifact_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from _artifact_strings(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _artifact_strings(child)


def test_report_run_artifacts_are_complete_and_allow_stage_skip(tmp_path):
    context = _make_context(tmp_path)
    result = report_stage.run(context)
    assert isinstance(result, StageResult)
    signature = "report-signature"
    context.state.begin_stage("report", signature)
    context.state.complete_stage(
        "report", signature, result.artifacts,
        status=result.status, summary=result.summary, context_data=result.data,
    )

    assert context.state._artifacts_complete(result.artifacts)
    assert context.state.should_skip("report", signature)
    relative_artifacts = set(_artifact_strings(result.artifacts))
    expected = {
        "results/main_results_report_zh.html",
        "results/main_results_report_en.html",
        "results/README.md",
        "results/tables/run_summary.json",
        "results/tables/delivery_manifest.json",
        "results/tables/placement_annotations.tsv",
    }
    assert expected <= relative_artifacts
    for relative in expected:
        assert (context.run_dir / relative).is_file(), relative

    (context.run_dir / "results" / "README.md").unlink()
    assert not context.state._artifacts_complete(result.artifacts)
    assert not context.state.should_skip("report", signature)


def test_risk_tables_export_none_as_blank(tmp_path):
    writer = getattr(report_stage, "_write_risk_tables", None)
    assert callable(writer)
    row = RiskRow(
        query_id="Q_UNKNOWN", original_id="", contig_id="", supergroup="SG",
        nn1="UNKNOWN", nn1_dist=None, family="", genus="",
        host_human=None, host_vertebrate=None, host_plant=None,
        pathogen_category="", pathogen_examples="", statutory_category="",
        statutory_examples="", hint="unknown",
        neighbors=[NeighborEvidence(2, "UNKNOWN", None, "", "", None, None, None)],
    )
    writer(tmp_path, [row])
    risk = read_tsv(tmp_path / "risk_warning.tsv")[0]
    assert risk["nn1_dist"] == ""
    assert risk["host_human"] == risk["host_vertebrate"] == risk["host_plant"] == ""
    neighbor = read_tsv(tmp_path / "neighbor_evidence.tsv")[0]
    assert neighbor["rank"] == "2"
    assert neighbor["distance"] == ""
    assert neighbor["human"] == neighbor["vertebrate"] == neighbor["plant"] == ""


def test_build_report_rejects_stale_tree_contract_and_asks_for_visualization_rerun(tmp_path):
    context = _make_context(tmp_path)
    tree_html = Path(context.data["tree_html"])
    tree_html.write_text(
        tree_html.read_text(encoding="utf-8") + "\n<!-- HOST_CLASSES 宿主大类 -->",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="重跑 visualization"):
        build_report(context)
    assert not (context.run_dir / "results" / "tables" / "delivery_manifest.json").exists()


def test_rerun_report_only_refuses_when_visualization_signature_is_stale(tmp_path):
    from scripts.rerun_report_only import rerun

    context = _make_context(tmp_path)
    context.state.begin_stage("visualization", "stale-signature")
    context.state.complete_stage(
        "visualization",
        "stale-signature",
        {"tree_html": str(context.data["tree_html"])},
        status="done",
        context_data={"tree_html": str(context.data["tree_html"])},
    )
    with pytest.raises(RuntimeError, match="visualization"):
        rerun(context.run_dir)


def test_refresh_align_removes_visualization_so_it_rebuilds(tmp_path):
    from scripts.refresh_stage_signatures import force_align

    context = _make_context(tmp_path)
    context.state.begin_stage("visualization", "aligned-signature")
    context.state.complete_stage(
        "visualization",
        "aligned-signature",
        {"tree_html": str(context.data["tree_html"])},
        status="done",
        context_data={"tree_html": str(context.data["tree_html"])},
    )
    force_align(context.run_dir)
    state = type(context.state)(context.run_dir)
    assert "visualization" not in state.snapshot["stages"]
    assert "report" not in state.snapshot["stages"]


def _custom_tree_with_empty_original_id(tmp_path):
    bundle = {
        "supergroup": "SG_TEST",
        "kind": "virid",
        "backbone_newick": "(PBV_TIP:0.1,OTHER_TIP:0.2);",
        "n_backbone": 2,
        "n_query": 1,
        "placements": [
            {
                "query_id": "Q_PBV",
                "original_id": "",
                "supergroup": "SG_TEST",
                "neighbors": [["PBV_TIP", 0.1]],
            }
        ],
        "truth": None,
        "metrics": {},
        "taxonomy": {"PBV_TIP": {"Family": "Picobirnaviridae"}, "OTHER_TIP": {"Family": "Otherviridae"}},
        "host_marks": {},
        "ranks": ["Family"],
        "value_counts": {"Family": {"Picobirnaviridae": 1, "Otherviridae": 1}},
        "provenance": {},
    }
    output = tmp_path / "proteins_tree.html"
    return build_tree_html([bundle], output, notice="proteins fixture")


def test_build_report_allows_proteins_placement_without_orf_map(tmp_path, monkeypatch):
    context = _make_context(tmp_path)
    context.data.pop("orf_contig_map", None)
    context.data["tree_html"] = str(_custom_tree_with_empty_original_id(tmp_path))

    result = build_report(context)
    results = context.run_dir / "results"
    assert (results / "tables" / "delivery_manifest.json").is_file()
    assert (results / "phylogenetic_explorer.html").is_file()
    assert result["summary"]["abundance_rows"] == []


def test_backbone_pending_warning_merges_counts_and_hides_zero():
    from cdc_pipeline.stages.report import backbone_pending_warning

    assert backbone_pending_warning(0, 0) is None
    zh, en = backbone_pending_warning(11, 10)
    assert zh == "21 条代表序列的骨架树正在建设中，暂未纳入系统发育放置"
    assert "21" in en and "backbone tree" in en and "nn1" not in en
    zh_only, en_only = backbone_pending_warning(0, 4)
    assert zh_only == "4 条代表序列的骨架树正在建设中，暂未纳入系统发育放置"
    assert "4" in en_only


def test_public_hint_replaces_latin_abbreviations():
    from cdc_pipeline.stages.report import _public_hint

    raw = (
        "Coronaviridae（nn1 所属科）的已知成员中，有名录收录成员涉及"
        "《病原微生物目录》第二类（示例冠状病毒）；需进一步核查，不构成新序列风险分级或病原确证。"
    )
    text = _public_hint(raw)
    assert "nn1" not in text
    assert "第 1 近邻所属科" in text
    assert _public_hint("nn1 所属科未识别；") == "第 1 近邻所属科未识别；"
    assert _public_hint("") == ""


def test_build_report_delivers_sanitized_public_hints(tmp_path):
    context = _make_context(tmp_path)
    build_report(context)
    results = context.run_dir / "results"
    risk = read_tsv(results / "tables" / "risk_warning.tsv")
    assert risk
    assert "nn1" not in risk[0]["hint"]
    delivered = json.loads((results / "tables" / "run_summary.json").read_text(encoding="utf-8"))
    for row in delivered["risk_rows"]:
        assert "nn1" not in row["hint"]
