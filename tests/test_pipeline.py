import importlib
import json

import pytest

from cdc_pipeline.pipeline import _taxonomy_dependencies, stage_dependencies, stages_for
from cdc_pipeline.stage_metrics import extract_sequence_funnel, extract_stage_metrics, extract_supergroup_outcome
from cdc_pipeline.stages.common import write_tsv


def test_stage_routes_cover_the_complete_flow():
    assert stages_for("proteins") == [
        "prefilter", "detect", "representatives", "supergroup", "classifier", "phylo", "visualization", "report"
    ]
    assert stages_for("contigs") == [
        "orf", "prefilter", "detect", "representatives", "supergroup", "classifier", "phylo", "visualization", "report"
    ]
    reads = stages_for("reads")
    assert reads[0:2] == ["preprocess", "orf"]
    assert reads[-1] == "report"
    assert reads.index("abundance") < reads.index("visualization")
    for route in ("proteins", "contigs", "reads"):
        stages = stages_for(route)
        assert stages[-1] == "report"
        assert "classifier" in stages


def test_every_route_stage_is_importable_and_callable():
    for route in ("proteins", "contigs", "reads"):
        for name in stages_for(route):
            module = importlib.import_module(f"cdc_pipeline.stages.{name}")
            assert callable(module.run), name


@pytest.mark.skipif(not _taxonomy_dependencies(), reason="外部超群 taxonomy 资源包未配置")
def test_report_and_visualization_signatures_track_taxonomy_dependencies():
    report_names = [path.name for path in stage_dependencies("report")]
    assert "html_report.py" in report_names
    assert "verify_tree_explorer.py" in report_names
    assert any(name.endswith(".csv") and name != "china_pathogen_catalog.csv" for name in report_names)

    visualization_names = [path.name for path in stage_dependencies("visualization")]
    assert "risk_assessment.py" in visualization_names
    assert "tree_explorer_template.html" in visualization_names
    assert "phylo_core.js" in visualization_names
    assert any(name.endswith(".csv") for name in visualization_names)
    assert stage_dependencies("phylo") == []


def _metrics_run(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    source = tmp_path / "input.faa"
    source.write_text(">input\nMPEPTIDE\n", encoding="utf-8")
    (run_dir / "run_manifest.json").write_text(json.dumps({
        "route": "proteins", "inputs": {"input": str(source)},
    }), encoding="utf-8")
    (run_dir / "state.json").write_text(json.dumps({"stages": {}}), encoding="utf-8")
    return run_dir


def _prediction_rows(prefix, count):
    return [{
        "query_id": f"{prefix}{index}", "analysis_set": "main", "call": "positive",
    } for index in range(count)]


def test_metrics_support_categorized_delivery_paths(tmp_path):
    run_dir = _metrics_run(tmp_path)
    tables = run_dir / "results" / "tables"
    sequences = run_dir / "results" / "sequences"
    tables.mkdir(parents=True)
    sequences.mkdir(parents=True)
    write_tsv(tables / "orf_predictions.tsv", ["query_id", "analysis_set", "call"], _prediction_rows("Q", 2))
    (sequences / "representative_orfs.faa").write_text(">Q0\nMA\n>Q1\nMB\n", encoding="utf-8")
    write_tsv(tables / "supergroup_assignments.tsv", ["query_id", "supergroup", "status"], [
        {"query_id": "Q0", "supergroup": "SG", "status": "assigned"},
        {"query_id": "Q1", "supergroup": "SG", "status": "assigned"},
    ])
    write_tsv(tables / "placements.tsv", ["query_id", "nn1"], [
        {"query_id": "Q0", "nn1": "R0"}, {"query_id": "Q1", "nn1": "R1"},
    ])

    funnel = {row["key"]: row["count"] for row in extract_sequence_funnel(run_dir)}
    assert funnel["positive"] == 2
    assert funnel["representatives"] == 2
    assert funnel["phylo"] == 2
    assert "host" not in funnel
    assert extract_supergroup_outcome(run_dir)["representatives"] == 2


def test_stage_metrics_zero_counts_keep_business_units_when_funnel_rows_missing(tmp_path):
    run_dir = _metrics_run(tmp_path)
    (run_dir / "state.json").write_text(json.dumps({
        "stages": {
            "representatives": {"status": "not_run"},
            "phylo": {"status": "not_run"},
            "report": {"status": "running"},
        },
    }), encoding="utf-8")
    metrics = extract_stage_metrics(run_dir)
    assert metrics
    rows = {row["stage_key"]: row for row in metrics}
    for key in ("representatives", "phylo"):
        assert key in rows
        assert rows[key]["input_count"] == 0
        assert rows[key]["output_count"] == 0
        assert rows[key]["output_unit"] == "条代表蛋白"
    assert rows["representatives"]["input_unit"] == "条蛋白 ORFs"
    assert rows["phylo"]["input_unit"] == "条代表蛋白"
    assert rows["representatives"]["status"] == "not_run"
    assert rows["orf"]["input_unit"] == "条蛋白序列"
    assert rows["orf"]["output_unit"] == "条蛋白 ORFs"
    for key in ("prefilter", "detect"):
        assert rows[key]["input_unit"] == "条蛋白 ORFs"
        assert rows[key]["output_unit"] == "条蛋白 ORFs"
    assert rows["report"]["input_unit"] == "条代表蛋白"
    assert rows["report"]["output_unit"] == "份报告"
    for row in metrics:
        assert row["input_unit"], row["stage_key"]
        assert row["output_unit"], row["stage_key"]


def test_metrics_prefer_original_stage_artifacts_over_delivery_copies(tmp_path):
    run_dir = _metrics_run(tmp_path)
    detect = run_dir / "step03_detect"
    representatives = run_dir / "step04_representatives"
    supergroup = run_dir / "step05_supergroup"
    detect.mkdir(parents=True)
    representatives.mkdir(parents=True)
    supergroup.mkdir(parents=True)
    write_tsv(detect / "orf_predictions.tsv", ["query_id", "analysis_set", "call"], _prediction_rows("ORIGINAL", 1))
    (representatives / "representatives.faa").write_text(">ORIGINAL0\nMA\n", encoding="utf-8")
    write_tsv(supergroup / "supergroup_assignments.tsv", ["query_id", "supergroup", "status"], [
        {"query_id": "ORIGINAL0", "supergroup": "SG", "status": "assigned"},
    ])
    legacy = run_dir / "results"
    legacy.mkdir()
    write_tsv(legacy / "orf_predictions.tsv", ["query_id", "analysis_set", "call"], _prediction_rows("COPY", 3))
    (legacy / "representatives.faa").write_text(">C0\nMA\n>C1\nMB\n>C2\nMC\n", encoding="utf-8")
    write_tsv(legacy / "supergroup_assignments.tsv", ["query_id", "supergroup", "status"], [
        {"query_id": f"C{index}", "supergroup": "SG", "status": "assigned"} for index in range(3)
    ])

    funnel = {row["key"]: row["count"] for row in extract_sequence_funnel(run_dir)}
    assert funnel["positive"] == 1
    assert funnel["representatives"] == 1
    assert extract_supergroup_outcome(run_dir)["representatives"] == 1
