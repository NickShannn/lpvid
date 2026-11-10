import html
import json
import re
import shutil
import subprocess
from dataclasses import replace
from types import SimpleNamespace

import pytest

from cdc_pipeline.config import load_config
from cdc_pipeline.html_report import build_html_report, length_statistics
from cdc_pipeline.html_template_v2 import render_html
from cdc_pipeline.risk_assessment import NeighborEvidence, RiskRow, summarize_risk
from cdc_pipeline.stages.common import write_tsv
from scripts.verify_report_v2 import check_file, main as verify_main


SECTION_IDS = ["overview", "risk", "sequences", "statistics", "abundance", "workflow", "links"]
ZH_TITLES = ["总体摘要", "风险提示", "代表序列", "详细序列统计", "相对丰度", "分析流程概览", "结果文件"]
EN_TITLES = [
    "Overview",
    "Prior-based risk notices",
    "Representative sequences",
    "Detailed sequence statistics",
    "Relative abundance",
    "Workflow overview",
    "Result files",
]
FORBIDDEN_PRODUCT_NAMES = ("Luca" + "Prot", "Luca" + "Virus", "Luca" + "Host")


def report_data(text):
    match = re.search(r"const REPORT_DATA = (\{.*?\});", text, re.S)
    assert match
    return json.loads(match.group(1))


def body_without_code(text):
    text = re.sub(r"<script.*?</script>", "", text, flags=re.S)
    return re.sub(r"<style.*?</style>", "", text, flags=re.S)


def section_titles(text):
    return [html.unescape(re.sub(r"<[^>]+>", "", value)).strip() for value in re.findall(r"<h2[^>]*>(.*?)</h2>", text, re.S)]


def neighbor(rank, reference, *, distance=None, family="", human=None, vertebrate=None):
    return NeighborEvidence(
        rank=rank,
        reference_id=reference,
        distance=distance,
        family=family,
        genus="",
        human=human,
        vertebrate=vertebrate,
        plant=None,
    )


def risk_row(query, contig, supergroup, family, neighbors, *, category="", statutory="", distance=None):
    return RiskRow(
        query_id=query,
        original_id=f"orf-{query}",
        contig_id=contig,
        supergroup=supergroup,
        nn1=neighbors[0].reference_id,
        nn1_dist=distance,
        family=family,
        genus="",
        host_human=neighbors[0].human,
        host_vertebrate=neighbors[0].vertebrate,
        host_plant=None,
        pathogen_category=category,
        pathogen_examples="示例冠状病毒" if category else "",
        statutory_category=statutory,
        statutory_examples="示例法定传染病" if statutory else "",
        hint="中文提示不得复用到英文页面",
        neighbors=neighbors,
    )


def rich_fixture(route="reads"):
    rows = [
        {"query_id": "Q1", "original_orf_id": "orf-Q1", "contig_id": "c1", "aa_len": 210, "call": "positive", "rdrp_prob": 0.91, "supergroup": "SG1", "classification_method": "reference match", "classification_status": "assigned", "reference_id": "ORIGINAL_REF_1", "pident": 88.0, "evalue": 1e-8, "nn1": "plain_ref_1", "nn1_dist": None, "rpm": 12.5, "mean_depth": 3.0, "coverage_breadth_pct": 80.0},
        {"query_id": "Q2", "original_orf_id": "orf-Q2", "contig_id": "c2", "aa_len": 260, "call": "positive", "rdrp_prob": 0.82, "supergroup": "SG1", "classification_method": "reference match", "classification_status": "assigned", "reference_id": "ORIGINAL_REF_2", "pident": 80.0, "evalue": 1e-6, "nn1": "human_ref_1", "nn1_dist": 0.2, "rpm": 8.0, "mean_depth": 2.0, "coverage_breadth_pct": 70.0},
        {"query_id": "Q3", "original_orf_id": "orf-Q3", "contig_id": "c3", "aa_len": 320, "call": "positive", "rdrp_prob": 0.77, "supergroup": "SG2", "classification_method": "generic classification", "classification_status": "classified_no_tree", "reference_id": "", "pident": "", "evalue": "", "nn1": "vertebrate_ref_1", "nn1_dist": 0.4, "rpm": 4.0, "mean_depth": 1.0, "coverage_breadth_pct": 60.0},
        {"query_id": "Q4", "original_orf_id": "orf-Q4", "contig_id": "c4", "aa_len": 390, "call": "positive", "rdrp_prob": 0.66, "supergroup": "SG3", "classification_method": "reference match", "classification_status": "classified_pending_merge", "reference_id": "ORIGINAL_REF_4", "pident": 70.0, "evalue": 1e-4, "nn1": "unknown_ref_1", "nn1_dist": None, "rpm": 2.0, "mean_depth": 0.5, "coverage_breadth_pct": 50.0},
        {"query_id": "Q5", "original_orf_id": "orf-Q5", "contig_id": "c5", "aa_len": 180, "call": "positive", "rdrp_prob": 0.61, "supergroup": "", "classification_method": "", "classification_status": "unassigned", "reference_id": "", "pident": "", "evalue": "", "nn1": "", "nn1_dist": None, "rpm": "", "mean_depth": "", "coverage_breadth_pct": ""},
    ]
    risks = [
        risk_row(
            "Q1", "c1", "SG1", "ZetaMatchedviridae",
            [neighbor(1, "plain_ref_1", distance=None, family="ZetaMatchedviridae", human=False, vertebrate=False),
             neighbor(4, "human_ref_4", distance=1.4, family="Otherviridae", human=True, vertebrate=True)],
            category="第二类", statutory="乙类", distance=None,
        ),
        risk_row(
            "Q2", "c2", "SG1", "AlphaUnmatchedviridae",
            [neighbor(1, "human_ref_1", distance=0.2, family="AlphaUnmatchedviridae", human=True, vertebrate=True)],
            distance=0.2,
        ),
        risk_row(
            "Q3", "c3", "SG2", "VertebrateOnlyviridae",
            [neighbor(1, "vertebrate_ref_1", distance=0.4, family="VertebrateOnlyviridae", human=False, vertebrate=True)],
            distance=0.4,
        ),
        risk_row(
            "Q4", "c4", "SG3", "",
            [neighbor(1, "unknown_ref_1", distance=None, family="", human=None, vertebrate=None)],
            distance=None,
        ),
    ]
    risk_summary = summarize_risk(risks)
    # Keep the real summary, then exercise a missing neighbor rank in rendering.
    risks[3] = replace(risks[3], neighbors=[replace(risks[3].neighbors[0], rank=None)])
    stats = [
        {"supergroup": "SG1", "representatives": 2, "positive_orfs": 4, "positive_contigs": 2, "length_statistics": length_statistics(rows[:2]), "classification_methods": {"reference": 2}, "tree_statuses": {"assigned": 2}, "abundance": {"n_contigs": 2, "mapped_reads": 30, "RPM": 20.5, "mean_depth_weighted": 2.5, "coverage_breadth_pct": 75.0}},
        {"supergroup": "SG2", "representatives": 1, "positive_orfs": 2, "positive_contigs": 2, "length_statistics": length_statistics(rows[2:3]), "classification_methods": {"generic classification": 1}, "tree_statuses": {"classified_no_tree": 1}, "abundance": {"n_contigs": 1, "mapped_reads": 10, "RPM": 4.0, "mean_depth_weighted": 1.0, "coverage_breadth_pct": 60.0}},
        {"supergroup": "SG3", "representatives": 1, "positive_orfs": 1, "positive_contigs": 1, "length_statistics": length_statistics(rows[3:4]), "classification_methods": {"reference": 1}, "tree_statuses": {"classified_pending_merge": 1}, "abundance": {"n_contigs": 1, "mapped_reads": 5, "RPM": 2.0, "mean_depth_weighted": 0.5, "coverage_breadth_pct": 50.0}},
        {"supergroup": "Unassigned", "representatives": 1, "positive_orfs": 1, "positive_contigs": 1, "length_statistics": length_statistics(rows[4:5]), "classification_methods": {"Unassigned": 1}, "tree_statuses": {"unassigned": 1}, "abundance": {"n_contigs": 0, "mapped_reads": 0, "RPM": 0.0, "mean_depth_weighted": 0.0, "coverage_breadth_pct": 0.0}},
    ]
    summary = {
        "run_id": "synthetic-run",
        "route": route,
        "sample": "synthetic-sample",
        "counts": {"positive_orfs": 8, "representatives": 5, "classified": 4, "placed": 4, "pending_tree": 4, "unassigned": 1, "positive_contigs": 5},
        "risk_rows": risks,
        "risk_summary": risk_summary,
        "supergroup_stats": stats,
        "multi_supergroup_contig_count": 1,
        "multi_supergroup_contigs": [{"contig_id": "c2", "supergroups": ["SG1", "SG2"]}],
        "abundance_rows": [] if route == "proteins" else [
            {"contig_id": "c1", "length": 1000, "mapped_reads": 20, "rpm": 12.5, "mean_depth": 3.0, "coverage_breadth_pct": 80.0},
            {"contig_id": "c2", "length": 1200, "mapped_reads": 10, "rpm": 8.0, "mean_depth": 2.0, "coverage_breadth_pct": 70.0},
        ],
        "funnel": [
            {"key": "input", "stage": "核酸 contigs", "count": 20, "unit": "条核酸 contigs"},
            {"key": "orf", "stage": "蛋白 ORFs", "count": 12, "unit": "条蛋白 ORFs"},
            {"key": "prefilter", "stage": "进入模型", "count": 0, "unit": "条蛋白 ORFs"},
            {"key": "positive", "stage": "RdRP 阳性", "count": 8, "unit": "条蛋白 ORFs"},
            {"key": "representatives", "stage": "代表蛋白", "count": 5, "unit": "条代表蛋白"},
        ],
        "stage_metrics": [
            {"stage_key": "orf", "stage": "蛋白 ORF 预测", "status": "done", "input_count": 20, "input_unit": "条核酸 contigs", "output_count": 12, "output_unit": "条蛋白 ORFs", "duration_minutes": 1.25},
            {"stage_key": "abundance", "stage": "reads 回比丰度", "status": "done", "input_count": 5, "input_unit": "条核酸 contigs", "output_count": 5, "output_unit": "条核酸 contigs", "duration_minutes": 0.5},
            {"stage_key": "report", "stage": "结果汇总与可视化", "status": "done", "input_count": 1, "input_unit": "份报告", "output_count": 1, "output_unit": "份报告", "duration_minutes": 0.25},
            {"stage_key": "phylo", "stage": "系统发育放置", "status": "skipped", "input_count": 0, "input_unit": "条代表蛋白", "output_count": 0, "output_unit": "条代表蛋白", "duration_minutes": 0},
        ],
        "total_minutes": 1.5,
        "warnings": ["合成警告"],
        "warnings_en": ["Synthetic warning"],
        "links": {
            "main_results_report_zh.html": "main_results_report_zh.html",
            "main_results_report_en.html": "main_results_report_en.html",
            "phylogenetic_explorer.html": "phylogenetic_explorer.html",
            "representative_map.tsv": "tables/representative_map.tsv",
            "representative_orfs.faa": "sequences/representative_orfs.faa",
            "rdrp_positive_orfs.faa": "sequences/rdrp_positive_orfs.faa",
            "risk_warning.tsv": "tables/risk_warning.tsv",
            "neighbor_evidence.tsv": "tables/neighbor_evidence.tsv",
            "run_summary.json": "tables/run_summary.json",
            "delivery_manifest.json": "tables/delivery_manifest.json",
        },
    }
    context = SimpleNamespace(config={"project": {"display_name": "lpvid"}}, route=route)
    return context, rows, summary


def render_fixture(tmp_path, lang="en", *, route="reads", empty=False):
    context, rows, summary = rich_fixture(route)
    if empty:
        rows = []
        summary.update({
            "counts": {"positive_orfs": 0, "representatives": 0, "classified": 0, "placed": 0, "pending_tree": 0, "unassigned": 0, "positive_contigs": 0},
            "risk_rows": [], "risk_summary": {"total": 0, "human_nn1": 0, "human_nn2_5": 0, "vertebrate_only": 0, "no_placed_label": 0},
            "supergroup_stats": [], "multi_supergroup_contigs": [], "multi_supergroup_contig_count": 0,
            "abundance_rows": [],
            "stage_metrics": [], "warnings": [], "warnings_en": [], "links": {},
        })
        for row in summary["funnel"]:
            row["count"] = 0
    text = render_html(context, rows, summary, lang=lang)
    path = tmp_path / f"main_results_report_{lang}.html"
    path.write_text(text, encoding="utf-8")
    return path, text


def render_negative_fixture(tmp_path, lang="en"):
    """Zero-positive run that still has non-empty stage metrics (real negative shape)."""
    context, rows, summary = rich_fixture("reads")
    rows = []
    summary.update({
        "counts": {"positive_orfs": 0, "representatives": 0, "classified": 0, "placed": 0, "pending_tree": 0, "unassigned": 0, "positive_contigs": 0},
        "risk_rows": [], "risk_summary": {"total": 0, "human_nn1": 0, "human_nn2_5": 0, "vertebrate_only": 0, "no_placed_label": 0},
        "supergroup_stats": [], "multi_supergroup_contigs": [], "multi_supergroup_contig_count": 0,
        "abundance_rows": [],
        "stage_metrics": [
            {"stage_key": "representatives", "stage": "阳性序列去冗余", "status": "not_run", "input_count": 0, "input_unit": "条蛋白 ORFs", "output_count": 0, "output_unit": "条代表蛋白", "duration_minutes": 0},
            {"stage_key": "phylo", "stage": "系统发育放置", "status": "not_run", "input_count": 0, "input_unit": "条代表蛋白", "output_count": 0, "output_unit": "条代表蛋白", "duration_minutes": 0},
        ],
        "warnings": [], "warnings_en": [], "links": {},
    })
    for row in summary["funnel"]:
        row["count"] = 0
    text = render_html(context, rows, summary, lang=lang)
    path = tmp_path / f"negative_report_{lang}.html"
    path.write_text(text, encoding="utf-8")
    return path, text


def status_stat(supergroup, statuses, representatives=1, *, abundance=None):
    return {
        "supergroup": supergroup,
        "representatives": representatives,
        "positive_orfs": representatives,
        "positive_contigs": representatives,
        "length_statistics": length_statistics([{"aa_len": 200}] * representatives),
        "classification_methods": {"reference": representatives},
        "tree_statuses": statuses,
        "abundance": abundance if abundance is not None else {
            "n_contigs": 1, "mapped_reads": 5, "RPM": 2.5,
            "mean_depth_weighted": 1.0, "coverage_breadth_pct": 50.0,
        },
    }


def status_fixture(route="reads"):
    context, rows, summary = rich_fixture(route)
    zero_abundance = {
        "n_contigs": 0, "mapped_reads": 0, "RPM": 0.0,
        "mean_depth_weighted": 0.0, "coverage_breadth_pct": 0.0,
    }
    summary["supergroup_stats"] = [
        status_stat("SG1", {"assigned": 2}, representatives=2),
        status_stat("SG_EMB_AVAILABLE", {"classified_available": 1}),
        status_stat("SG_NO_TREE", {"classified_no_tree": 1}),
        status_stat("SG_PENDING_MERGE", {"classified_pending_merge": 1}),
        status_stat("SG_PENDING_RENAME", {"classified_pending_rename": 1}),
        status_stat("SG_BARE_PENDING", {"pending_merge": 1}),
        status_stat("SG_BARE_RENAME", {"pending_rename": 1}),
        status_stat("SG_MIXED", {"assigned": 1, "classified_no_tree": 1}, representatives=2),
        status_stat("SG_EMPTY_STATUS", {}),
        status_stat("SG_UNKNOWN", {"unknown_status": 1}),
        status_stat("SG_MISSING_COVERAGE", {"assigned": 1}, abundance={
            "n_contigs": 1, "mapped_reads": 2, "RPM": 3.0,
            "mean_depth_weighted": 1.5, "coverage_breadth_pct": None,
        }),
        status_stat("Unassigned", {"unassigned": 1}, abundance=zero_abundance),
    ]
    return context, rows, summary


def render_status_fixture(lang="en", *, route="reads"):
    context, rows, summary = status_fixture(route)
    return render_html(context, rows, summary, lang=lang)


def supergroup_card(text, heading):
    body = body_without_code(text)
    for card in re.findall(r'<article class="supergroup-card">.*?</article>', body, re.S):
        title = re.search(r'<h3>(.*?)</h3>', card, re.S)
        assert title
        plain = html.unescape(re.sub(r"<[^>]+>", "", title.group(1))).strip()
        if plain == heading:
            return card
    raise AssertionError(f"未找到超群卡片: {heading}")


def div_span(text, marker):
    """Return the balanced <div> span that starts at marker."""
    start = text.index(marker)
    depth = 0
    for match in re.finditer(r"<div\b|</div>", text[start:]):
        if match.group(0) == "</div>":
            depth -= 1
            if depth == 0:
                return start, start + match.end()
        else:
            depth += 1
    raise AssertionError(f"未闭合的 div 容器: {marker}")


def tag_span(text, marker, tag):
    """Return the balanced element span that starts at marker."""
    start = text.index(marker)
    depth = 0
    for match in re.finditer(rf"<{tag}\b|</{tag}>", text[start:]):
        if match.group(0) == f"</{tag}>":
            depth -= 1
            if depth == 0:
                return start, start + match.end()
        else:
            depth += 1
    raise AssertionError(f"未闭合的 {tag} 容器: {marker}")


def opening_tag(text, marker):
    start = text.index(marker)
    match = re.match(r"<[a-z]+[^>]*>", text[start:])
    assert match
    return match.group(0)


def family_card(text, family):
    marker = f'<details class="risk-family-card" data-family="{family}"'
    start, end = tag_span(text, marker, "details")
    return text[start:end]


def visible_text(text):
    """Rendered text only: drop scripts/styles/comments/tags, keep entity-free prose."""
    body = body_without_code(text)
    body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
    return re.sub(r"<[^>]+>", " ", body)


def test_exact_bilingual_section_titles_and_dom_navigation_order(tmp_path):
    for lang, titles in (("zh", ZH_TITLES), ("en", EN_TITLES)):
        _, text = render_fixture(tmp_path, lang)
        assert re.findall(r'<section id="([^"]+)"', text) == SECTION_IDS
        nav = re.search(r"<nav[^>]*>(.*?)</nav>", text, re.S).group(1)
        assert re.findall(r'href="#([^"]+)"', nav) == SECTION_IDS
        assert section_titles(text) == titles
        for retired in ("管理层摘要", "结果解读", "chartCalls", "chartFunnel", "riskTable", "supergroupTable"):
            assert retired not in text
        assert "echarts" not in text.lower()


def test_risk_uses_summary_real_ranks_nulls_and_prioritized_family_cards(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    data = report_data(text)
    assert data["riskSummary"]["human_nn2_5"] == 1
    for key in ("human-nn1", "human-nn2-5", "vertebrate-only", "no-placed-label"):
        assert f'id="risk-{key}">1</strong>' in text
    banner = re.search(r'<div id="risk-banner".*?</div>', text, re.S).group(0)
    assert "1 query / 1 family matched at least one list" in banner
    assert "4 placed queries = 1 matched + 2 unmatched + 1 with undetermined family classification" in banner
    assert "3 identified non-empty families" in banner
    assert text.index('data-family="ZetaMatchedviridae"') < text.index('data-family="AlphaUnmatchedviridae"')
    assert "orf-Q1" in text and "SG1" in text
    assert "Neighbor 4" in text and "human_ref_4" in text
    assert "Distance unavailable" in text
    assert "Human label unavailable" in text
    assert "No placed human/vertebrate label" in text
    assert "does not mean host absence or no risk" in text
    assert "not a risk classification or pathogen confirmation" in text
    assert "does not exclude unlisted or unidentified viruses" in text
    assert 'id="riskTable"' not in text


def test_english_risk_is_structural_and_never_reuses_chinese_examples_or_hints(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    assert "Category II" in text
    assert "Class B" in text
    for chinese in ("第二类", "乙类", "示例冠状病毒", "示例法定传染病", "中文提示不得复用到英文页面"):
        assert chinese not in text
    assert not re.search(r"[一-鿿]", body_without_code(text))
    _, zh_text = render_fixture(tmp_path, "zh")
    assert "第二类" in zh_text and "乙类" in zh_text and "示例冠状病毒" in zh_text


def test_representative_flow_reference_column_and_or_search(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    data = report_data(text)
    assert 'id="flow-prefilter">0</strong>' in text
    assert len(re.findall(r'class="flow-node"', text)) == 5
    ref_column = next(column for column in data["sequenceColumns"] if column["field"] == "reference_id")
    assert ref_column["title"] == "Original reference-sequence ID"
    assert any(row["reference_id"] == "ORIGINAL_REF_1" for row in data["sequences"])
    assert '<option value="">All supergroups</option>' in text
    assert "Object.values(row).some" in text
    assert "field:'query_id',type:'like'" not in text
    _, zh_text = render_fixture(tmp_path, "zh")
    assert "原始参考序列 ID" in zh_text
    assert '<option value="">全部超群</option>' in zh_text


def test_statistics_use_reconciled_stats_dynamic_height_and_data_twins(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    data = report_data(text)
    sg1 = next(row for row in data["supergroupRows"] if row["supergroup"] == "SG1")
    assert (sg1["representatives"], sg1["positive_orfs"], sg1["positive_contigs"]) == (2, 4, 2)
    assert 'id="supergroupChartFrame" style="height:248px"' in text
    assert 'id="lengthBinData"' in text and "<summary" in text
    assert re.search(r'<details id="lengthBinData".*?<table', text, re.S)
    assert 'id="stageTable"' in text and 'id="chartTime"' in text
    assert "Backbone tree under construction" in text
    assert "counted in more than one supergroup" in text
    assert "c2" in text and "SG1" in text and "SG2" in text
    assert "20.5" in text
    app_js = text.rsplit("<script>", 1)[1].split("</script>", 1)[0]
    assert "yAxes" not in app_js and "y2" not in app_js
    assert "\n      scales:horizontal" in app_js
    assert "\n    scales:horizontal" not in app_js


def test_sequence_table_drops_legacy_probability_columns(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    data = report_data(text)
    fields = {column["field"] for column in data["sequenceColumns"]}
    for legacy in ("host_human", "host_mammal", "host_vertebrate", "host_suppression_reason"):
        assert legacy not in fields
        assert legacy not in text
    seq = next(row for row in data["sequences"] if row["query_id"] == "Q1")
    assert seq["rpm"] == 12.5
    assert "suppression" not in text.lower()


def test_abundance_section_tables_and_route_notes(tmp_path):
    _, en = render_fixture(tmp_path, "en")
    assert 'id="abundanceTable"' in en
    data = report_data(en)
    assert data["abundanceRows"] == [
        {"contig_id": "c1", "length": 1000, "mapped_reads": 20, "rpm": 12.5, "mean_depth": 3.0, "coverage_breadth_pct": 80.0},
        {"contig_id": "c2", "length": 1200, "mapped_reads": 10, "rpm": 8.0, "mean_depth": 2.0, "coverage_breadth_pct": 70.0},
    ]
    assert "Relative abundance by contig" in en
    assert "c1" in en and "c2" in en
    for legacy in ("hostTable", "host-coverage", "pbv-suppression"):
        assert legacy not in en
    _, zh = render_fixture(tmp_path, "zh")
    assert "相对丰度（按 contig）" in zh

    _, proteins = render_fixture(tmp_path, "en", route="proteins")
    assert "Abundance mapping is not applicable to the proteins input route" in proteins
    assert 'id="abundanceTable"' in proteins
    _, zh_proteins = render_fixture(tmp_path, "zh", route="proteins")
    assert "路线不适用" in zh_proteins


def test_workflow_links_tree_preview_and_forbidden_product_names(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    workflow = re.search(r'<section id="workflow".*?</section>', text, re.S).group(0)
    assert 'id="stageTable"' in workflow and 'id="chartTime"' in workflow
    assert "chartFunnel" not in workflow
    assert workflow.count('id="stageTable"') == 1
    for forbidden in FORBIDDEN_PRODUCT_NAMES:
        assert forbidden not in text
    links = re.search(r'<section id="links".*?</section>', text, re.S).group(0)
    for anchor in re.findall(r"<a\b[^>]*>", links):
        assert 'target="_blank"' in anchor and 'rel="noopener"' in anchor
    assert "Structured top-5 neighbor evidence" in links
    assert "Prior-based family notices" in links
    start, end = div_span(text, '<div class="tree-preview no-print">')
    preview = text[start:end]
    assert preview.endswith("</div>")
    assert '<iframe id="treePreview" loading="lazy" src="phylogenetic_explorer.html"' in preview
    standalone = re.search(r'<a id="treeStandalone"[^>]*>', preview).group(0)
    assert 'target="_blank"' in standalone and 'rel="noopener"' in standalone
    assert "per-query" not in text.lower()


def test_dynamic_tables_empty_results_and_route_specific_abundance(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    assert "function makeTableOptions" in text
    assert "rowCount > 20" in text and "rowCount > 10" in text
    assert "height:'300px'" not in text and 'height:"300px"' not in text
    _, empty = render_fixture(tmp_path, "en", empty=True)
    assert 'id="treePreview"' not in empty and 'id="treeStandalone"' not in empty
    assert empty.count('class="flow-node"') == 5
    assert 'id="flow-input">0</strong>' in empty
    assert "No representative sequences are available" in empty
    assert "No sequences were placed" in empty
    assert "No supergroup statistics are available" in empty
    assert "No stage metrics are available" in empty
    assert "No abundance rows are available for this reads run" in empty
    _, proteins = render_fixture(tmp_path, "en", route="proteins")
    assert "Abundance mapping is not applicable to the proteins input route" in proteins


def test_supergroup_tree_status_business_mapping_and_raw_audit():
    zh = render_status_fixture("zh")
    en = render_status_fixture("en")
    zh_expect = {
        "SG1": "有可用骨架树",
        "SG_EMB_AVAILABLE": "有可用骨架树",
        "SG_NO_TREE": "骨架树正在建设中",
        "SG_PENDING_MERGE": "骨架树正在建设中",
        "SG_PENDING_RENAME": "骨架树正在建设中",
        "SG_BARE_PENDING": "骨架树正在建设中",
        "SG_BARE_RENAME": "骨架树正在建设中",
        "SG_MIXED": "部分有可用骨架树，其余正在建设中",
        "SG_EMPTY_STATUS": "不可用",
        "SG_UNKNOWN": "不可用",
        "未分类": "未分类，不进入系统发育放置",
    }
    en_expect = {
        "SG1": "Backbone tree available",
        "SG_EMB_AVAILABLE": "Backbone tree available",
        "SG_NO_TREE": "Backbone tree under construction",
        "SG_PENDING_MERGE": "Backbone tree under construction",
        "SG_PENDING_RENAME": "Backbone tree under construction",
        "SG_BARE_PENDING": "Backbone tree under construction",
        "SG_BARE_RENAME": "Backbone tree under construction",
        "SG_MIXED": "Mixed availability",
        "SG_EMPTY_STATUS": "Unavailable",
        "SG_UNKNOWN": "Unavailable",
        "Unclassified": "Unclassified; not placed",
    }
    for text, expected in ((zh, zh_expect), (en, en_expect)):
        for heading, status_text in expected.items():
            card = supergroup_card(text, heading)
            assert status_text in card
        available = supergroup_card(text, "SG1")
        assert "建设中" not in available and "under construction" not in available.lower()
        unassigned = supergroup_card(text, "未分类" if text is zh else "Unclassified")
        assert "建设中" not in unassigned and "under construction" not in unassigned.lower()
        assert "有可用骨架树" not in unassigned and "Backbone tree available" not in unassigned
        for heading in ("SG_EMPTY_STATUS", "SG_UNKNOWN"):
            empty = supergroup_card(text, heading)
            assert "建设中" not in empty and "under construction" not in empty.lower()
            assert "有可用骨架树" not in empty and "Backbone tree available" not in empty

    data = report_data(zh)
    raw_rows = {row["supergroup"]: row for row in data["supergroupRows"]}
    assert raw_rows["SG1"]["tree_status"] == "有可用骨架树"
    assert raw_rows["SG_MIXED"]["tree_status"] == "部分有可用骨架树，其余正在建设中"
    assert raw_rows["SG_EMPTY_STATUS"]["tree_status"] == "不可用"
    assert raw_rows["SG_UNKNOWN"]["tree_status"] == "不可用"
    assert raw_rows["Unassigned"]["tree_status"] == "未分类，不进入系统发育放置"
    assert raw_rows["SG_MIXED"]["tree_statuses"] == {"assigned": 1, "classified_no_tree": 1}
    assert raw_rows["SG_EMPTY_STATUS"]["tree_statuses"] == {}
    assert raw_rows["SG_UNKNOWN"]["tree_statuses"] == {"unknown_status": 1}
    assert raw_rows["Unassigned"]["tree_statuses"] == {"unassigned": 1}
    for row in raw_rows.values():
        assert "host" not in row


def test_abundance_zero_contigs_is_unavailable_without_zero_values():
    zh = render_status_fixture("zh")
    zh_unassigned = supergroup_card(zh, "未分类")
    assert "不可用（reads 路线无可用数据）" in zh_unassigned
    assert "RPM" not in zh_unassigned
    for pattern in (r"RPM\s*0", r"reads\s*0", r"深度\s*0", r"覆盖度\s*0"):
        assert not re.search(pattern, zh_unassigned)

    en = render_status_fixture("en")
    en_unassigned = supergroup_card(en, "Unclassified")
    assert "Unavailable (no usable data on the reads route)" in en_unassigned
    assert "RPM" not in en_unassigned
    for pattern in (r"RPM\s*0", r"reads\s*0", r"depth\s*0", r"coverage\s*0"):
        assert not re.search(pattern, en_unassigned, re.I)

    zh_cov = supergroup_card(zh, "SG_MISSING_COVERAGE")
    assert "覆盖度 不可用" in zh_cov and "不可用%" not in zh_cov
    en_cov = supergroup_card(en, "SG_MISSING_COVERAGE")
    assert "coverage unavailable" in en_cov and "unavailable%" not in en_cov

    zh_proteins = render_status_fixture("zh", route="proteins")
    zh_protein_card = supergroup_card(zh_proteins, "SG1")
    assert "路线不适用（非 reads 输入）" in zh_protein_card
    assert "RPM" not in zh_protein_card
    en_proteins = render_status_fixture("en", route="proteins")
    en_protein_card = supergroup_card(en_proteins, "SG1")
    assert "Not applicable for this input route" in en_protein_card
    assert "RPM" not in en_protein_card


def test_report_uses_natural_data_labels_and_no_internal_twin_jargon(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    _, en = render_fixture(tmp_path, "en")
    assert "twin" not in zh.lower()
    assert "twin" not in en.lower()
    zh_body = body_without_code(zh)
    en_body = body_without_code(en)
    assert "图表对应数据" in zh_body and "数据明细" in zh_body
    assert "长度分布数据表" in zh_body
    assert "长度 bin 数据表" not in zh_body
    assert "图表数据 twin" not in zh_body and "耗时图数据 twin" not in zh_body
    assert "chart data" in en_body.lower() and "data details" in en_body.lower()
    assert "Length distribution data table" in en_body
    assert "Length-bin data table" not in en_body


def test_representative_table_localizes_codes_methods_and_reference_note(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    zh_data = report_data(zh)
    zh_seq = {row["query_id"]: row for row in zh_data["sequences"]}
    assert zh_seq["Q1"]["classification_status"] == "有可用骨架树"
    assert zh_seq["Q1"]["classification_method"] == "参考序列匹配"
    assert zh_seq["Q3"]["classification_status"] == "骨架树正在建设中"
    assert zh_seq["Q3"]["classification_method"] == "分类模型"
    assert zh_seq["Q4"]["classification_status"] == "骨架树正在建设中"
    assert zh_seq["Q5"]["supergroup"] == "未分类"
    assert zh_seq["Q5"]["classification_status"] == "未分类，不进入系统发育放置"
    zh_ref = next(column for column in zh_data["sequenceColumns"] if column["field"] == "reference_id")
    assert "不等于当前超群归属" in zh_ref["headerTooltip"]
    assert "不等于当前超群归属" in body_without_code(zh)
    for code in ("assigned", "classified_no_tree", "classified_pending_merge", "unassigned", "Unassigned"):
        assert code not in body_without_code(zh)
        assert all(row["classification_status"] != code for row in zh_data["sequences"])
        assert all(row["supergroup"] != code for row in zh_data["sequences"])

    _, en = render_fixture(tmp_path, "en")
    en_data = report_data(en)
    en_seq = {row["query_id"]: row for row in en_data["sequences"]}
    assert en_seq["Q1"]["classification_status"] == "Backbone tree available"
    assert en_seq["Q1"]["classification_method"] == "Reference match"
    assert en_seq["Q3"]["classification_status"] == "Backbone tree under construction"
    assert en_seq["Q3"]["classification_method"] == "Classification model"
    assert en_seq["Q5"]["supergroup"] == "Unclassified"
    assert en_seq["Q5"]["classification_status"] == "Unclassified; not placed"
    en_ref = next(column for column in en_data["sequenceColumns"] if column["field"] == "reference_id")
    assert "not the current supergroup assignment" in en_ref["headerTooltip"]
    assert "not the current supergroup assignment" in body_without_code(en)
    for code in ("assigned", "classified_no_tree", "classified_pending_merge", "unassigned", "Unassigned"):
        assert code not in body_without_code(en)


def test_link_descriptions_include_orf_files_and_manifest_status_only(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    _, en = render_fixture(tmp_path, "en")
    zh_links = re.search(r'<section id="links".*?</section>', body_without_code(zh), re.S).group(0)
    assert "代表 ORF 蛋白序列" in zh_links
    assert "RdRP 阳性 ORF 蛋白序列" in zh_links
    assert "交付完成状态与文件清单" in zh_links
    assert "校验和" not in zh_links and "校验信息" not in zh_links
    en_links = re.search(r'<section id="links".*?</section>', body_without_code(en), re.S).group(0)
    assert "Representative ORF protein sequences" in en_links
    assert "RdRP-positive ORF protein sequences" in en_links
    assert "Delivery completion status and file list" in en_links
    assert "checksum" not in en_links.lower()


def test_workflow_stage_statuses_and_minutes_label_are_localized(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    zh_data = report_data(zh)
    assert [row["status"] for row in zh_data["stageRows"]] == ["完成", "完成", "完成", "跳过"]
    assert zh_data["labels"]["minutes"] == "分钟"
    _, en = render_fixture(tmp_path, "en")
    en_data = report_data(en)
    assert [row["status"] for row in en_data["stageRows"]] == ["Completed", "Completed", "Completed", "Skipped"]
    assert en_data["labels"]["minutes"] == "min"
    app_js = en.rsplit("<script>", 1)[1].split("</script>", 1)[0]
    assert "${Number(context.parsed.x).toLocaleString()} ${L.minutes}" in app_js
    assert "${fmtUnavailable(cell.getValue())} ${L.minutes}" in app_js


def test_risk_known_family_without_list_match_is_not_unavailable(tmp_path):
    _, en = render_fixture(tmp_path, "en")
    en_body = body_without_code(en)
    q2 = re.search(r'<div class="risk-query"><h4>Query Q2</h4>.*?</div>', en_body, re.S).group(0)
    assert q2.count("No match in supplied lists") == 2
    assert "unavailable" not in q2
    assert "does not exclude unlisted" in q2
    q4 = re.search(r'<div class="risk-query"><h4>Query Q4</h4>.*?</div>', en_body, re.S).group(0)
    assert q4.count("taxonomy unavailable") == 2
    assert "No match in supplied lists" not in q4

    _, zh = render_fixture(tmp_path, "zh")
    zh_body = body_without_code(zh)
    q2_zh = re.search(r'<div class="risk-query"><h4>代表序列 Q2</h4>.*?</div>', zh_body, re.S).group(0)
    assert q2_zh.count("未匹配名录") == 2
    assert "不可用" not in q2_zh
    assert "不排除名录外或未识别病毒" in zh_body
    q4_zh = re.search(r'<div class="risk-query"><h4>代表序列 Q4</h4>.*?</div>', zh_body, re.S).group(0)
    assert q4_zh.count("分类学信息不可用") == 2
    assert "未匹配名录" not in q4_zh


def test_missing_neighbor_rank_never_renders_nn_none(tmp_path):
    _, en = render_fixture(tmp_path, "en")
    en_body = body_without_code(en)
    assert "nnNone" not in en_body
    assert "rank unavailable" in en_body
    _, zh = render_fixture(tmp_path, "zh")
    zh_body = body_without_code(zh)
    assert "nnNone" not in zh_body
    assert "排序位次不可用" in zh_body


def test_english_warnings_never_fall_back_to_chinese(tmp_path):
    for warnings_en in (None, []):
        context, rows, summary = rich_fixture("reads")
        if warnings_en is None:
            summary.pop("warnings_en")
        else:
            summary["warnings_en"] = warnings_en
        text = render_html(context, rows, summary, lang="en")
        assert "合成警告" not in text
        assert "<h3>Warnings</h3>" not in text

    context, rows, summary = rich_fixture("reads")
    en = render_html(context, rows, summary, lang="en")
    assert "Synthetic warning" in en and "合成警告" not in en
    zh = render_html(context, rows, summary, lang="zh")
    assert "合成警告" in zh and "Synthetic warning" not in zh


def test_workflow_stage_statuses_cover_not_run_failed_running_stale():
    statuses = ["not_run", "failed", "running", "stale"]
    stage_keys = [("detect", "RdRP 识别模型"), ("prefilter", "宿主/细胞同源粗筛"),
                  ("abundance", "reads 回比丰度"), ("supergroup", "超群分类")]
    context, rows, summary = rich_fixture("reads")
    summary["stage_metrics"] = [
        {"stage_key": key, "stage": label, "status": status, "input_count": 0,
         "input_unit": "条蛋白 ORFs", "output_count": 0, "output_unit": "条代表蛋白",
         "duration_minutes": 0.0}
        for (key, label), status in zip(stage_keys, statuses)
    ]
    zh = render_html(context, rows, summary, lang="zh")
    assert [row["status"] for row in report_data(zh)["stageRows"]] == ["未运行", "失败", "运行中", "过时"]
    en = render_html(context, rows, summary, lang="en")
    assert [row["status"] for row in report_data(en)["stageRows"]] == ["Not run", "Failed", "Running", "Stale"]


def test_verifier_accepts_legitimate_empty_and_negative_bilingual_reports(tmp_path):
    for lang in ("zh", "en"):
        empty_path, _ = render_fixture(tmp_path, lang, empty=True)
        assert check_file(empty_path) == []
        negative_path, _ = render_negative_fixture(tmp_path, lang)
        assert check_file(negative_path) == []


def test_verifier_rejects_risk_kpi_count_mismatch(tmp_path):
    path, _ = render_fixture(tmp_path, "en")
    text = path.read_text(encoding="utf-8")
    text = text.replace('<strong id="risk-human-nn1">1</strong>', '<strong id="risk-human-nn1">9</strong>', 1)
    path.write_text(text, encoding="utf-8")
    assert any("风险 KPI" in problem for problem in check_file(path))


def test_verifier_rejects_iframe_outside_tree_preview_container(tmp_path):
    path, text = render_fixture(tmp_path, "en")
    start, end = div_span(text, '<div class="tree-preview no-print">')
    container = text[start:end]
    iframe = re.search(r'<iframe id="treePreview"[^>]*></iframe>', container).group(0)
    broken = text[:start] + container.replace(iframe, "") + iframe + text[end:]
    path.write_text(broken, encoding="utf-8")
    assert any("iframe" in problem and "tree-preview" in problem for problem in check_file(path))


def make_context(tmp_path):
    representative = write_tsv(
        tmp_path / "rep.tsv",
        ["query_id", "representative_id", "aa_len", "motif_class_count", "motif_classes", "is_representative"],
        [{"query_id": "Q1<script>", "representative_id": "Q1<script>", "aa_len": 200, "motif_class_count": 2, "motif_classes": "A;C", "is_representative": 1}],
    )
    assignments = write_tsv(
        tmp_path / "sg.tsv", ["query_id", "supergroup", "status"],
        [{"query_id": "Q1<script>", "supergroup": "SG1</script><script>alert(1)</script>", "status": "assigned"}],
    )
    return SimpleNamespace(
        data={"representative_map": str(representative), "supergroup_assignments": str(assignments)},
        config=load_config(), run_dir=tmp_path / "run", route="proteins", sample="<sample>",
    )


def test_html_is_selfcontained_and_escapes_dynamic_text(tmp_path):
    context = make_context(tmp_path)
    summary = {
        "run_id": "run<script>", "route": "proteins", "sample": "<sample>",
        "counts": {"positive_orfs": 1, "representatives": 1, "assigned": 1, "unassigned": 0, "placed": 0, "positive_contigs": 0},
        "warnings": ["test warning"], "links": {"tree_explorer.html": "tree_explorer.html"},
        "stages": {"report": {"status": "done"}},
    }
    output = build_html_report(context, tmp_path / "report.html", summary, lang="en")
    text = output.read_text()
    assert "lpvid" in text and "REPORT_DATA" in text
    assert "&lt;sample&gt;" in text and "<sample>" not in text
    assert "Q1\\u003cscript\\u003e" in text
    assert "SG1\\u003c/script\\u003e\\u003cscript\\u003ealert(1)\\u003c/script\\u003e" in text
    assert "</script><script>alert(1)</script>" not in text
    assert 'id="globalFilter"' in text and 'id="supergroupFilter"' in text
    assert "Tabulator" in text and "echarts" not in text.lower()
    assert "<script src=" not in text and '<link rel="stylesheet" href=' not in text


def test_verifier_accepts_only_explicit_synthetic_paths_and_catches_retired_markup(tmp_path, capsys):
    zh_path, _ = render_fixture(tmp_path, "zh")
    en_path, en_text = render_fixture(tmp_path, "en")
    assert check_file(zh_path) == []
    assert check_file(en_path) == []
    assert verify_main(["verify_report_v2.py"]) == 2
    assert "explicit" in capsys.readouterr().err
    assert verify_main(["verify_report_v2.py", str(zh_path), str(en_path)]) == 0
    assert "2/2" in capsys.readouterr().out
    en_path.write_text(en_text.replace('id="chartTime"', 'id="chartCalls"'), encoding="utf-8")
    problems = check_file(en_path)
    assert any("chartTime" in problem for problem in problems)
    assert any("chartCalls" in problem for problem in problems)


def _risk_evidence(query, family, *, category="", statutory="", human=False, vertebrate=False):
    evidence = neighbor(1, f"nbr-{query}", distance=0.1, family=family, human=human, vertebrate=vertebrate)
    return risk_row(
        query, f"contig-{query}", "SG1", family, [evidence],
        category=category, statutory=statutory, distance=0.1,
    )


def large_risk_fixture(route="reads"):
    context, rows, summary = rich_fixture(route)
    risks = (
        [_risk_evidence(f"ML{i}", "MatchedLargeviridae", category="第二类") for i in range(1, 5)]
        + [_risk_evidence(f"MS{i}", "MatchedSmallviridae", statutory="乙类", human=True) for i in range(1, 3)]
        + [_risk_evidence(f"UA{i}", "AlphaUnmatchedviridae", vertebrate=True) for i in range(1, 3)]
        + [_risk_evidence(f"UB{i}", "BetaUnmatchedviridae") for i in range(1, 4)]
    )
    summary["risk_rows"] = risks
    summary["risk_summary"] = summarize_risk(risks)
    summary["counts"] = {**summary["counts"], "placed": len(risks)}
    return context, rows, summary


def small_coronavirus_fixture():
    context, rows, summary = rich_fixture("reads")
    risks = [_risk_evidence("FBZJ02", "Coronaviridae", category="第二类", human=True, vertebrate=True)]
    summary["risk_rows"] = risks
    summary["risk_summary"] = summarize_risk(risks)
    summary["counts"] = {**summary["counts"], "placed": 1}
    return context, rows, summary


def test_risk_family_cards_stay_collapsed_with_full_evidence_inside(tmp_path):
    context, rows, summary = large_risk_fixture()
    zh = render_html(context, rows, summary, lang="zh")
    large = family_card(zh, "MatchedLargeviridae")
    assert " open" not in opening_tag(large, '<details class="risk-family-card"')
    assert "4 条代表序列" in large
    for index in range(1, 5):
        assert f"代表序列 ML{index}" in large
        assert f"nbr-ML{index}" in large
    assert large.count('<div class="risk-query">') == 4

    small = family_card(zh, "MatchedSmallviridae")
    assert " open" not in opening_tag(small, '<details class="risk-family-card"')
    assert "2 条代表序列" in small
    for index in range(1, 3):
        assert f"代表序列 MS{index}" in small
        assert f"nbr-MS{index}" in small

    en = render_html(context, rows, summary, lang="en")
    large_en = family_card(en, "MatchedLargeviridae")
    assert " open" not in opening_tag(large_en, '<details class="risk-family-card"')
    assert "4 representative sequences" in large_en
    for index in range(1, 5):
        assert f"Query ML{index}" in large_en
        assert f"nbr-ML{index}" in large_en
    small_en = family_card(en, "MatchedSmallviridae")
    assert " open" not in opening_tag(small_en, '<details class="risk-family-card"')
    assert "2 representative sequences" in small_en


def test_risk_unmatched_families_are_grouped_closed_with_counts(tmp_path):
    context, rows, summary = large_risk_fixture()
    zh = render_html(context, rows, summary, lang="zh")
    marker = '<details class="risk-unmatched"'
    start, end = tag_span(zh, marker, "details")
    unmatched = zh[start:end]
    assert " open" not in opening_tag(zh, marker)
    assert "未匹配名录：2 个已识别科 · 5 条代表序列；科级分类未确定：0 条代表序列" in unmatched
    for family in ("AlphaUnmatchedviridae", "BetaUnmatchedviridae"):
        assert f'data-family="{family}"' in unmatched
    for family in ("MatchedLargeviridae", "MatchedSmallviridae"):
        family_index = zh.index(f'data-family="{family}"')
        assert not start <= family_index < end
    assert '<div id="risk-banner"' in zh
    assert zh.count('class="kpi risk-kpi"') == 4
    assert zh.count('<div class="risk-query">') == 11
    assert zh.count('<ul class="neighbor-list">') == 11
    risk_section = re.search(r'<section id="risk".*?</section>', body_without_code(zh), re.S).group(0)
    assert risk_section.count('class="evidence-grid"') == 11
    assert risk_section.count('class="risk-family-overview"') == 4
    for query in [*(f"ML{i}" for i in range(1, 5)), *(f"MS{i}" for i in range(1, 3)),
                  *(f"UA{i}" for i in range(1, 3)), *(f"UB{i}" for i in range(1, 4))]:
        assert f"代表序列 {query}" in zh
        assert f"nbr-{query}" in zh
    assert "<table" not in risk_section
    assert "不排除名录外或未识别病毒" in risk_section

    en = render_html(context, rows, summary, lang="en")
    en_start, en_end = tag_span(en, marker, "details")
    en_unmatched = en[en_start:en_end]
    assert " open" not in opening_tag(en, marker)
    assert "No catalog/list match: 2 families · 5 queries; undetermined family classification: 0 queries" in en_unmatched
    assert "does not exclude unlisted or unidentified viruses" in en


def test_risk_all_unmatched_keeps_group_counts_and_non_exclusion_note(tmp_path):
    context, rows, summary = large_risk_fixture()
    summary["risk_rows"] = [row for row in summary["risk_rows"] if not (row.pathogen_category or row.statutory_category)]
    summary["risk_summary"] = summarize_risk(summary["risk_rows"])
    zh = render_html(context, rows, summary, lang="zh")
    start, end = tag_span(zh, '<details class="risk-unmatched"', "details")
    assert "未匹配名录：2 个已识别科 · 5 条代表序列；科级分类未确定：0 条代表序列" in zh[start:end]
    assert " open" not in opening_tag(zh, '<details class="risk-unmatched"')
    assert "不排除名录外或未识别病毒" in zh
    assert "未匹配名录" in zh[start:end]
    en = render_html(context, rows, summary, lang="en")
    en_start, en_end = tag_span(en, '<details class="risk-unmatched"', "details")
    assert "No catalog/list match: 2 families · 5 queries; undetermined family classification: 0 queries" in en[en_start:en_end]
    assert "does not exclude unlisted or unidentified viruses" in en


def test_risk_small_matched_family_keeps_evidence_while_collapsed(tmp_path):
    context, rows, summary = small_coronavirus_fixture()
    zh = render_html(context, rows, summary, lang="zh")
    card = family_card(zh, "Coronaviridae")
    assert " open" not in opening_tag(card, '<details class="risk-family-card"')
    assert "1 条代表序列" in card
    assert "代表序列 FBZJ02" in card and "nbr-FBZJ02" in card
    assert '<dl class="evidence-grid">' in card and '<ul class="neighbor-list">' in card
    assert '<div class="risk-family-overview">' in card
    assert "不排除名录外或未识别病毒" in body_without_code(zh)


def mixed_family_count_fixture():
    context, rows, summary = rich_fixture("reads")
    risks = (
        [_risk_evidence(f"UP{i}", "Picobirnaviridae", vertebrate=True) for i in range(1, 6)]
        + [_risk_evidence(f"UNK{i}", "") for i in range(1, 3)]
        + [_risk_evidence("MATCH1", "Coronaviridae", category="第二类", human=True, vertebrate=True)]
    )
    summary["risk_rows"] = risks
    summary["risk_summary"] = summarize_risk(risks)
    summary["counts"] = {**summary["counts"], "placed": len(risks)}
    return context, rows, summary


def test_family_counts_separate_identified_families_and_unknown_bucket(tmp_path):
    context, rows, summary = mixed_family_count_fixture()
    zh = render_html(context, rows, summary, lang="zh")
    banner = re.search(r'<div id="risk-banner".*?</div>', zh, re.S).group(0)
    assert "1 条代表序列 / 1 个科匹配" in banner
    assert "已放置代表序列 8 条 = 1 匹配 + 5 未匹配 + 2 科级分类未确定" in banner
    assert "已识别科 2 个" in banner
    assert "已识别科 3 个" not in banner
    start, end = tag_span(zh, '<details class="risk-unmatched"', "details")
    unmatched = zh[start:end]
    assert "未匹配名录：1 个已识别科 · 5 条代表序列；科级分类未确定：2 条代表序列" in unmatched
    assert 'data-family="Picobirnaviridae"' in unmatched
    pico_card = family_card(zh, "Picobirnaviridae")
    assert " open" not in opening_tag(pico_card, '<details class="risk-family-card"')
    assert 'data-family=""' in unmatched
    assert "第 1 近邻所属科不可用" in unmatched
    assert summary["risk_summary"]["total"] == 8

    en = render_html(context, rows, summary, lang="en")
    en_banner = re.search(r'<div id="risk-banner".*?</div>', en, re.S).group(0)
    assert "1 query / 1 family matched at least one list" in en_banner
    assert "8 placed queries = 1 matched + 5 unmatched + 2 with undetermined family classification" in en_banner
    assert "2 identified non-empty families" in en_banner
    assert "3 identified non-empty families" not in en_banner
    en_start, en_end = tag_span(en, '<details class="risk-unmatched"', "details")
    en_unmatched = en[en_start:en_end]
    assert "No catalog/list match: 1 family · 5 queries; undetermined family classification: 2 queries" in en_unmatched
    assert 'data-family="Picobirnaviridae"' in en_unmatched
    assert "Nearest-neighbor family unavailable" in en_unmatched


def small_unmatched_fixture():
    context, rows, summary = rich_fixture("reads")
    risks = [
        _risk_evidence("UA1", "AlphaUnmatchedviridae"),
        _risk_evidence("UA2", "AlphaUnmatchedviridae"),
        _risk_evidence("UB1", "BetaUnmatchedviridae"),
    ]
    summary["risk_rows"] = risks
    summary["risk_summary"] = summarize_risk(risks)
    summary["counts"] = {**summary["counts"], "placed": len(risks)}
    return context, rows, summary


def test_risk_small_unmatched_sample_stays_collapsed_by_default(tmp_path):
    context, rows, summary = small_unmatched_fixture()
    zh = render_html(context, rows, summary, lang="zh")
    assert " open" not in opening_tag(zh, '<details class="risk-unmatched"')
    for family in ("AlphaUnmatchedviridae", "BetaUnmatchedviridae"):
        card = family_card(zh, family)
        assert " open" not in opening_tag(card, '<details class="risk-family-card"')
        assert '<dl class="evidence-grid">' in card and '<ul class="neighbor-list">' in card
    for query in ("UA1", "UA2", "UB1"):
        assert f"代表序列 {query}" in zh and f"nbr-{query}" in zh

    en = render_html(context, rows, summary, lang="en")
    assert " open" not in opening_tag(en, '<details class="risk-unmatched"')
    for family in ("AlphaUnmatchedviridae", "BetaUnmatchedviridae"):
        card = family_card(en, family)
        assert " open" not in opening_tag(card, '<details class="risk-family-card"')
    assert "No catalog/list match: 2 families · 3 queries; undetermined family classification: 0 queries" in en


def _extract_js_function(text, name):
    app_js = text.rsplit("<script>", 1)[1].split("</script>", 1)[0]
    start = app_js.index(f"function {name}(")
    depth = 0
    for index in range(start, len(app_js)):
        char = app_js[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return app_js[start:index + 1]
    raise AssertionError(f"未找到完整函数: {name}")


def test_print_expands_all_details_and_restores_screen_state(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    app_js = text.rsplit("<script>", 1)[1].split("</script>", 1)[0]
    assert "function setupPrintDetails()" in app_js
    assert "addEventListener('beforeprint'" in app_js
    assert "addEventListener('afterprint'" in app_js
    assert "setupPrintDetails();" in app_js

    css = re.search(r"<style>(.*?)</style>", text, re.S).group(1)
    print_rule = css[css.index("@media print{"):]
    assert "display:block!important" in print_rule
    for selector in ("risk-family-card", "risk-unmatched", "workflow-flow"):
        assert selector in print_rule
    for retired_selector in ("host-group-details", "host-suppression-details"):
        assert retired_selector not in css

    node = shutil.which("node")
    if node is None:
        pytest.skip("node 不可用，无法执行打印折叠/恢复逻辑测试")
    source_path = tmp_path / "print_details.js"
    source_path.write_text(_extract_js_function(text, "setupPrintDetails"), encoding="utf-8")
    harness_path = tmp_path / "print_details_harness.js"
    harness_path.write_text(
        """
const fs = require('fs');
const source = fs.readFileSync(process.argv[2], 'utf8');
const details = [{open:false},{open:true},{open:false}];
const listeners = {};
global.window = {addEventListener(type, callback){listeners[type]=callback;}};
global.document = {querySelectorAll(){return details;}};
eval(source + '\\nsetupPrintDetails();');
if (typeof listeners.beforeprint !== 'function' || typeof listeners.afterprint !== 'function') {
  throw new Error('print listeners missing');
}
listeners.beforeprint();
if (!details.every(item => item.open)) throw new Error('beforeprint did not expand all details');
listeners.afterprint();
if (details.map(item => item.open).join(',') !== 'false,true,false') {
  throw new Error('afterprint did not restore state');
}
details[0].open = true;
listeners.beforeprint();
listeners.afterprint();
if (details.map(item => item.open).join(',') !== 'true,true,false') {
  throw new Error('second print cycle did not save current state');
}
console.log('PRINT_DETAILS_OK');
""",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [node, str(harness_path), str(source_path)],
        capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    assert "PRINT_DETAILS_OK" in completed.stdout


def test_closed_risk_details_css_hides_content(tmp_path):
    _, text = render_fixture(tmp_path, "en")
    css = re.search(r"<style>(.*?)</style>", text, re.S).group(1)
    for selector in (
        r"\.risk-family-card:not\(\[open\]\)",
        r"\.risk-unmatched:not\(\[open\]\)",
    ):
        rule = re.search(selector + r"[^{}]*\{([^}]*)\}", css)
        assert rule and "display:none" in rule.group(1), selector
    assert ".host-group" not in css and ".host-coverage" not in css


def test_public_report_labels_avoid_internal_algorithm_terms(tmp_path):
    import json as _json
    body_terms = ["luca" + "prot", "luca" + "virus", "embed" + "ding", "m" + "lp", "stage" + "-2", "top" + "-hit", "diamond"]
    for lang in ("zh", "en"):
        _, text = render_fixture(tmp_path, lang)
        body = body_without_code(text).lower()
        serialized = _json.dumps(report_data(text), ensure_ascii=False).lower()
        for term in body_terms:
            assert term not in body, (lang, term)
            assert term not in serialized, (lang, term)


def test_zh_visible_text_avoids_latin_abbreviations_and_reads_naturally(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    text = visible_text(zh)
    for banned in ("nn1", "nn2", "nn3", "nn4", "nn5", "query", "Query", "rank", "top-5", "PBV"):
        assert banned not in text, banned
    for required in ("第 1 近邻", "第 2 至第 5 近邻", "代表序列", "排序位次", "前 5 个近邻"):
        assert required in text, required
    assert "病原微生物目录" in text and "法定传染病名单" in text


def test_zh_neighbor_items_show_all_five_structured_fields(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    risk = re.search(r'<section id="risk".*?</section>', body_without_code(zh), re.S).group(0)
    items = re.findall(r'<li class="neighbor-item">.*?</li>', risk, re.S)
    assert len(items) >= 5
    first = items[0]
    for label in ("第 1 近邻", "参考序列", "距离", "人源标签", "脊椎动物标签", "科"):
        assert label in first, label


def test_english_risk_uses_professional_neighbor_wording(tmp_path):
    _, en = render_fixture(tmp_path, "en")
    text = visible_text(en)
    assert not re.search(r"\bnn[1-5]\b", text)
    for phrase in ("Nearest neighbor", "Human label", "Vertebrate label", "Family"):
        assert phrase in text, phrase


def test_risk_family_cards_are_collapsed_details_with_badged_summaries(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    matched = family_card(zh, "ZetaMatchedviridae")
    assert " open" not in opening_tag(matched, '<details class="risk-family-card"')
    assert 'data-matched="true"' in matched
    summary = re.search(r"<summary>(.*?)</summary>", matched, re.S).group(1)
    assert "ZetaMatchedviridae" in summary and "1 条代表序列" in summary
    assert '<span class="badge badge-sev-mid"><strong>病原微生物目录</strong>第二类</span>' in summary
    assert '<span class="badge badge-sev-mid"><strong>法定传染病名单</strong>乙类</span>' in summary

    unmatched = family_card(zh, "AlphaUnmatchedviridae")
    assert 'data-matched="false"' in unmatched
    unmatched_summary = re.search(r"<summary>(.*?)</summary>", unmatched, re.S).group(1)
    assert '<span class="badge badge-sev-none">未匹配名录</span>' in unmatched_summary

    _, en = render_fixture(tmp_path, "en")
    en_matched = family_card(en, "ZetaMatchedviridae")
    en_summary = re.search(r"<summary>(.*?)</summary>", en_matched, re.S).group(1)
    assert "1 representative sequence" in en_summary
    assert '<span class="badge badge-sev-mid"><strong>Pathogen catalog</strong>Category II</span>' in en_summary
    assert '<span class="badge badge-sev-mid"><strong>Statutory disease list</strong>Class B</span>' in en_summary


def test_risk_badges_color_categories_by_severity(tmp_path):
    context, rows, summary = rich_fixture("reads")
    risks = [
        _risk_evidence("C1", "Highviridae", category="第一类", statutory="甲类"),
        _risk_evidence("C2", "Midviridae", category="第二类", statutory="乙类"),
        _risk_evidence("C3", "Lowviridae", category="第三类", statutory="丙类"),
        _risk_evidence("C4", "Cat4viridae", category="第四类"),
    ]
    summary["risk_rows"] = risks
    summary["risk_summary"] = summarize_risk(risks)
    summary["counts"] = {**summary["counts"], "placed": len(risks)}
    zh = render_html(context, rows, summary, lang="zh")
    high = re.search(r"<summary>(.*?)</summary>", family_card(zh, "Highviridae"), re.S).group(1)
    assert '<span class="badge badge-sev-high"><strong>病原微生物目录</strong>第一类</span>' in high
    assert '<span class="badge badge-sev-high"><strong>法定传染病名单</strong>甲类</span>' in high
    mid = re.search(r"<summary>(.*?)</summary>", family_card(zh, "Midviridae"), re.S).group(1)
    assert '<span class="badge badge-sev-mid"><strong>病原微生物目录</strong>第二类</span>' in mid
    assert '<span class="badge badge-sev-mid"><strong>法定传染病名单</strong>乙类</span>' in mid
    low = re.search(r"<summary>(.*?)</summary>", family_card(zh, "Lowviridae"), re.S).group(1)
    assert '<span class="badge badge-sev-neutral"><strong>病原微生物目录</strong>第三类</span>' in low
    assert '<span class="badge badge-sev-neutral"><strong>法定传染病名单</strong>丙类</span>' in low
    cat4 = re.search(r"<summary>(.*?)</summary>", family_card(zh, "Cat4viridae"), re.S).group(1)
    assert '<span class="badge badge-sev-neutral"><strong>病原微生物目录</strong>第四类</span>' in cat4
    css = re.search(r"<style>(.*?)</style>", zh, re.S).group(1)
    for severity in ("badge-sev-high", "badge-sev-mid", "badge-sev-neutral", "badge-sev-none"):
        assert f".{severity}{{" in css or f".{severity} {{" in css, severity


def test_family_overview_precedes_representative_evidence(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    card = family_card(zh, "ZetaMatchedviridae")
    assert card.index('class="risk-family-overview"') < card.index('class="risk-family-evidence"')
    overview = re.search(r'<div class="risk-family-overview">(.*?)</div>', card, re.S).group(1)
    assert "代表序列数" in overview and "1 条" in overview
    assert "涉及超群" in overview and "SG1" in overview
    assert "病原微生物目录" in overview and "法定传染病名单" in overview


def test_supergroup_cards_use_header_badge_metric_chips_and_evidence_grid(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    card = supergroup_card(zh, "SG1")
    assert re.search(r'<article class="supergroup-card"><h3>SG1</h3><span class="badge badge-sev-ok">有可用骨架树</span>', card)
    assert card.count('class="metric-chip"') == 3
    for label in ("代表序列", "全部阳性 ORFs", "唯一阳性 contigs"):
        assert label in card, label
    assert 'class="evidence-grid"' in card
    for dt_label in ("长度", "分类方法", "丰度"):
        assert f"<dt>{dt_label}" in card, dt_label
    assert "host" not in card.lower() and "宿主" not in card

    _, en = render_fixture(tmp_path, "en")
    en_card = supergroup_card(en, "SG3")
    assert re.search(r'<span class="badge badge-sev-mid">Backbone tree under construction</span>', en_card)
    unassigned = supergroup_card(zh, "未分类")
    assert re.search(r'<span class="badge badge-sev-neutral">未分类，不进入系统发育放置</span>', unassigned)


def test_workflow_flow_chart_chains_stages_with_io_and_total_time(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    workflow = re.search(r'<section id="workflow".*?</section>', zh, re.S).group(0)
    assert 'class="workflow-flow"' in workflow and 'id="workflowFlow"' in workflow
    assert 'id="workflowTotal"' in workflow and "总耗时" in workflow
    assert "1.5 分钟" in workflow
    assert workflow.count('class="workflow-flow-node"') == 4
    orf_node = re.search(r'<div class="workflow-flow-node" data-stage-key="orf">.*?</div>', workflow, re.S).group(0)
    assert "蛋白 ORF 预测" in orf_node
    assert "输入 20 条核酸 contigs" in orf_node
    assert "输出 12 条蛋白 ORFs" in orf_node
    assert "耗时 1.25 分钟" in orf_node
    assert workflow.count('class="workflow-flow-arrow"') == 3
    assert 'id="stageTable"' in workflow and 'id="chartTime"' in workflow

    _, en = render_fixture(tmp_path, "en")
    en_workflow = re.search(r'<section id="workflow".*?</section>', en, re.S).group(0)
    assert "Total time" in en_workflow and "1.5 min" in en_workflow
    assert "Input 20 nucleotide contigs" in en_workflow
    assert "Output 12 protein ORFs" in en_workflow
    assert "Time 1.25 min" in en_workflow

    css = re.search(r"<style>(.*?)</style>", zh, re.S).group(1)
    assert ".workflow-flow" in css and ".workflow-flow-node" in css
    responsive = css[css.index("@media(max-width:700px"):]
    assert "workflow-flow-arrow" in responsive
    print_rule = css[css.index("@media print{"):]
    assert "workflow-flow" in print_rule


def test_abundance_intro_explains_reads_route_scope(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    assert "仅适用于 reads 输入路线" in zh
    assert "回比" in zh and "RPM" in zh and "平均深度" in zh and "覆盖度" in zh
    _, en = render_fixture(tmp_path, "en")
    assert "applies only to the reads input route" in en


def test_links_zh_visible_text_avoids_rank_and_top5(tmp_path):
    _, zh = render_fixture(tmp_path, "zh")
    links = re.search(r'<section id="links".*?</section>', zh, re.S).group(0)
    stripped = re.sub(r"<!--.*?-->", "", links)
    visible = re.sub(r"<[^>]+>", " ", stripped)
    assert "top-5" not in visible and "rank" not in visible
    assert "前 5 个近邻" in visible and "排序位次" in visible
    _, en = render_fixture(tmp_path, "en")
    en_links = re.search(r'<section id="links".*?</section>', en, re.S).group(0)
    assert "Structured top-5 neighbor evidence" in en_links


def test_warnings_section_deduplicates_repeated_notices(tmp_path):
    context, rows, summary = rich_fixture("reads")
    summary["warnings"] = [
        "21 条代表序列的骨架树正在建设中，暂未纳入系统发育放置",
        "21 条代表序列的骨架树正在建设中，暂未纳入系统发育放置",
        "另一条警告",
    ]
    zh = render_html(context, rows, summary, lang="zh")
    assert zh.count("21 条代表序列的骨架树正在建设中") == 1
    assert zh.count("另一条警告") == 1
    context, rows, summary = rich_fixture("reads")
    summary["warnings"] = []
    summary["warnings_en"] = []
    zh = render_html(context, rows, summary, lang="zh")
    assert "<h3>警告</h3>" not in zh
