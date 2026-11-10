"""最终结果归集、静态图表和 Markdown/JSON 汇总。"""

from __future__ import annotations

import csv
import html
import shutil
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cdc_pipeline.html_report import build_html_report
from cdc_pipeline.pipeline import PipelineContext, StageResult
from cdc_pipeline.config import PROJECT_ROOT
from cdc_pipeline.stages.common import artifact_path, read_tsv, write_tsv
from cdc_pipeline.state import atomic_write_json

BLUE = "#2a78d6"
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
SECONDARY = "#52514e"
GRID = "#e3e2df"
DELIVERY_MANIFEST = "delivery_manifest.json"
DELIVERY_MANIFEST_FORMAT_VERSION = 1
DELIVERY_SUMMARY_FIELDS = (
    "run_id", "route", "sample", "counts", "risk_rows", "risk_summary",
    "supergroup_stats", "multi_supergroup_contigs", "multi_supergroup_contig_count",
    "abundance_rows",
    "funnel", "stage_metrics", "timings", "total_minutes", "supergroup_outcome",
    "links", "warnings", "warnings_en",
)


def _delivery_summary(summary: dict) -> dict:
    """交付视图仅暴露业务字段；内部 state snapshot、绝对路径留在阶段审计。"""
    return {key: summary[key] for key in DELIVERY_SUMMARY_FIELDS if key in summary}


def _copy_delivery(source: str | Path, target: Path) -> Path:
    source = Path(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != target.resolve():
        shutil.copy2(source, target)
    return target


def _placed_orf_mappings(path: str | Path | None, placed_ids: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    """单次流式扫描大映射表，仅保留已放置代表序列的原始 ORF 与 contig。"""
    original_ids: dict[str, str] = {}
    contig_map: dict[str, str] = {}
    if not path or not Path(path).is_file() or not placed_ids:
        return original_ids, contig_map
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query_id = row.get("query_id", "")
            if query_id not in placed_ids:
                continue
            original_ids[query_id] = row.get("orf_id") or row.get("original_orf_id") or ""
            contig_map[query_id] = row.get("contig_id", "")
    return original_ids, contig_map


def _blank(value):
    return "" if value is None else value


def _flag(value):
    return "" if value is None else int(value)


def backbone_pending_warning(pending_merge: int, no_tree: int) -> tuple[str, str] | None:
    """骨架树在建警告：合并为一条；合计为 0 时不产生警告，避免重复提示。"""
    total = int(pending_merge or 0) + int(no_tree or 0)
    if total <= 0:
        return None
    zh = f"{total} 条代表序列的骨架树正在建设中，暂未纳入系统发育放置"
    en = (
        f"{total} representative sequence(s) await backbone tree construction "
        "and are not included in phylogenetic placement yet"
    )
    return zh, en


_INTERNAL_NEAREST = "nn1"
_PUBLIC_HINT_SUBSTITUTIONS = (
    (f"（{_INTERNAL_NEAREST} 所属科）", "（第 1 近邻所属科）"),
    (f"{_INTERNAL_NEAREST} 所属科未识别", "第 1 近邻所属科未识别"),
    (f"{_INTERNAL_NEAREST} 所属科", "第 1 近邻所属科"),
    (_INTERNAL_NEAREST, "第 1 近邻"),
)


def _public_hint(hint: str) -> str:
    """把上游内部缩写提示改写为对外中文表述；上游模块不在本次整改范围。"""
    text = str(hint or "")
    for old, new in _PUBLIC_HINT_SUBSTITUTIONS:
        text = text.replace(old, new)
    return text


def _write_risk_tables(output_dir: str | Path, risk_rows: list) -> tuple[Path, Path]:
    """写风险提示宽表和逐真实近邻排序位次的长表；未知值保持为空。"""
    output_dir = Path(output_dir)
    risk_fields = [
        "query_id", "original_id", "contig_id", "supergroup", "nn1", "nn1_dist",
        "family", "genus", "host_human", "host_vertebrate", "host_plant",
        "pathogen_category", "pathogen_examples", "statutory_category", "statutory_examples", "hint",
    ]
    risk_records = []
    neighbor_records = []
    for row in risk_rows:
        risk_records.append({
            "query_id": row.query_id, "original_id": row.original_id, "contig_id": row.contig_id,
            "supergroup": row.supergroup, "nn1": row.nn1, "nn1_dist": _blank(row.nn1_dist),
            "family": row.family, "genus": row.genus,
            "host_human": _flag(row.host_human), "host_vertebrate": _flag(row.host_vertebrate),
            "host_plant": _flag(row.host_plant), "pathogen_category": row.pathogen_category,
            "pathogen_examples": row.pathogen_examples, "statutory_category": row.statutory_category,
            "statutory_examples": row.statutory_examples, "hint": row.hint,
        })
        for neighbor in row.neighbors:
            neighbor_records.append({
                "query_id": row.query_id, "rank": neighbor.rank, "reference": neighbor.reference_id,
                "distance": _blank(neighbor.distance), "family": neighbor.family, "genus": neighbor.genus,
                "human": _flag(neighbor.human), "vertebrate": _flag(neighbor.vertebrate),
                "plant": _flag(neighbor.plant),
            })
    risk_path = write_tsv(output_dir / "risk_warning.tsv", risk_fields, risk_records)
    neighbor_path = write_tsv(
        output_dir / "neighbor_evidence.tsv",
        ["query_id", "rank", "reference", "distance", "family", "genus", "human", "vertebrate", "plant"],
        neighbor_records,
    )
    return risk_path, neighbor_path


def _remove_known_file(path: Path) -> None:
    if path.is_file() or path.is_symlink():
        path.unlink()


def _remove_empty_dir(path: Path) -> None:
    try:
        path.rmdir()
    except OSError:
        pass


def _delivery_relative_files(result_dir: Path, paths) -> list[str]:
    """交付清单中的路径统一为 results 根目录下的相对 POSIX 路径。"""
    relative = set()
    root = result_dir.resolve()
    for path in paths:
        try:
            relative.add(Path(path).resolve().relative_to(root).as_posix())
        except ValueError:
            continue
    return sorted(relative)


def _cleanup_legacy_delivery(result_dir: Path) -> None:
    """新包完整生成后，仅删除已知的旧交付项，保留所有未知用户文件。"""
    legacy_root = {
        "prefilter_calls.tsv", "orf_predictions.tsv", "contig_level_calls.tsv",
        "rdrp_positive.faa", "rdrp_positive_contigs.fna", "representatives.faa",
        "representative_map.tsv", "motifs.tsv", "supergroup_assignments.tsv",
        "placements.tsv", "abundance.tsv",
        "supergroup_abundance.tsv", "placement_annotations.tsv", "risk_warning.tsv",
        "neighbor_evidence.tsv", "summary.md", "summary.json", "run_summary.json",
        "tree_explorer.html",
    }
    for name in legacy_root:
        _remove_known_file(result_dir / name)

    report_dir = result_dir / "report"
    legacy_report = {
        "main_results_report_zh.html", "main_results_report_en.html", "phylogenetic_explorer.html",
        "tree_explorer.html", "summary.md", "summary.json", "run_summary.json", "README.md",
        "supergroup_details.tsv", "placement_details.tsv", "stage_timings.tsv",
        "sequence_funnel.tsv", "stage_metrics.tsv", "risk_warning.tsv", "neighbor_evidence.tsv",
    }
    for name in legacy_report:
        _remove_known_file(report_dir / name)
    figure_names = {
        f"{stem}.{extension}"
        for stem in ("workflow_overview", "supergroup_counts", "placement_status", "abundance_by_supergroup")
        for extension in ("png", "svg")
    }
    for figure_dir in (report_dir / "figures", result_dir / "figures"):
        for name in figure_names:
            _remove_known_file(figure_dir / name)
        _remove_empty_dir(figure_dir)
    _remove_empty_dir(report_dir)


def _delivery_description(relative: str) -> str:
    name = Path(relative).name
    descriptions = {
        "main_results_report_zh.html": "中文主报告；单位：报告。<br>Chinese main report; unit: report.",
        "main_results_report_en.html": "英文主报告；单位：报告。<br>English main report; unit: report.",
        "phylogenetic_explorer.html": "离线系统发育树交互浏览器；单位：报告。<br>Offline interactive phylogenetic explorer; unit: report.",
        "README.md": "本双语交付索引；单位：索引。<br>This bilingual delivery index; unit: index.",
        "delivery_manifest.json": "交付完成哨兵与文件清单；单位：运行。<br>Delivery-complete sentinel and file inventory; unit: run.",
        "run_summary.json": "结构化运行汇总、风险行与近邻证据；单位：运行。<br>Structured run, risk and neighbor summary; unit: run.",
        "prefilter_calls.tsv": "蛋白粗筛判定；单位：蛋白 ORF。<br>Protein prefilter calls; unit: protein ORF.",
        "orf_predictions.tsv": "RdRP 识别模型逐 ORF 分数与判定；单位：蛋白 ORF。<br>Per-ORF RdRP model scores and calls; unit: protein ORF.",
        "contig_level_calls.tsv": "核酸 contig 级检出汇总；单位：contig。<br>Contig-level detection summary; unit: contig.",
        "representative_map.tsv": "代表与全部成员的映射；单位：蛋白 ORF。<br>Representative-to-member mapping; unit: protein ORF.",
        "motifs.tsv": "RdRP motif 位置；单位：蛋白 ORF。<br>RdRP motif positions; unit: protein ORF.",
        "supergroup_assignments.tsv": "超群分配及同源证据；单位：代表蛋白。<br>Supergroup assignments and evidence; unit: representative protein.",
        "placements.tsv": "系统发育放置的最近邻与距离；单位：代表蛋白。<br>Phylogenetic nearest neighbors and distances; unit: representative protein.",
        "abundance.tsv": "reads 回比、RPM、深度和覆盖度；单位：源 contig。<br>Read mapping, RPM, depth and coverage; unit: source contig.",
        "supergroup_abundance.tsv": "超群相对丰度汇总；单位：超群。<br>Relative abundance by supergroup; unit: supergroup.",
        "placement_annotations.tsv": "放置注释及前 5 个近邻；单位：代表蛋白。<br>Placement annotations and top-5 neighbors; unit: representative protein.",
        "risk_warning.tsv": "类群归属与名录提示，不构成风险判定；单位：代表蛋白。<br>Taxonomic and catalog notices, not risk calls; unit: representative protein.",
        "neighbor_evidence.tsv": "每条代表序列的每个真实近邻排序位次一行；单位：近邻证据。<br>One row per real neighbor rank per query; unit: neighbor evidence.",
        "supergroup_details.tsv": "超群分类明细；单位：代表蛋白。<br>Supergroup classification details; unit: representative protein.",
        "placement_details.tsv": "系统发育放置明细；单位：代表蛋白。<br>Phylogenetic placement details; unit: representative protein.",
        "stage_timings.tsv": "各阶段状态与耗时；单位：阶段。<br>Stage status and timing; unit: stage.",
        "sequence_funnel.tsv": "各步骤序列数量与类型；单位：步骤。<br>Sequence counts and types by step; unit: step.",
        "stage_metrics.tsv": "各阶段输入、输出和耗时；单位：阶段。<br>Stage inputs, outputs and timing; unit: stage.",
        "rdrp_positive_orfs.faa": "RdRP 阳性 ORF 蛋白序列；单位：蛋白 ORF。<br>RdRP-positive ORF proteins; unit: protein ORF.",
        "rdrp_positive_contigs.fna": "RdRP 阳性源核酸 contigs；单位：contig。<br>RdRP-positive source nucleotide contigs; unit: contig.",
        "representative_orfs.faa": "去冗余代表 ORF 蛋白序列；单位：代表蛋白。<br>Deduplicated representative ORF proteins; unit: representative protein.",
    }
    if relative.startswith("trees/"):
        return "本次可视化使用的骨架树 Newick；单位：树。<br>Backbone Newick used by this visualization; unit: tree."
    return descriptions.get(name, "交付数据文件。<br>Delivered data file.")


def _write_delivery_readme(
    result_dir: Path,
    delivered_files: list[Path],
    warnings: list[str] | None = None,
    warnings_en: list[str] | None = None,
) -> Path:
    output = result_dir / "README.md"
    files = {output.resolve(): output, **{path.resolve(): path for path in delivered_files}}
    lines = [
        "# 结果文件索引 / File index",
        "",
        "本目录包含本次运行的最终交付文件；中间汇总和静态图保留在 report 阶段目录。",
        "This directory contains the final delivery for this run. Intermediate summaries and static figures remain in the report stage directory.",
        "风险提示只描述骨架 taxonomy 的类群归属与名录先验，不构成新序列的风险分级或病原确证。",
        "Risk notices describe backbone-taxonomy group membership and catalog priors only; they are not risk classifications or pathogen confirmations for new sequences.",
    ]
    warning_lines = list(dict.fromkeys(
        [f"- {warning}" for warning in warnings or []]
        + [f"- {warning}" for warning in warnings_en or []]
    ))
    if warning_lines:
        lines += ["", "## 警告与限制 / Warnings and limitations", "", *warning_lines]
    lines += [
        "",
        "| 文件 / File | 内容与单位 / Content and unit |",
        "|---|---|",
    ]
    for path in sorted(files.values(), key=lambda item: item.relative_to(result_dir).as_posix()):
        relative = path.relative_to(result_dir).as_posix()
        escaped = html.escape(relative, quote=True)
        lines.append(f'| <a href="{escaped}" target="_blank" rel="noopener">{escaped}</a> | {_delivery_description(relative)} |')
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def plot_bars(
    output_prefix: str | Path,
    categories: list[str],
    values: list[float],
    title: str,
    value_label: str,
    config: dict,
    horizontal: bool = False,
) -> list[Path]:
    """单一系列条形图；单类别由表格表达，避免无意义的一根柱。"""
    if len(categories) < 2:
        return []
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "text.color": TEXT,
        "axes.labelcolor": SECONDARY,
        "xtick.color": SECONDARY,
        "ytick.color": SECONDARY,
        "svg.fonttype": "none",
    })
    width = max(5.5, min(10, 1.2 * len(categories) + 3))
    height = max(3.2, 0.5 * len(categories) + 1.8) if horizontal else 3.8
    figure, axis = plt.subplots(figsize=(width, height), facecolor=SURFACE)
    axis.set_facecolor(SURFACE)
    positions = list(range(len(categories)))
    if horizontal:
        axis.barh(positions, values, height=0.55, color=BLUE, edgecolor="none")
        axis.set_yticks(positions, categories)
        axis.invert_yaxis()
        axis.set_xlabel(value_label)
        axis.xaxis.grid(True, color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        maximum = max(values) if values else 0
        axis.set_xlim(0, maximum * 1.22 if maximum else 1)
        for position, value in zip(positions, values):
            axis.text(value + maximum * 0.02, position, f"{value:.3g}", va="center", color=TEXT)
    else:
        axis.bar(positions, values, width=0.55, color=BLUE, edgecolor="none")
        axis.set_xticks(positions, categories, rotation=18 if max(map(len, categories)) > 12 else 0, ha="right" if max(map(len, categories)) > 12 else "center")
        axis.set_ylabel(value_label)
        axis.yaxis.grid(True, color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        maximum = max(values) if values else 0
        axis.set_ylim(0, maximum * 1.22 if maximum else 1)
        for position, value in zip(positions, values):
            axis.text(position, value + maximum * 0.02, f"{value:.3g}", ha="center", color=TEXT)
    axis.set_title(title, loc="left", fontsize=13, fontweight="semibold", pad=12)
    for side in ["top", "right", "left"]:
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color(GRID)
    axis.spines["bottom"].set_linewidth(0.8)
    axis.tick_params(axis="both", length=0)
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for fmt in config["report"]["formats"]:
        path = output_prefix.with_suffix(f".{fmt}")
        figure.savefig(path, dpi=config["report"]["dpi"], bbox_inches="tight", facecolor=SURFACE)
        outputs.append(path)
    plt.close(figure)
    return outputs


def summary_markdown(summary: dict) -> str:
    counts = summary["counts"]
    lines = [
        "# lpvid RNA 病毒 RdRP 分析报告",
        "",
        f"- 运行：`{summary['run_id']}`",
        f"- 样本：`{summary['sample']}`",
        f"- 输入路线：`{summary['route']}`",
        f"- **总耗时：{summary.get('total_minutes', 0):.1f} 分钟**",
        "",
        "## 一、分析步骤与序列数量",
        "",
        "| 处理步骤 | 状态 | 输入 | 输出 | 耗时（分钟） |",
        "|---|---|---|---|---:|",
    ]
    for row in summary.get("stage_metrics", []):
        lines.append(
            f"| {row.get('stage','')} | {row.get('status','')} | "
            f"{row.get('input_count',0):,} {row.get('input_unit','')} | "
            f"{row.get('output_count',0):,} {row.get('output_unit','')} | "
            f"{row.get('duration_minutes',0):.2f} |"
        )
    lines += [
        f"| **总计** | — | — | — | **{summary.get('total_minutes', 0):.2f}** |",
        "",
        "## 二、序列数量漏斗（明确序列类型）",
        "",
        "| 处理步骤 | 数量 | 单位 | 说明 |",
        "|---|---:|---|---|",
    ]
    for f in summary.get("funnel", []):
        lines.append(f"| {f['stage']} | {f['count']:,} | {f.get('unit','')} | {f.get('note','')} |")
    lines += [
        "",
        "## 三、核心结果",
        "",
        "| 指标 | 数量 |",
        "|---|---:|",
        f"| RdRP 阳性蛋白 ORF（主分析集） | {counts.get('positive_orfs', 0)} |",
        f"| 去冗余代表蛋白 | {counts.get('representatives', 0)} |",
        f"| 已分类代表蛋白 | {counts.get('classified', 0)} |",
        f"| 有可用骨架树的代表蛋白 | {counts.get('tree_ready', 0)} |",
        f"| 骨架树建设中的代表蛋白 | {counts.get('pending_tree', 0)} |",
        f"| 未分类代表蛋白 | {counts.get('unassigned', 0)} |",
        f"| 完成系统发育放置 | {counts.get('placed', 0)} |",
        f"| RdRP 阳性核酸 contigs | {counts.get('positive_contigs', 0)} |",
        "",
        "## 四、输出索引",
        "",
    ]
    for name, link in summary.get("links", {}).items():
        lines.append(f"- [{name}]({link})")
    lines += [
        "",
        "## 五、超群分配详情",
        "",
        "| 代表序列 | 超群 | 分类方法 | 状态 | 参考序列 | 比对 E-value | 一致性 (%) |",
        "|---|---|---|---|---|---|---:|",
    ]
    from cdc_pipeline.stage_metrics import extract_supergroup_details, extract_placement_details
    sg_details = extract_supergroup_details(summary.get("_run_dir", ""))
    pl_details = extract_placement_details(summary.get("_run_dir", ""))
    if sg_details:
        for d in sg_details[:30]:
            lines.append(
                f"| {d['query_id']} | {d['supergroup']} | {d.get('method','')} "
                f"| {d['status']} | {d.get('reference_id','')[:30]} "
                f"| {d.get('evalue','')} | {d.get('pident','')} |"
            )
        if len(sg_details) > 30:
            lines.append(f"| ... 共 {len(sg_details)} 条，详见 supergroup_details.tsv | | | | | | |")
    else:
        lines.append("| （无超群分配结果） | | | | | | |")

    if pl_details:
        lines += [
            "",
            "## 六、系统发育放置详情",
            "",
            "| Query | 超群 | 最近邻 | 距离 |",
            "|---|---|---|---|",
        ]
        for d in pl_details[:20]:
            lines.append(
                f"| {d['query_id']} | {d['supergroup']} | {d['nearest_neighbor'][:30]} "
                f"| {d['nn1_distance']} |"
            )
        if len(pl_details) > 20:
            lines.append(f"| ... 共 {len(pl_details)} 条，详见 placement_details.tsv | | | |")

    lines += [
        "",
        "## 警告、跳过与限制",
        "",
    ]
    warnings = summary.get("warnings") or ["无"]
    lines += [f"- {warning}" for warning in warnings]
    outcome = summary.get("supergroup_outcome", {})
    lines += [
        "",
        "## 五、结果解读（按分析顺序）",
        "",
    ]
    for item in summary.get("funnel", []):
        lines.append(f"- **{item['stage']}**：{item['count']:,} {item.get('unit','')}。{item.get('note','')}")
    lines += [
        f"- **骨架树状态**：有可用树 {outcome.get('tree_available', 0):,} 条；骨架树正在建设中 {outcome.get('pending_merge', 0) + outcome.get('no_tree', 0):,} 条；未分类 {outcome.get('unresolved', 0):,} 条。",
        "",
        "## 六、使用边界",
        "",
        "- 参考序列的一致性百分比（Pident）是同源比对指标，不是分类概率。",
        "- LucaPhylo 结果表示最近邻骨架序列及其预测距离，不代表枝级置信度。",
        "- 去冗余成员全部保留在 representative_map 中；超群和系统发育分析针对代表蛋白。",
    ]
    if summary["route"] == "reads":
        lines.append("- RPM 是文库内相对丰度；双端分母为 R1+R2 reads，不是 copies/mL 绝对载量。")
    lines += ["", "## 七、阶段状态", "", "| 阶段 | 状态 |", "|---|---|"]
    for name, stage in summary.get("stages", {}).items():
        lines.append(f"| {name} | {stage.get('status', '')} |")
    lines.append("")
    return "\n".join(lines)


def build_report(context: PipelineContext) -> dict:
    from cdc_pipeline.stage_metrics import extract_supergroup_outcome

    data = context.data
    detect = data.get("detect_counts", {})
    representatives = data.get("representative_counts", {})
    supergroups = data.get("supergroup_counts", {})
    phylo = data.get("phylo_counts", {})
    result_dir = context.run_dir / "results"
    tables_dir = result_dir / "tables"
    sequences_dir = result_dir / "sequences"
    trees_dir = result_dir / "trees"
    for directory in (tables_dir, sequences_dir, trees_dir):
        directory.mkdir(parents=True, exist_ok=True)

    manifest_path = tables_dir / DELIVERY_MANIFEST
    _remove_known_file(manifest_path)
    _remove_known_file(manifest_path.with_name(manifest_path.name + ".tmp"))
    source_files = {
        "prefilter_calls": ("tables", "prefilter_calls.tsv"),
        "orf_predictions": ("tables", "orf_predictions.tsv"),
        "contig_calls": ("tables", "contig_level_calls.tsv"),
        "positive_faa": ("sequences", "rdrp_positive_orfs.faa"),
        "positive_contigs_fna": ("sequences", "rdrp_positive_contigs.fna"),
        "representatives_faa": ("sequences", "representative_orfs.faa"),
        "representative_map": ("tables", "representative_map.tsv"),
        "motifs": ("tables", "motifs.tsv"),
        "supergroup_assignments": ("tables", "supergroup_assignments.tsv"),
        "placements_tsv": ("tables", "placements.tsv"),
        "abundance_tsv": ("tables", "abundance.tsv"),
        "supergroup_abundance_tsv": ("tables", "supergroup_abundance.tsv"),
        "placement_annotations": ("tables", "placement_annotations.tsv"),
    }
    links = {}
    result_files = {}
    for key, (category, name) in source_files.items():
        source = data.get(key)
        if source and Path(source).is_file():
            target = _copy_delivery(source, result_dir / category / name)
            result_files[key] = target
            links[name] = f"{category}/{name}"
    if data.get("tree_html") and Path(data["tree_html"]).is_file():



        from scripts.verify_tree_explorer import check_file as check_tree_contract

        tree_problems = check_tree_contract(
            Path(data["tree_html"]),
            require_original_ids=bool(data.get("orf_contig_map")),
        )
        if tree_problems:
            raise RuntimeError(
                "树浏览器为旧产物或与当前 schema 不一致，请先重跑 visualization: "
                + "; ".join(tree_problems[:5])
            )
        tree_target = _copy_delivery(data["tree_html"], result_dir / "phylogenetic_explorer.html")
        result_files["phylo_html"] = tree_target
        links["phylogenetic_explorer.html"] = "phylogenetic_explorer.html"
    delivered_trees = []
    for supergroup, tree_path in (data.get("tree_files") or {}).items():
        tree_path = Path(tree_path)
        if tree_path.is_file():
            target = _copy_delivery(tree_path, trees_dir / tree_path.name)
            delivered_trees.append(target)
            result_files[f"tree_{supergroup}"] = target
            links[f"trees/{target.name}"] = f"trees/{target.name}"

    outcome = extract_supergroup_outcome(context.run_dir)
    counts = {
        "positive_orfs": detect.get("main_rdrp_positive_orfs", 0),
        "representatives": representatives.get("representatives", 0),
        "classified": outcome["representatives"] - outcome["unresolved"],
        "tree_ready": outcome["tree_available"],
        "pending_tree": outcome["pending_merge"] + outcome["no_tree"],
        "placed": phylo.get("placed", 0),
        "unassigned": outcome["unresolved"],
        "positive_contigs": detect.get("positive_contigs", 0),
    }
    warnings = []
    warnings_en = []
    phylo_errors = [f"{sg}: {error}" for sg, error in data.get("phylo_errors", {}).items()]
    warnings += phylo_errors
    warnings_en += phylo_errors
    pending_warning = backbone_pending_warning(outcome["pending_merge"], outcome["no_tree"])
    if pending_warning:
        warnings.append(pending_warning[0])
        warnings_en.append(pending_warning[1])
    if outcome["unresolved"]:
        warnings.append(f"{outcome['unresolved']} 条代表序列未能分配超群")
        warnings_en.append(f"{outcome['unresolved']} representative(s) could not be assigned to a supergroup")
    if detect.get("qa_rdrp_positive_orfs"):
        warnings.append(f"DISCARD QA 中发现 {detect['qa_rdrp_positive_orfs']} 条阳性 ORF；未混入主结果")
        warnings_en.append(f"{detect['qa_rdrp_positive_orfs']} positive ORF(s) found in DISCARD QA; kept out of the main results")


    from cdc_pipeline.risk_assessment import build_risk_rows, load_pathogen_catalog, load_taxonomy, summarize_risk
    risk_rows = []
    risk_summary = {}
    if result_files.get("placements_tsv"):
        phylo_root = context.config.get("assets", {}).get("phylo_root", str(PROJECT_ROOT / "resources" / "phylo"))
        taxonomy = load_taxonomy(Path(phylo_root) / "taxonomy")
        catalog = load_pathogen_catalog(PROJECT_ROOT / "assets" / "risk" / "china_pathogen_catalog.csv")
        placements = read_tsv(result_files["placements_tsv"])
        placed_ids = {row.get("query_id", "") for row in placements if row.get("query_id")}
        original_ids, contig_map = _placed_orf_mappings(data.get("orf_contig_map"), placed_ids)
        risk_rows = build_risk_rows(
            placements, taxonomy, catalog,
            original_ids=original_ids, contig_map=contig_map,
        )
        risk_rows = [replace(row, hint=_public_hint(row.hint)) for row in risk_rows]
        risk_summary = summarize_risk(risk_rows)
    risk_tsv, neighbor_tsv = _write_risk_tables(tables_dir, risk_rows)
    result_files["risk_warning"] = risk_tsv
    result_files["neighbor_evidence"] = neighbor_tsv
    links["risk_warning.tsv"] = "tables/risk_warning.tsv"
    links["neighbor_evidence.tsv"] = "tables/neighbor_evidence.tsv"

    figures = []
    figures += plot_bars(
        context.stage_dir / "figures" / "workflow_overview",
        ["Positive ORFs", "Representatives", "Classified", "Placed"],
        [counts["positive_orfs"], counts["representatives"], counts["classified"], counts["placed"]],
        "RdRP workflow overview", "Sequences", context.config,
    )
    sg_counts = {key: value for key, value in supergroups.items() if key not in {"unassigned", "representatives"}}
    sg_counts = dict(sorted(sg_counts.items(), key=lambda item: (-item[1], item[0]))[:10])
    figures += plot_bars(
        context.stage_dir / "figures" / "supergroup_counts",
        list(sg_counts), list(sg_counts.values()), "Representatives by supergroup", "Representative sequences",
        context.config, horizontal=True,
    )
    figures += plot_bars(
        context.stage_dir / "figures" / "placement_status",
        ["Placed", "Unassigned"], [counts["placed"], counts["unassigned"]],
        "Placement status", "Representative sequences", context.config,
    )
    if data.get("supergroup_abundance_tsv") and Path(data["supergroup_abundance_tsv"]).is_file():
        abundance_rows = read_tsv(data["supergroup_abundance_tsv"])
        abundance_rows = sorted(abundance_rows, key=lambda row: -float(row["RPM"]))[:10]
        figures += plot_bars(
            context.stage_dir / "figures" / "abundance_by_supergroup",
            [row["supergroup"] for row in abundance_rows], [float(row["RPM"]) for row in abundance_rows],
            "Relative abundance by supergroup", "RPM", context.config, horizontal=True,
        )
    stages = context.state.snapshot.get("stages", {})
    stages["report"] = {**stages.get("report", {}), "status": "done", "error": None}

    summary = {
        "run_id": context.run_dir.name,
        "route": context.route,
        "sample": context.sample,
        "counts": counts,
        "risk_rows": risk_rows,
        "risk_summary": risk_summary,
        "warnings": warnings,
        "warnings_en": warnings_en,
        "links": links,
        "stages": stages,
        "figures": [str(path) for path in figures],
    }
    summary["links"]["main_results_report_zh.html"] = "main_results_report_zh.html"
    summary["links"]["main_results_report_en.html"] = "main_results_report_en.html"
    summary["links"]["run_summary.json"] = "tables/run_summary.json"
    summary["links"]["delivery_manifest.json"] = f"tables/{DELIVERY_MANIFEST}"
    summary["links"]["README.md"] = "README.md"


    from cdc_pipeline.stage_metrics import extract_stage_metrics, extract_stage_timings, extract_sequence_funnel
    timings, total_minutes = extract_stage_timings(context.run_dir)
    funnel = extract_sequence_funnel(context.run_dir)
    stage_metrics = extract_stage_metrics(context.run_dir, timings)
    summary["timings"] = timings
    summary["total_minutes"] = total_minutes
    summary["funnel"] = funnel
    summary["stage_metrics"] = stage_metrics
    summary["supergroup_outcome"] = outcome
    summary["_run_dir"] = str(context.run_dir)


    from cdc_pipeline.stage_metrics import extract_supergroup_details, extract_placement_details
    sg_details = extract_supergroup_details(context.run_dir)
    pl_details = extract_placement_details(context.run_dir)

    generated_details = set()
    sg_tsv = context.stage_dir / "supergroup_details.tsv"
    if sg_details:
        write_tsv(sg_tsv, list(sg_details[0].keys()), sg_details)
        result_files["supergroup_details"] = sg_tsv
        generated_details.add("supergroup_details")

    pl_tsv = context.stage_dir / "placement_details.tsv"
    if pl_details:
        write_tsv(pl_tsv, list(pl_details[0].keys()), pl_details)
        result_files["placement_details"] = pl_tsv
        generated_details.add("placement_details")


    timing_tsv = context.stage_dir / "stage_timings.tsv"
    write_tsv(timing_tsv, ["stage", "status", "duration_minutes"], timings)
    result_files["stage_timings"] = timing_tsv


    funnel_tsv = context.stage_dir / "sequence_funnel.tsv"
    write_tsv(funnel_tsv, ["stage", "count", "unit", "note"], funnel)
    result_files["sequence_funnel"] = funnel_tsv


    stage_metrics_tsv = context.stage_dir / "stage_metrics.tsv"
    write_tsv(stage_metrics_tsv, [
        "stage_key", "stage", "status", "input_count", "input_unit",
        "output_count", "output_unit", "duration_minutes", "note",
    ], stage_metrics)
    result_files["stage_metrics"] = stage_metrics_tsv


    for key, source, name in [
        ("supergroup_details", sg_tsv, "supergroup_details.tsv"),
        ("placement_details", pl_tsv, "placement_details.tsv"),
        ("stage_timings", timing_tsv, "stage_timings.tsv"),
        ("sequence_funnel", funnel_tsv, "sequence_funnel.tsv"),
        ("stage_metrics", stage_metrics_tsv, "stage_metrics.tsv"),
    ]:

        if key in {"supergroup_details", "placement_details"} and key not in generated_details:
            continue
        if source and Path(source).is_file():
            target = _copy_delivery(source, tables_dir / name)
            result_files[key] = target
            summary["links"][name] = f"tables/{name}"


    report_zh = build_html_report(context, result_dir / "main_results_report_zh.html", summary, lang="zh")
    report_en = build_html_report(context, result_dir / "main_results_report_en.html", summary, lang="en")
    result_files["report_zh"] = report_zh
    result_files["report_en"] = report_en

    markdown = context.stage_dir / "summary.md"
    markdown.write_text(summary_markdown(summary), encoding="utf-8")
    summary_json = context.stage_dir / "summary.json"
    atomic_write_json(summary_json, summary)
    delivery_summary = tables_dir / "run_summary.json"
    atomic_write_json(delivery_summary, _delivery_summary(summary))
    result_files["run_summary_json"] = delivery_summary
    result_files["delivery_manifest"] = manifest_path
    readme = _write_delivery_readme(
        result_dir,
        list(result_files.values()),
        warnings=summary["warnings"],
        warnings_en=summary["warnings_en"],
    )
    result_files["readme"] = readme


    _cleanup_legacy_delivery(result_dir)

    for key, name in (
        ("supergroup_details", "supergroup_details.tsv"),
        ("placement_details", "placement_details.tsv"),
    ):
        if key not in generated_details:
            _remove_known_file(context.stage_dir / name)
            _remove_known_file(tables_dir / name)


    manifest_files = _delivery_relative_files(result_dir, [readme, *result_files.values()])
    atomic_write_json(manifest_path, {
        "status": "complete",
        "format_version": DELIVERY_MANIFEST_FORMAT_VERSION,
        "run_id": context.run_dir.name,
        "files": manifest_files,
    })

    return {"summary": summary, "markdown": markdown, "summary_json": summary_json, "figures": figures, "result_files": result_files, "report_zh": report_zh, "report_en": report_en}


def run(context: PipelineContext) -> StageResult:
    result = build_report(context)
    return StageResult(
        artifacts={
            "summary_md": artifact_path(context, result["markdown"]),
            "summary_json": artifact_path(context, result["summary_json"]),
            "tree_html": artifact_path(context, result["result_files"]["phylo_html"]) if "phylo_html" in result["result_files"] else None,
            "report_zh": artifact_path(context, result["report_zh"]),
            "report_en": artifact_path(context, result["report_en"]),
            "figures": [artifact_path(context, path) for path in result["figures"]],
            "results": {key: artifact_path(context, path) for key, path in result["result_files"].items()},
        },
        summary=result["summary"]["counts"],
        data={"summary_md": str(result["markdown"]), "summary_json": str(result["summary_json"]), "report_zh": str(result["report_zh"]), "report_en": str(result["report_en"]), "report_figures": [str(path) for path in result["figures"]]},
    )
