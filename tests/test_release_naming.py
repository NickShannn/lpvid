"""Release-wide naming contract for the public distribution."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import re

import yaml

from cdc_pipeline.config import PROJECT_ROOT, load_config
from cdc_pipeline.pipeline import stages_for

TEXT_SUFFIXES = {
    ".py", ".yaml", ".yml", ".json", ".md", ".tsv", ".txt", ".csv",
    ".sh", ".html", ".js", ".fai", ".fa", ".faa", ".fna",
}
SKIP_DIRS = {"third_party", ".git", ".pytest_cache", "__pycache__"}
SKIP_FILES = {"test_release_naming.py"}


def _term(*parts: str) -> str:
    return "".join(parts)


def _forbidden_terms() -> list[str]:
    return [
        _term("luca", "prot"),
        _term("luca", "virus"),
        _term("embed", "ding"),
        _term("m", "lp"),
        _term("stage", "-2"),
        _term("stage", "_2"),
        _term("top", "-hit"),
        _term("双", "模型"),
        _term("对外", "演示"),
        _term("演", "示版"),
        _term("luca", "host"),
        _term("luca", "pathogen"),
        _term("supp", "lementary"),
    ]


def _text_files():
    for path in sorted(PROJECT_ROOT.rglob("*")):
        if not path.is_file() or path.name in SKIP_FILES:
            continue
        if any(domain in path.parts for domain in SKIP_DIRS):
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            yield path, path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue


def _legacy_hits() -> list[str]:
    hits = []
    for path, text in _text_files():
        lower = text.lower()
        for term in _forbidden_terms():
            if term == _term("m", "lp"):
                tokens = [token.strip(".,;:()[]{}<>'\"`") for token in lower.replace("_", " ").replace("-", " ").split()]
                if term in tokens:
                    hits.append(f"{path.relative_to(PROJECT_ROOT)}: {term}")
                continue
            if term in lower:
                hits.append(f"{path.relative_to(PROJECT_ROOT)}: {term}")
    return hits


def test_no_first_party_file_uses_retired_model_or_internal_scheme_terms():
    hits = _legacy_hits()
    assert hits == [], "\n".join(hits)


def test_retired_files_and_old_install_names_are_removed():
    retired = [
        _term("luca", "prot_runner.py"),
        _term("luca", "virus_runner.py"),
        _term("embed", "ding_classifier.py"),
        _term("test_luca", "prot_runner.py"),
        _term("test_luca", "virus_runner.py"),
        _term("prepare_luca", "virus_compat.py"),
        _term("supp", "lementary_classifier.py"),
        "RESOURCE_PACKAGE_MANIFEST.json",
    ]
    for name in retired:
        assert not list(PROJECT_ROOT.rglob(name)), name
    assert not list(PROJECT_ROOT.rglob(_term("resources/models/", "luca", "prot")))
    assert not list(PROJECT_ROOT.rglob(_term("resources/models/", "luca", "virus", "*")))


def test_generic_stage_modules_and_tests_exist():
    stage_dir = PROJECT_ROOT / "cdc_pipeline" / "stages"
    assert (stage_dir / "rdrp_runner.py").is_file()
    assert (stage_dir / "classifier.py").is_file()
    assert (PROJECT_ROOT / "scripts" / "prepare_runtime_libs.py").is_file()
    assert (PROJECT_ROOT / "tests" / "test_rdrp_runner.py").is_file()


def test_default_and_local_configs_use_generic_names():
    default = yaml.safe_load((PROJECT_ROOT / "config" / "default.yaml").read_text(encoding="utf-8"))
    local = yaml.safe_load((PROJECT_ROOT / "config" / "local.example.yaml").read_text(encoding="utf-8"))
    assert default["detect"]["rdrp_truncation"] == 4096
    for raw in (default, local):
        assert set(raw["environments"]) == {"lucaphylo", "qc", "rdrp"}
        assert raw["assets"]["rdrp_model_dir"] == "resources/models/rdrp_detector"
        assert raw["assets"]["runtime_lib_dir"] == "runtime/libs"
        assert raw["assets"]["rdrp_180_db"] == "resources/rdrp/rdrp_reference.dmnd"
        assert raw["assets"]["rdrp_180_csv"] == "resources/rdrp/rdrp_reference.csv"
        assert raw["assets"]["classifier_pipeline"] == "resources/models/classifier/pipeline"
        assert raw["assets"]["classifier_encoder"] == "resources/models/classifier/encoder/checkpoint-18000"
        assert raw["assets"]["classifier_model"] == "resources/models/classifier/model/Ref_Tree_Orf_All_orf_tree_processed_Pool_Mean"
    config = load_config(PROJECT_ROOT / "config" / "default.yaml")
    assert Path(config["assets"]["runtime_lib_dir"]) == PROJECT_ROOT / "runtime" / "libs"
    assert Path(config["assets"]["classifier_pipeline"]).is_absolute()


def test_pipeline_routes_use_the_generic_classifier_stage_name():
    for route in ("proteins", "contigs", "reads"):
        stages = stages_for(route)
        assert "classifier" in stages
        assert _term("embed", "ding_classifier") not in stages
        assert _term("supp", "lementary_classifier") not in stages


def test_resource_manifest_contains_only_generic_public_paths():
    manifest = PROJECT_ROOT / "resources" / "RESOURCE_MANIFEST.tsv"
    text = manifest.read_text(encoding="utf-8")
    lower = text.lower()
    for term in _forbidden_terms():
        assert term not in lower, term
    assert "RESOURCE_PACKAGE_MANIFEST" not in text
    assert "resources/rdrp/rdrp_reference.dmnd" in text
    assert "resources/rdrp/rdrp_reference.fas" in text
    assert "resources/rdrp/rdrp_reference.csv" in text
    assert "resources/models/rdrp_detector" in text
    assert "resources/models/classifier/pipeline" in text
    assert "resources/models/classifier/encoder/checkpoint-18000" in text
    assert "resources/models/classifier/model/Ref_Tree_Orf_All_orf_tree_processed_Pool_Mean" in text


def test_small_runtime_resources_are_bundled_with_checksums():
    bundled_files = [
        "resources/rdrp/rdrp_reference.dmnd",
        "resources/rdrp/rdrp_reference.fas",
        "resources/rdrp/rdrp_reference.csv",
        "resources/phylo/backbone_dir/Nido/Nido.fas",
        "resources/phylo/taxonomy/Nido.csv",
        "resources/models/classifier/pipeline/classify.py",
        "resources/models/classifier/model/Ref_Tree_Orf_All_orf_tree_processed_Pool_Mean/model.pt",
        "resources/blastdb/refseq_protein/TAXID_viruses.txt",
        "resources/blastdb/refseq_protein/TAXID_cellular.txt",
    ]
    for relative in bundled_files:
        assert (PROJECT_ROOT / relative).is_file(), relative

    with (PROJECT_ROOT / "resources/RESOURCE_MANIFEST.tsv").open(newline="") as handle:
        rows = {row["install_path"]: row for row in csv.DictReader(handle, delimiter="\t")}
    for relative in bundled_files:
        row = rows[relative]
        assert row["source_hint"] == "bundled in repository", relative
        assert len(row["sha256"]) == 64, relative

    assert not (PROJECT_ROOT / "resources/phylo/adapter_dir").exists()
    assert not (PROJECT_ROOT / "resources/phylo/taxonomy/na-daxi.csv").exists()
    assert not list((PROJECT_ROOT / "resources").rglob("*.bak"))
    assert not (PROJECT_ROOT / "resources/rrna").exists()
    assert not (PROJECT_ROOT / "resources/models/rdrp_detector").exists()
    assert not (PROJECT_ROOT / "resources/models/lucaphylo").exists()
    assert not (PROJECT_ROOT / "resources/models/classifier/encoder").exists()
    assert not (PROJECT_ROOT / "resources/blastdb/refseq_protein/refseq_protein.dmnd").exists()


def test_release_uses_apache2_and_exact_dependency_versions():
    assert (PROJECT_ROOT / "LICENSE").read_text(encoding="utf-8").startswith("Apache License")
    for document in ("README.md", "README.zh-CN.md", "docs/INSTALL.md", "docs/INSTALL.zh-CN.md"):
        assert "Apache License 2.0" in (PROJECT_ROOT / document).read_text(encoding="utf-8"), document

    environment = yaml.safe_load((PROJECT_ROOT / "environment.yml").read_text(encoding="utf-8"))
    conda_packages = [item for item in environment["dependencies"] if isinstance(item, str)]
    for package in conda_packages:
        if package == "pip":
            continue
        assert "=" in package, package
    pip_packages = next(item for item in environment["dependencies"] if isinstance(item, dict))["pip"]
    for package in pip_packages:
        assert "==" in package, package

    for line in (PROJECT_ROOT / "config/requirements-lpvid.txt").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        assert "==" in line, line


def _manifest_module():
    path = PROJECT_ROOT / "resources" / "src" / "build_resource_manifest.py"
    spec = importlib.util.spec_from_file_location("build_resource_manifest_under_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_manifest_generator_scans_new_layout_without_the_internal_package_index():
    source = (PROJECT_ROOT / "resources" / "src" / "build_resource_manifest.py").read_text(encoding="utf-8")
    lower = source.lower()
    for term in _forbidden_terms():
        assert term not in lower, term
    assert "RESOURCE_PACKAGE_MANIFEST" not in source
    assert "assets/RdRP_blast_ref" not in source
    for expected in (
        "resources/rdrp/rdrp_reference.dmnd",
        "resources/rdrp/rdrp_reference.fas",
        "resources/rdrp/rdrp_reference.csv",
        "resources/models/rdrp_detector",
        "resources/models/classifier/pipeline",
        "resources/models/classifier/encoder",
        "resources/models/classifier/model",
    ):
        assert expected in source, expected


def test_manifest_generator_classifies_new_paths_directly():
    module = _manifest_module()
    assert module.classify_file("resources/rdrp/rdrp_reference.dmnd")[0] == "rdrp_ref_db"
    assert module.classify_file("resources/rdrp/rdrp_reference.fas")[0] == "rdrp_ref_db"
    assert module.classify_file("resources/rdrp/rdrp_reference.csv")[0] == "rdrp_ref_db"
    assert module.classify_file("resources/phylo/backbone_dir/Nido/Nido.fas")[0] == "lucaphylo_backbone"
    assert module.classify_bundle("resources/models/rdrp_detector")[0] == "rdrp_model"
    assert module.classify_bundle("resources/models/classifier/pipeline")[0] == "classifier_pipeline"
    assert module.classify_bundle("resources/models/classifier/encoder")[0] == "classifier_encoder"
    assert module.classify_bundle("resources/models/classifier/model")[0] == "classifier_model"
    assert module.classify_file("assets/RdRP_blast_ref/legacy.dmnd") is None
    assert module.classify_bundle("resources/models/unknown_model") is None


def test_report_sources_keep_latin_abbreviations_out_of_chinese_prose():
    """中文文案里不得出现 nn1/query/rank/top-5 这类拉丁缩写（数据字段名除外）。"""
    patterns = (
        r"[一-鿿]\s*(?:nn\d|query|rank|top-5)",
        r"(?:nn\d|query|rank|top-5)\s*[一-鿿]",
    )
    for relative in (
        Path("cdc_pipeline") / "html_template_v2.py",
        Path("cdc_pipeline") / "stages" / "report.py",
    ):
        text = (PROJECT_ROOT / relative).read_text(encoding="utf-8")
        # HTML 注释中的兼容校验锚点不参与渲染，不计入中文正文。
        text = re.sub(r"<!--.*?-->", "", text)
        for pattern in patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            assert match is None, f"{relative}: {match.group(0) if match else ''}"
