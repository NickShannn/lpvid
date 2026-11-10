from pathlib import Path

import yaml

from cdc_pipeline.config import DEFAULT_CONFIG, PROJECT_ROOT, deep_merge, load_config


def test_default_config_has_confirmed_scientific_defaults():
    config = load_config(DEFAULT_CONFIG)
    assert config["min_input_aa"] == 150
    assert config["min_analysis_aa"] == 150
    assert config["model_threshold"] == 0.5
    assert config["prefilter"]["enabled"] is True
    assert config["prefilter"]["qa_fraction"] == 0.01
    assert config["cluster"]["identity"] == 0.9
    assert config["cluster"]["coverage"] == 0.8
    assert config["phylo"]["max_seq_length"] == 2048
    assert config["phylo"]["top_k"] == 5


def test_deep_merge_does_not_replace_unrelated_nested_values():
    base = {"a": {"b": 1, "c": 2}, "d": 3}
    merged = deep_merge(base, {"a": {"b": 9}})
    assert merged == {"a": {"b": 9, "c": 2}, "d": 3}
    assert base["a"]["b"] == 1


def test_load_config_resolves_relative_paths(tmp_path):
    override = tmp_path / "override.yaml"
    override.write_text(yaml.safe_dump({"paths": {"cache_dir": "custom-cache"}}))
    config = load_config(override)
    assert Path(config["paths"]["cache_dir"]) == PROJECT_ROOT / "custom-cache"


def test_internal_assets_are_portable_relative_paths():
    raw = yaml.safe_load(DEFAULT_CONFIG.read_text())
    assert raw["tools"]["fastorf"] == "third_party/fastorf/bin/fastorf"
    assert raw["assets"]["phylo_root"] == "resources/phylo"
    assert raw["assets"]["runtime_lib_dir"] == "runtime/libs"
    config = load_config(DEFAULT_CONFIG)
    assert Path(config["tools"]["fastorf"]) == PROJECT_ROOT / "third_party/fastorf/bin/fastorf"
    assert Path(config["assets"]["phylo_root"]) == PROJECT_ROOT / "resources/phylo"
    assert Path(config["assets"]["runtime_lib_dir"]) == PROJECT_ROOT / "runtime/libs"
    for key, value in config["assets"].items():
        if isinstance(value, str) and value.startswith(("resources/", "runtime/", "third_party/")):
            assert Path(value).is_absolute(), key


def test_engineering_threshold_is_explicitly_preserved():
    config = load_config(PROJECT_ROOT / "config" / "engineering_smoke.yaml")
    assert config["min_input_aa"] == 100
    assert config["min_analysis_aa"] == 100
