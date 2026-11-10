from cdc_pipeline.config import PROJECT_ROOT, load_config
from cdc_pipeline.project import project_metadata


def test_project_identity_and_single_command():
    metadata = project_metadata()
    assert metadata["display_name"] == "lpvid"
    assert metadata["distribution"] == "lpvid"
    assert metadata["command"] == "lpvid"
    assert metadata["version"] == "0.2"
    assert (PROJECT_ROOT / "lpvid").is_file()
    assert not (PROJECT_ROOT / "cdc-rdrp").exists()


def test_current_user_entrypoints_do_not_reference_legacy_command():
    for relative in ["README.md", "docs/RUNBOOK.md", "scripts/validate.sh", "setup.py"]:
        assert "cdc-rdrp" not in (PROJECT_ROOT / relative).read_text(), relative


def test_config_exposes_same_project_identity():
    config = load_config()
    assert config["project"] == project_metadata()
    assert config["project"]["display_name"] == "lpvid"
