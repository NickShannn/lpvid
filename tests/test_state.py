import json

import pytest

from cdc_pipeline.state import RunState, atomic_write_json, compute_signature, file_fingerprint


def test_resume_requires_same_signature_and_complete_artifact(tmp_path):
    state = RunState(tmp_path)
    state.start({"route": "proteins"}, {"x": 1})
    artifact = tmp_path / "result.tsv"
    artifact.write_text("ok\n")
    state.begin_stage("detect", "sig-1")
    state.complete_stage("detect", "sig-1", {"table": "result.tsv"})
    assert state.should_skip("detect", "sig-1") is True
    assert state.should_skip("detect", "sig-2") is False
    artifact.unlink()
    assert state.should_skip("detect", "sig-1") is False


def test_state_is_json_and_failure_is_persistent(tmp_path):
    state = RunState(tmp_path)
    state.start({"route": "contigs"}, {"x": 1})
    state.begin_stage("orf", "sig")
    state.fail_stage("orf", "sig", "boom")
    data = json.loads((tmp_path / "state.json").read_text())
    assert data["status"] == "failed"
    assert data["stages"]["orf"]["error"] == "boom"


def test_invalidate_marks_downstream_stale_without_deleting_outputs(tmp_path):
    state = RunState(tmp_path)
    state.start({"route": "proteins"}, {})
    for name in ("detect", "phylo"):
        artifact = tmp_path / f"{name}.tsv"
        artifact.write_text(name)
        state.begin_stage(name, "old")
        state.complete_stage(name, "old", {"table": f"{name}.tsv"})
    state.invalidate_stages(["phylo"], reason="upstream_changed")
    assert not state.should_skip("phylo", "new")
    assert (tmp_path / "phylo.tsv").exists()
    assert state.should_skip("phylo", "old") is True
    assert state.snapshot["stages"]["phylo"]["recovery"] == "checkpoint_restored"


def test_interrupted_running_stage_restores_matching_checkpoint(tmp_path):
    state = RunState(tmp_path)
    state.start({"route": "proteins"}, {})
    artifact = tmp_path / "result.tsv"
    artifact.write_text("ok\n")
    state.begin_stage("detect", "sig-1")
    state.complete_stage("detect", "sig-1", {"table": "result.tsv"}, context_data={"value": 1})
    state.begin_stage("detect", "interrupted-signature")
    assert state.should_skip("detect", "sig-1") is True
    assert state.snapshot["stages"]["detect"]["context_data"] == {"value": 1}


def test_signature_and_large_file_fingerprint_are_stable(tmp_path):
    path = tmp_path / "x.txt"
    path.write_text("abc")
    assert compute_signature({"b": 2, "a": 1}) == compute_signature({"a": 1, "b": 2})
    first = file_fingerprint(path)
    assert first["sha256"] == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    assert file_fingerprint(path) == first


def test_atomic_write_json_recurses_to_dict_and_rejects_nan(tmp_path):
    class Child:
        def __init__(self, value):
            self.value = value

        def to_dict(self):
            return {"value": self.value}

    class Parent:
        def to_dict(self):
            return {"children": [Child(1), Child(None)]}

    output = tmp_path / "summary.json"
    atomic_write_json(output, {"rows": [Parent()]})
    assert json.loads(output.read_text(encoding="utf-8")) == {
        "rows": [{"children": [{"value": 1}, {"value": None}]}]
    }

    with pytest.raises(ValueError):
        atomic_write_json(output, {"invalid": float("nan")})
    assert not output.with_name(output.name + ".tmp").exists()
    assert "rows" in json.loads(output.read_text(encoding="utf-8"))
