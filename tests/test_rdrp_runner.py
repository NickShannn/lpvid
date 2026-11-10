from pathlib import Path

from cdc_pipeline.fasta import FastaRecord
from cdc_pipeline.stages.rdrp_runner import prepare_model_bridge, parse_rdrp_text, screen_rdrp_input


def test_toxic_sequence_only_excludes_rdrp_input():
    toxic = FastaRecord("toxic", "ACGTN" * 20)
    starts_m = FastaRecord("starts_m", "M" + "A" * 99)
    normal = FastaRecord("normal", "M" + "ACDEFGHIKLMNPQRSTVWY" * 5)
    clean, excluded = screen_rdrp_input([toxic, starts_m, normal])
    assert [r.record_id for r in clean] == ["starts_m", "normal"]
    assert [r.record_id for r in excluded] == ["toxic"]


def test_rdrp_parser_uses_named_columns():
    rows = parse_rdrp_text("seq_id,seq,prob,label_index,label\nQ1,ACD,0.8,1,positive\n")
    assert rows["Q1"].prob == 0.8
    assert rows["Q1"].label == "positive"


def test_model_bridge_links_resolve_when_repo_path_is_relative(tmp_path, monkeypatch):
    actual = tmp_path / "actual-repo"
    checkpoint = actual / "models/RdRP/protein/binary_class_v2/luca_base/matrix/protv2full_20261009100409/checkpoint-95936"
    checkpoint.mkdir(parents=True)
    (checkpoint / "training_args.bin").write_bytes(b"args")
    (checkpoint.parent / "best/label.txt").parent.mkdir(parents=True)
    (checkpoint.parent / "best/label.txt").write_text("label\n")
    (actual / "requirements.txt").write_text("dependencies\n")
    resources = tmp_path / "resources"
    resources.mkdir()
    linked_repo = resources / "rdrp_detector"
    linked_repo.symlink_to(actual)
    bridge = tmp_path / "run/model_bridge"
    monkeypatch.chdir(tmp_path)
    prepare_model_bridge(Path("resources/rdrp_detector"), bridge)
    assert (bridge / "requirements.txt").resolve() == actual / "requirements.txt"
    assert (bridge / "models/RdRP/protein/binary_class/luca_base/matrix/protv2full_20261009100409/checkpoint-95936/training_args.bin").is_file()
