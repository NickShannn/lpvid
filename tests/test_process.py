import json
import sys

import pytest

from cdc_pipeline.process import CommandError, run_command


def test_run_command_writes_log_and_command_record(tmp_path):
    log = tmp_path / "run.log"
    run_command([sys.executable, "-c", "print('hello')"], env_name="lucaphylo", cwd=tmp_path, log_path=log)
    assert "hello" in log.read_text()
    record = json.loads((tmp_path / "commands.jsonl").read_text().splitlines()[0])
    assert record["returncode"] == 0
    assert record["argv"][0] == sys.executable


def test_run_command_can_capture_stdout_without_shell(tmp_path):
    output = tmp_path / "stats.tsv"
    run_command(
        [sys.executable, "-c", "print('a\\tb')"],
        env_name="lucaphylo", cwd=tmp_path, log_path=tmp_path / "capture.log", stdout_path=output,
    )
    assert output.read_text().strip() == "a\tb"
    record = json.loads((tmp_path / "commands.jsonl").read_text().splitlines()[0])
    assert record["stdout_path"] == str(output)


def test_run_command_fails_loudly(tmp_path):
    log = tmp_path / "run.log"
    with pytest.raises(CommandError, match="exit 7"):
        run_command([sys.executable, "-c", "raise SystemExit(7)"], env_name="lucaphylo", cwd=tmp_path, log_path=log)
    assert json.loads((tmp_path / "commands.jsonl").read_text().splitlines()[0])["returncode"] == 7
