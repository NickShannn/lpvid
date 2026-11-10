"""外部命令执行、环境隔离、日志和 GPU 选择。"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

from cdc_pipeline.state import utc_now

ENV_NAMES = frozenset({"lucaphylo", "qc", "rdrp"})


_ENV_BIN_CACHE: dict[str, str] = {}


class CommandError(RuntimeError):
    pass


def env_bin(env_name: str) -> str:
    if env_name not in ENV_NAMES:
        raise KeyError(f"未知运行环境: {env_name}")
    if env_name not in _ENV_BIN_CACHE:
        override = os.environ.get(f"LPVID_ENV_BIN_{env_name.upper().replace('-', '_')}")
        _ENV_BIN_CACHE[env_name] = override or str(Path(sys.executable).resolve().parent)
    return _ENV_BIN_CACHE[env_name]


def build_env(env_name: str, extra: dict[str, str] | None = None) -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = env_bin(env_name) + os.pathsep + env.get("PATH", "")
    env["PYTHONNOUSERSITE"] = "1"
    env["TOKENIZERS_PARALLELISM"] = "false"
    if extra:
        env.update({key: str(value) for key, value in extra.items()})
    return env


def run_command(
    argv: Iterable[str | Path],
    *,
    env_name: str,
    cwd: str | Path,
    log_path: str | Path,
    env: dict[str, str] | None = None,
    echo: bool = False,
    allowed_returncodes: tuple[int, ...] = (0,),
    stdout_path: str | Path | None = None,
    discard_stdout: bool = False,
) -> None:
    argv = [str(part) for part in argv]
    cwd = Path(cwd)
    log_path = Path(log_path)
    cwd.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if stdout_path:
        stdout_path = Path(stdout_path)
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    commands_path = log_path.parent / "commands.jsonl"
    started = time.monotonic()
    record = {
        "argv": argv,
        "command": shlex.join(argv),
        "cwd": str(cwd),
        "env_name": env_name,
        "started_at": utc_now(),
        "stdout_path": str(stdout_path) if stdout_path else None,
    }
    returncode = 127
    try:
        with log_path.open("a", encoding="utf-8") as log:
            log.write("\n$ " + shlex.join(argv) + "\n")
            log.flush()
            stdout_handle = stdout_path.open("w", encoding="utf-8") if stdout_path else None
            try:
                process = subprocess.Popen(
                    argv,
                    cwd=cwd,
                    env=build_env(env_name, env),
                    stdout=subprocess.DEVNULL if discard_stdout else (stdout_handle or subprocess.PIPE),
                    stderr=subprocess.PIPE if (stdout_handle or discard_stdout) else subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
            except OSError as exc:
                log.write(f"启动失败: {exc}\n")
                raise CommandError(f"命令启动失败: {argv[0]}: {exc}") from exc
            try:
                stream = process.stderr if (stdout_handle or discard_stdout) else process.stdout
                assert stream is not None
                for line in stream:
                    log.write(line)
                    if echo:
                        print(line, end="", flush=True)
                returncode = process.wait()
            finally:
                if stdout_handle:
                    stdout_handle.close()
    finally:
        record.update({
            "returncode": returncode,
            "ended_at": utc_now(),
            "duration_seconds": round(time.monotonic() - started, 3),
            "log": str(log_path),
        })
        with commands_path.open("a", encoding="utf-8") as commands:
            commands.write(json.dumps(record, ensure_ascii=False) + "\n")
    if returncode not in allowed_returncodes:
        raise CommandError(f"命令 exit {returncode}: {shlex.join(argv)}；日志: {log_path}")


def pick_free_gpu(requested: str | int = "auto") -> str:
    if str(requested).lower() != "auto":
        return str(requested)
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,memory.free", "--format=csv,noheader,nounits"],
            check=True,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError("自动选卡失败；请用 --gpu-id 显式指定 GPU") from exc
    candidates = []
    for line in result.stdout.splitlines():
        index, free = [part.strip() for part in line.split(",")]
        candidates.append((int(free), int(index)))
    if not candidates:
        raise RuntimeError("nvidia-smi 未返回可用 GPU")
    return str(max(candidates)[1])
