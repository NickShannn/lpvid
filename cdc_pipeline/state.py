"""运行清单、原子状态写入、阶段签名和断点恢复。"""

from __future__ import annotations

import copy
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

SMALL_FILE_HASH_LIMIT = 64 * 1024 * 1024
TERMINAL_SKIP_STATUSES = {"done", "skipped", "not_applicable"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_default(value: Any) -> Any:
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return to_dict()
    return str(value)


def atomic_write_json(path: str | Path, payload: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(
                payload, handle, ensure_ascii=False, indent=2, sort_keys=True,
                default=_json_default, allow_nan=False,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return path


def file_fingerprint(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    stat = path.stat()
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if path.is_file() and stat.st_size <= SMALL_FILE_HASH_LIMIT:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    return result


def compute_signature(payload: Any) -> str:
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _artifact_paths(value: Any) -> Iterable[str]:
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from _artifact_paths(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _artifact_paths(child)


class RunState:
    def __init__(self, run_dir: str | Path):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.state_path = self.run_dir / "state.json"
        self.checkpoint_path = self.run_dir / "checkpoint.json"
        self.manifest_path = self.run_dir / "run_manifest.json"
        self._data: dict[str, Any] | None = None
        if self.state_path.exists():
            self._data = json.loads(self.state_path.read_text(encoding="utf-8"))

    @property
    def manifest(self) -> dict[str, Any]:
        if not self.manifest_path.exists():
            raise FileNotFoundError(f"运行清单不存在: {self.manifest_path}")
        return json.loads(self.manifest_path.read_text(encoding="utf-8"))

    @property
    def snapshot(self) -> dict[str, Any]:
        return copy.deepcopy(self._data or {"status": "pending", "stages": {}})

    def start(self, metadata: dict[str, Any], config: dict[str, Any]) -> None:
        if self.state_path.exists():
            raise FileExistsError(f"运行已存在；请使用 resume: {self.run_dir}")
        created = utc_now()
        run_id = self.run_dir.name
        atomic_write_json(self.manifest_path, {
            "run_id": run_id,
            "created_at": created,
            **metadata,
            "config": config,
        })
        self._data = {
            "run_id": run_id,
            "status": "running",
            "created_at": created,
            "updated_at": created,
            "stages": {},
        }
        self._flush()
        atomic_write_json(self.checkpoint_path, {"run_id": run_id, "stages": {}})

    def _require_started(self) -> dict[str, Any]:
        if self._data is None:
            raise RuntimeError(f"运行尚未初始化: {self.run_dir}")
        return self._data

    def _flush(self) -> None:
        data = self._require_started()
        data["updated_at"] = utc_now()
        atomic_write_json(self.state_path, data)

    def begin_stage(self, name: str, signature: str) -> None:
        data = self._require_started()
        previous = data["stages"].get(name, {})
        data["stages"][name] = {
            **previous,
            "status": "running",
            "signature": signature,
            "started_at": utc_now(),
            "ended_at": None,
            "error": None,
        }
        data["status"] = "running"
        self._flush()

    def complete_stage(
        self,
        name: str,
        signature: str,
        artifacts: dict[str, Any],
        status: str = "done",
        summary: dict[str, Any] | None = None,
        context_data: dict[str, Any] | None = None,
    ) -> None:
        if status not in TERMINAL_SKIP_STATUSES:
            raise ValueError(f"非法完成状态: {status}")
        data = self._require_started()
        data["stages"][name] = {
            **data["stages"].get(name, {}),
            "status": status,
            "signature": signature,
            "artifacts": artifacts,
            "summary": summary or {},
            "context_data": context_data or {},
            "ended_at": utc_now(),
            "error": None,
        }
        self._flush()
        checkpoint = json.loads(self.checkpoint_path.read_text(encoding="utf-8")) if self.checkpoint_path.exists() else {"stages": {}}
        checkpoint["stages"][name] = {
            "signature": signature,
            "status": status,
            "artifacts": artifacts,
            "updated_at": utc_now(),
        }
        atomic_write_json(self.checkpoint_path, checkpoint)

    def fail_stage(self, name: str, signature: str, error: str) -> None:
        data = self._require_started()
        data["stages"][name] = {
            **data["stages"].get(name, {}),
            "status": "failed",
            "signature": signature,
            "error": error,
            "ended_at": utc_now(),
        }
        data["status"] = "failed"
        self._flush()

    def should_skip(self, name: str, signature: str) -> bool:
        data = self._require_started()
        stage = data["stages"].get(name)
        if not stage:
            return False
        if stage.get("status") in TERMINAL_SKIP_STATUSES and stage.get("signature") == signature \
                and self._artifacts_complete(stage.get("artifacts", {})):
            return True

        if stage.get("status") in {"running", "stale"} and self.checkpoint_path.exists():
            checkpoint = json.loads(self.checkpoint_path.read_text(encoding="utf-8"))
            saved = checkpoint.get("stages", {}).get(name)
            if saved and saved.get("status") in TERMINAL_SKIP_STATUSES \
                    and saved.get("signature") == signature \
                    and self._artifacts_complete(saved.get("artifacts", {})):
                stage.update({
                    "status": saved["status"],
                    "signature": signature,
                    "artifacts": saved.get("artifacts", {}),
                    "error": None,
                    "recovery": "checkpoint_restored",
                })
                self._flush()
                return True
        return False

    def _artifacts_complete(self, artifacts: dict) -> bool:
        for artifact in _artifact_paths(artifacts):
            path = Path(artifact)
            if not path.is_absolute():
                path = self.run_dir / path
            if not path.exists():
                return False
        return True

    def invalidate_stages(self, names: Iterable[str], reason: str) -> None:
        data = self._require_started()
        changed = False
        for name in names:
            if name in data["stages"]:
                data["stages"][name]["status"] = "stale"
                data["stages"][name]["stale_reason"] = reason
                changed = True
        if changed:
            self._flush()

    def finish(self, status: str = "done") -> None:
        data = self._require_started()
        data["status"] = status
        data["ended_at"] = utc_now()
        self._flush()
