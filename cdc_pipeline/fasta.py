"""轻量 FASTA 读写和跨阶段稳定 ID。"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Iterator, TextIO


@dataclass(frozen=True)
class FastaRecord:
    record_id: str
    sequence: str
    description: str = ""
    query_id: str = ""

    @property
    def header_id(self) -> str:
        return self.query_id or self.record_id

    @property
    def length(self) -> int:
        return len(self.sequence)


def _open(path: str | Path) -> TextIO:
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    """逐条读取 FASTA；record_id 取首 token，description 保留完整原始标题。"""
    header = None
    chunks: list[str] = []
    with _open(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _make_record(header, chunks)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: FASTA 标题为空")
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: 序列出现在首个 FASTA 标题之前")
                chunks.append("".join(line.split()).upper())
        if header is not None:
            yield _make_record(header, chunks)


def _make_record(header: str, chunks: list[str]) -> FastaRecord:
    return FastaRecord(
        record_id=header.split()[0],
        sequence="".join(chunks),
        description=header,
    )


def write_fasta(path: str | Path, records: Iterable[FastaRecord], line_width: int = 80) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            header_id = record.header_id
            description = record.description.strip()
            if not description:
                header = header_id
            elif description == header_id or description.startswith(header_id + " "):
                header = description
            else:
                header = f"{header_id} {description}"
            handle.write(f">{header}\n")
            for start in range(0, len(record.sequence), line_width):
                handle.write(record.sequence[start:start + line_width] + "\n")
            if not record.sequence:
                handle.write("\n")
    return path


def stable_id(sample: str, number: int) -> str:
    """运行内稳定 query ID；样本名始终作为独立联合键保存。"""
    if not sample:
        raise ValueError("sample 不能为空")
    if number < 1:
        raise ValueError("FASTA ID 序号从 1 开始")
    return f"Q{number:06d}"


def assign_stable_ids(records: Iterable[FastaRecord], sample: str) -> tuple[list[FastaRecord], list[dict[str, str]]]:
    assigned = []
    mapping = []
    for number, record in enumerate(records, 1):
        query_id = stable_id(sample, number)
        assigned.append(replace(record, query_id=query_id))
        mapping.append({
            "sample": sample,
            "query_id": query_id,
            "original_id": record.record_id,
            "original_description": record.description,
        })
    return assigned, mapping


def write_id_map(path: str | Path, mapping: Iterable[dict[str, str]]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "query_id", "original_id", "original_description"]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in mapping:
            handle.write("\t".join(str(row.get(field, "")).replace("\t", " ") for field in fields) + "\n")
    return path
