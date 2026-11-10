#!/usr/bin/env python3
"""Build examples/contigs/demo_contigs.fa from one validation run.

Input is a contigs-route run directory of lpvid. The script takes
the first RdRP-positive contig, adds a fixed number of negative background
contigs, renames every record to demo_* identifiers, and writes a FASTA file
small enough to ship with the repository.

Usage:
    python3 extract_demo_contigs.py \
        --positive-fna RUN/03_detect/rdrp_positive_contigs.fna \
        --calls-tsv RUN/03_detect/contig_level_calls.tsv \
        --source-contigs SAMPLE.contigs.fa \
        --out examples/contigs/demo_contigs.fa

Only the first positive contig is used. Background contigs are the first N
records called both_neg in the detection table whose length falls in
[min_len, max_len].
"""

import argparse
from pathlib import Path


def read_fasta(path):
    record_id = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(chunks)
                record_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if record_id is not None:
        yield record_id, "".join(chunks)


def negative_contig_ids(calls_tsv, limit, min_len_hint):

    ids = []
    with open(calls_tsv) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col_id = header.index("contig_id")
        col_call = header.index("call")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[col_call] == "both_neg":
                ids.append(fields[col_id])
                if len(ids) >= limit * 20:
                    break
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positive-fna", required=True)
    parser.add_argument("--calls-tsv", required=True)
    parser.add_argument("--source-contigs", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--background-count", type=int, default=8)
    parser.add_argument("--min-len", type=int, default=700)
    parser.add_argument("--max-len", type=int, default=2500)
    args = parser.parse_args()

    positives = list(read_fasta(args.positive_fna))
    if not positives:
        raise SystemExit("no RdRP-positive contigs in input")
    positive_id, positive_seq = positives[0]
    if len(positive_seq) > 40000:
        raise SystemExit("positive contig longer than 40 kb, pick a smaller one")

    wanted = set(negative_contig_ids(args.calls_tsv, args.background_count, args.min_len))
    background = []
    for record_id, sequence in read_fasta(args.source_contigs):
        if record_id == positive_id:
            continue
        if record_id in wanted and args.min_len <= len(sequence) <= args.max_len:
            background.append((record_id, sequence))
            if len(background) >= args.background_count:
                break
    if len(background) < args.background_count:
        raise SystemExit(
            f"only {len(background)} negative contigs in [{args.min_len}, {args.max_len}]"
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        handle.write(f">demo_rdrp_positive_1 len={len(positive_seq)}\n")
        for i in range(0, len(positive_seq), 60):
            handle.write(positive_seq[i : i + 60] + "\n")
        for index, (_, sequence) in enumerate(background, start=1):
            handle.write(f">demo_background_{index} len={len(sequence)}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")

    total = len(positive_seq) + sum(len(s) for _, s in background)
    print(f"wrote {1 + len(background)} records, {total} bp to {out_path}")


if __name__ == "__main__":
    main()
