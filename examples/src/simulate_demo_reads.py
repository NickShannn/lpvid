#!/usr/bin/env python3
"""Simulate paired-end Illumina reads (2x150) from demo contigs.

Fragments are drawn along each contig with a normal(300, 30) length
distribution clipped to [230, 450]; R1 is the forward prefix, R2 the reverse
complement of the suffix, so short fragments overlap like real PE data.
Base errors are substitutions at a fixed rate with Illumina-like quality
strings that degrade towards the read tail.

Usage:
    python3 simulate_demo_reads.py \
        --contigs examples/contigs/demo_contigs.fa \
        --r1 examples/reads/demo_R1.fastq.gz \
        --r2 examples/reads/demo_R2.fastq.gz
"""

import argparse
import gzip
import random
from pathlib import Path

READ_LEN = 150
FRAG_MEAN = 300
FRAG_SD = 30
FRAG_MIN = 230
FRAG_MAX = 450
ERROR_RATE = 0.008

QUAL_HEAD = "I" * 80
QUAL_TAIL = "".join(chr(33 + max(20, 40 - int((i - 80) * 20 / 70))) for i in range(80, READ_LEN))
QUAL_TEMPLATE = QUAL_HEAD + QUAL_TAIL
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


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


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def simulate_read_pair(rng, contig, start, frag_len):
    r1 = list(contig[start : start + READ_LEN])
    r2 = list(revcomp(contig[start + frag_len - READ_LEN : start + frag_len]))
    quals = []
    for read in (r1, r2):
        qual = list(QUAL_TEMPLATE)
        for i, base in enumerate(read):
            if rng.random() < ERROR_RATE:
                read[i] = rng.choice([b for b in "ACGT" if b != read[i]])
                qual[i] = chr(33 + max(14, ord(qual[i]) - 33 - 12))
        quals.append("".join(qual))
    return "".join(r1), quals[0], "".join(r2), quals[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--positive-id", default="demo_rdrp_positive_1")
    parser.add_argument("--positive-coverage", type=float, default=25.0)
    parser.add_argument("--background-coverage", type=float, default=8.0)
    parser.add_argument("--seed", type=int, default=20261110)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    contigs = list(read_fasta(args.contigs))

    pairs = []
    for record_id, sequence in contigs:
        coverage = args.positive_coverage if record_id == args.positive_id else args.background_coverage
        n_pairs = int(round(coverage * len(sequence) / FRAG_MEAN))
        for _ in range(n_pairs):
            frag_len = max(FRAG_MIN, min(FRAG_MAX, round(rng.gauss(FRAG_MEAN, FRAG_SD))))
            if frag_len > len(sequence):
                frag_len = len(sequence)
            start = rng.randrange(0, len(sequence) - frag_len + 1)
            pairs.append(simulate_read_pair(rng, sequence, start, frag_len))
    rng.shuffle(pairs)

    for path, mate in ((args.r1, 0), (args.r2, 1)):
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out, "wt", compresslevel=9) as handle:
            for index, pair in enumerate(pairs, start=1):
                seq, qual = pair[mate * 2], pair[mate * 2 + 1]
                mate_tag = "1:N:0:1" if mate == 0 else "2:N:0:1"
                handle.write(f"@demo:{index}:{mate_tag}\n{seq}\n+\n{qual}\n")
        print(f"wrote {len(pairs)} reads to {out}")


if __name__ == "__main__":
    main()
