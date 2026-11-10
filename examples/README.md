# Bundled example data

Version: 0.2

Small test inputs for verifying an lpvid installation. They exercise the full
workflow on a known positive and are intended for connectivity checks only, not
for sensitivity, specificity or runtime evaluation.

> Chinese version: [README.zh-CN.md](README.zh-CN.md)

## Files

| File | Content |
|---|---|
| `contigs/demo_contigs.fa` | 9 contigs, 39,577 bp total: `demo_rdrp_positive_1` (27,452 bp, RdRP-positive) and `demo_background_1` to `demo_background_8` (703 to 2,238 bp, non-viral) |
| `reads/demo_R1.fastq.gz` | 2,612 paired-end reads, 2x150, about 114 KB gzipped |
| `reads/demo_R2.fastq.gz` | mate of R1 |
| `src/extract_demo_contigs.py` | script that produced the contigs file |
| `src/simulate_demo_reads.py` | read simulator (fragment N(300, 30), 0.8% substitution errors, Illumina-like qualities, seed 20261110) |

The positive contig is derived from a real metatranscriptomic assembly validated
in an lpvid run. It carries an RdRP-positive ORF, is assigned to the Nido
supergroup, and is placed near alphacoronavirus references. Background contigs
are negative, non-viral contigs from the same assembly. Record identifiers were
renamed; no sample identifiers are retained.

## Usage

```bash
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo --gpu-id 0

lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo --gpu-id 0
```

## Expected outcome

- `demo_rdrp_positive_1`: at least one ORF called RdRP-positive; the contig
  appears in `results/sequences/rdrp_positive_contigs.fna`; its representative
  is assigned to supergroup Nido by reference match and placed next to
  alphacoronavirus backbone references in `results/tables/placements.tsv`.
- `demo_background_*`: no positive calls expected.
- reads route: the positive genome re-assembles with about 25x coverage and gets
  RPM, depth and coverage rows in `results/tables/abundance.tsv`.

See [docs/QUICKSTART.md](../docs/QUICKSTART.md) for the full walkthrough.

## Regenerating

Both scripts use only the Python standard library:

```bash
python3 src/extract_demo_contigs.py \
    --positive-fna <validation_run>/sequences/rdrp_positive_contigs.fna \
    --calls-tsv <validation_run>/tables/contig_level_calls.tsv \
    --source-contigs <assembly>.contigs.fa --out contigs/demo_contigs.fa

python3 src/simulate_demo_reads.py --contigs contigs/demo_contigs.fa \
    --r1 reads/demo_R1.fastq.gz --r2 reads/demo_R2.fastq.gz
```

`extract_demo_contigs.py` requires access to the original validation run and
assembly. It is shipped for provenance and auditability, not for redelivery of
the source data.
