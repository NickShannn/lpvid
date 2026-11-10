# lpvid: Quickstart

Version: 0.2

A five-minute walkthrough on the bundled example data. It assumes that
installation is complete and that `lpvid doctor` reports all rows as `[OK]`;
otherwise start with [INSTALL.md](INSTALL.md).

> Chinese version: [QUICKSTART.zh-CN.md](QUICKSTART.zh-CN.md)

The example data lives in [`examples/`](../examples/README.md): nine contigs
(one RdRP-positive and eight negative backgrounds) and simulated paired-end
2x150 reads. It exercises the full workflow but is not a performance benchmark.

## 1. Contigs route

```bash
cd lpvid          # repository root
conda activate lpvid

lpvid run-contigs \
  --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs \
  --sample demo \
  --gpu-id 0
```

The run prints `<status>\t<run_dir>` when finished. Follow progress in another
shell with `lpvid status runs/demo_contigs`.

### Expected results

| Check | Where | Expected |
|---|---|---|
| RdRP-positive ORFs | `results/sequences/rdrp_positive_orfs.faa` | non-empty; all from `demo_rdrp_positive_1` |
| Positive source contig | `results/sequences/rdrp_positive_contigs.fna` | `demo_rdrp_positive_1` only |
| Contig calls | `results/tables/contig_level_calls.tsv` | `demo_rdrp_positive_1` positive; `demo_background_*` negative |
| Supergroup | `results/tables/supergroup_assignments.tsv` | representative of the positive ORF assigned to supergroup Nido by reference match |
| Placement | `results/tables/placements.tsv` | top neighbors inside the Nido backbone; nearest neighbor is an alphacoronavirus RdRP reference |
| Reports | `results/main_results_report_en.html` | the positive contig appears with its supergroup, placement and catalog notices |

Open `results/main_results_report_en.html` and
`results/phylogenetic_explorer.html` in a browser. Both are single offline
files.

## 2. Reads route

```bash
lpvid run-reads \
  --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads \
  --sample demo \
  --gpu-id 0
```

This route adds quality control, rRNA depletion, assembly and abundance
estimation, then proceeds through the same identification, classification and
placement steps.

Expected results:

- `results/tables/abundance.tsv`: the re-assembled counterpart of
  `demo_rdrp_positive_1`, with mapped reads, RPM, mean depth and coverage
  breadth. The example reads give roughly 25x coverage on the positive genome;
  identifiers come from the assembly and differ from the input contig names.
- `results/tables/supergroup_abundance.tsv`: one row for the detected
  supergroup.
- Contigs from `demo_background_*` may assemble partially; they are expected to
  stay RdRP-negative.

## 3. Proteins route (optional)

The proteins route expects an amino-acid FASTA, with each sequence at least
150 aa for the main analysis set. Run it with your own protein input:

```bash
lpvid run-proteins \
  --input proteins.faa \
  --out runs/protein_demo --sample demo --gpu-id 0
```

Protein input has no source nucleic acid, so contig-level tables stay empty by
design.

## 4. Resuming and inspecting

```bash
lpvid status runs/demo_contigs --json   # per-stage status
lpvid resume runs/demo_contigs          # re-run stale stages only
```

## 5. Scope of the example data

- The example inputs verify end-to-end connectivity, configuration and reporting
  on a small known-positive input. They are not a sensitivity, specificity or
  runtime benchmark; use real evaluation cohorts for that.
- Simulated reads carry about 0.8% substitution errors and 300 bp mean
  fragments; assembly quality is intentionally modest.
- Background contigs are non-viral assembly products. Occasional
  `KEEP_unknown` screening calls for them are expected.

Next steps: [RUNBOOK.md](RUNBOOK.md) for production routes and defaults,
[RESULTS.md](RESULTS.md) for table definitions.
