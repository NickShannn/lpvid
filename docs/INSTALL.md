# lpvid: Installation

Version: 0.2

This guide takes a fresh Linux machine from zero to a validated lpvid
installation. Source code and large resources (model weights, databases) are
installed separately.

> Chinese version: [INSTALL.zh-CN.md](INSTALL.zh-CN.md)

## 1. System requirements

| Component | Requirement |
|---|---|
| OS | Linux x86_64 (glibc-based; Ubuntu 20.04+ and CentOS 8+ tested) |
| CPU | 8 or more cores recommended; assembly and database tools are multi-threaded |
| RAM | 64 GB recommended; 32 GB minimum for small samples |
| GPU | NVIDIA GPU with at least 16 GB VRAM; 24 GB recommended |
| NVIDIA driver | 550.x or newer; verify with `nvidia-smi` |
| CUDA runtime | 12.8, provided by the PyTorch wheel; no local CUDA toolkit is required |
| Disk | Repository and environment about 15 GB, plus resources |
| Software | conda or mamba (Miniconda or Miniforge), git, pip |

Resources require substantial additional disk space, dominated by the protein
reference database. The exact file inventory and sizes are recorded in
[resources/RESOURCE_MANIFEST.tsv](../resources/RESOURCE_MANIFEST.tsv) for the
release.

Node.js is needed only for the optional tree-explorer JavaScript tests, not for
running the workflow.

## 2. Clone the repository

```bash
git clone https://github.com/NickShannn/lpvid.git
cd lpvid
git status        # branch or tag: 0.2
```

## 3. Create the conda environment

The unified environment `lpvid` contains Python 3.10, the required
bioinformatics tools and the Python dependencies, including PyTorch 2.8.0
(CUDA 12.8).

```bash
# conda
conda env create -f environment.yml -n lpvid
conda activate lpvid

# or mamba (faster solver)
mamba env create -f environment.yml -n lpvid
conda activate lpvid
```

The environment includes `fastp`, `bowtie2`, `MEGAHIT`, `samtools`, `DIAMOND`,
`CD-HIT` and `seqkit`. The `fastorf` binary is bundled under
`third_party/fastorf/bin/`.

`environment.yml` already defines `name: lpvid`; `-n lpvid` makes it explicit.
All direct conda and Python dependencies are pinned to exact versions. A fully
pinned transitive lock file is kept in `config/requirements-lpvid.txt`; install
it when bit-for-bit reproduction is required:

```bash
pip install -r config/requirements-lpvid.txt
```

Verify the GPU stack and a core tool:

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# expect: 2.8.0 True
diamond --version
```

If `torch.cuda.is_available()` is `False`, check `nvidia-smi` and update the
NVIDIA driver. If your package mirror serves CPU-only wheels, install PyTorch
from the CUDA 12.8 index:

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128
```

## 4. Install the workflow

```bash
pip install --no-deps -e .
lpvid --version      # 0.2
lpvid --help
```

`pip install -e .` registers the `lpvid` console script while keeping the source
tree editable. The repository-root `./lpvid` launcher also works. Do not install
models or databases into site-packages; they are referenced through
configuration.

## 5. Place resources

Small runtime resources are already bundled in the git repository: the RdRP
reference database, backbone trees and taxonomy tables, the classifier launcher,
labels and runtime model files, and the two RefSeq taxid whitelists. Only the
large model checkpoints, the 24 LoRA placement adapters, the full RefSeq protein
DIAMOND database and the SILVA rRNA indexes must be obtained separately. Place
them under `resources/` following the layout below. The authoritative inventory
is [resources/RESOURCE_MANIFEST.tsv](../resources/RESOURCE_MANIFEST.tsv); the
downloadable archives, their contents and SHA-256 values are listed in
[resources/EXTERNAL_RESOURCES.md](../resources/EXTERNAL_RESOURCES.md).

```text
resources/
  phylo/
    backbone_dir/<supergroup>/{*.nwk,*.fas}   # bundled
    adapter_dir/<supergroup>/...              # Zenodo resource record
    taxonomy/<supergroup>.csv                  # bundled
  rdrp/
    rdrp_reference.dmnd                        # bundled
    rdrp_reference.fas                         # bundled
    rdrp_reference.csv                         # bundled
  models/
    rdrp_detector/                             # Zenodo resource record
    lucaphylo/                                 # Zenodo resource record
    classifier/
      pipeline/                                # bundled
      encoder/                                 # Zenodo resource record
      model/                                   # bundled
  blastdb/refseq_protein/
    refseq_protein.dmnd                        # separate distribution (large)
    TAXID_viruses.txt                          # bundled
    TAXID_cellular.txt                         # bundled
  rrna/
    silva_ssu.*                                # Zenodo resource record
    silva_lsu.*                                # Zenodo resource record
```

Verify the downloaded archives against the SHA-256 values in
`resources/EXTERNAL_RESOURCES.md`, then verify files against the manifest.
Bundled files carry SHA-256 values; files inside the archives are marked
`verify_after_download` in the manifest.

## 6. Configure local paths

```bash
cp config/local.example.yaml config/local.yaml
```

`config/local.yaml` is machine-private and git-ignored. The template assumes the
unified environment created above; if conda lives elsewhere, edit the
`/opt/conda/envs/lpvid` prefix (`conda info --base` shows the root). It also
assumes the `resources/` layout from step 5. Override scientific or runtime
defaults in the same file rather than editing `config/default.yaml`.

Run `lpvid` from the repository root so that relative paths in the template
resolve; alternatively put absolute paths in `config/local.yaml`.

If the precompiled `third_party/fastorf/bin/fastorf` does not run on your
platform, build it with Rust and refresh the binary:

```bash
cd third_party/fastorf
cargo build --release
mkdir -p bin && cp target/release/fastorf bin/fastorf
cd ../..
```

## 7. Validate the deployment

```bash
lpvid doctor      # environments, tools, models and databases; every row must be [OK]
lpvid prepare     # requires backbone trees, taxonomy and adapters; builds local reference assets
./scripts/validate.sh   # unit tests and tree-explorer core tests
```

`prepare` is deterministic: identical inputs are not rebuilt. It does not
download anything from the internet.

## 8. First runs

```bash
# contigs, bundled example data
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo --gpu-id 0

# paired-end reads, bundled example data
lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo --gpu-id 0
```

Expected results and a five-minute walkthrough:
[QUICKSTART.md](QUICKSTART.md).

### Understanding the output

Each run directory contains:

```text
run_manifest.json              inputs, merged configuration, resource fingerprints
state.json / checkpoint.json    per-stage status and signatures
step*/                         intermediate stage outputs and logs
results/                       delivery: HTML reports, tables/, sequences/, trees/
```

Key files in `results/`:

| Path | Content |
|---|---|
| `main_results_report_zh.html` / `main_results_report_en.html` | Self-contained bilingual main reports |
| `phylogenetic_explorer.html` | Interactive placement tree with SVG, PNG, PDF and Newick export |
| `tables/orf_predictions.tsv` | Per-protein RdRP probability and detection call |
| `tables/contig_level_calls.tsv` | Per-contig calls; any positive ORF makes the contig positive |
| `tables/supergroup_assignments.tsv` | 180-supergroup assignment with reference-match or classifier-classification evidence |
| `tables/placements.tsv` | LucaPhylo top-5 backbone neighbors and distances |
| `tables/abundance.tsv` | reads route only: mapped reads, RPM, depth and coverage |
| `sequences/rdrp_positive_*.fna` / `.faa` | RdRP-positive source contigs and ORFs |
| `results/README.md` | Generated bilingual index of delivered files |

Column-by-column reference: [RESULTS.md](RESULTS.md).

## 9. Updating

The distribution is versioned by tag; the current release is `0.2`.

```bash
git fetch origin
git pull origin main        # mainline updates
# or switch to a newer tag: git checkout <tag>
conda env update -f environment.yml -n lpvid
pip install --no-deps -e .
lpvid doctor
```

Bundled resources are updated with the source repository. Downloadable model,
adapter and index resources are published as a versioned resource record;
update them according to the release notes of the matching version.

## 10. Troubleshooting

| Symptom | Fix |
|---|---|
| `lpvid: command not found` | Activate the `lpvid` environment, re-run `pip install --no-deps -e .`, or use `./lpvid` from the repository root |
| `doctor` reports `missing` | Correct the matching entry in `config/local.yaml`; do not disable the check |
| `torch.cuda.is_available()` is `False` | Update the NVIDIA driver to 550 or newer; reinstall PyTorch from the CUDA 12.8 index |
| `fastp`, `bowtie2`, `MEGAHIT` or `samtools` not found | They live in the `lpvid` environment; ensure it is active and the environment paths in `config/local.yaml` point into it |
| Reference database jobs are slow | Large database scans are I/O-bound; a small query can still take about 10 minutes on the full database. Do not run several full-library jobs concurrently |
| `fastorf` binary fails to execute | Rebuild it from source as shown in step 6, and check the glibc version |
| GPU out of memory | Reduce the model batch sizes in `config/local.yaml`, or select a larger GPU with `--gpu-id` |
| A run fails mid-way | `lpvid status <run_dir> --json`; fix the cause; `lpvid resume <run_dir>` re-runs stale stages only |

## 11. Uninstall

```bash
conda deactivate
conda env remove -n lpvid
pip uninstall lpvid
```

Runs, resources and caches are independent of the installation and remain where
you placed them.

## 12. License and citation

lpvid is released under the Apache License 2.0. The license text is in the
repository root as `LICENSE`; bundled components and their license texts are
listed in `THIRD_PARTY_NOTICES.md`. Model weights, adapters, reference data and
databases are data assets with terms stated on the accompanying resource
record.

Citation metadata is provided in `CITATION.cff`.
