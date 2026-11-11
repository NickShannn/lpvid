# lpvid

Version: 0.2

lpvid is a command-line workflow for RNA virus discovery and characterization
from metatranscriptomic data. It accepts sequencing reads, assembled contigs or
protein sequences, identifies RdRP-containing sequences, classifies them in the
180-supergroup scheme, places representative sequences on curated phylogenetic
backbone trees, estimates abundance for reads input, and writes self-contained
bilingual HTML reports.

## Capabilities

- Three input routes: paired-end or single-end reads, assembled nucleotide
  contigs, and amino-acid FASTA.
- RdRP identification with per-protein and per-contig calls.
- Supergroup classification in the 180-supergroup scheme.
- Phylogenetic placement with LucaPhylo on curated backbone trees, plus a
  standalone interactive tree explorer.
- Abundance analysis for reads input: mapped reads, RPM, mean depth and
  coverage breadth.
- Bilingual (English and Chinese) offline HTML reports with sortable tables,
  catalog-based notices, and SVG, PNG, PDF and Newick export.
- Resumable runs: each stage is fingerprinted, and `lpvid resume` re-executes
  only the stages that changed.

## Inputs

| Route | Input | Command |
|---|---|---|
| reads | Paired-end or single-end FASTQ | `lpvid run-reads` |
| contigs | Assembled nucleotide FASTA | `lpvid run-contigs` |
| proteins | Amino-acid FASTA | `lpvid run-proteins` |

## Installation summary

```bash
git clone https://github.com/<org>/lpvid.git
cd lpvid

conda env create -f environment.yml -n lpvid
conda activate lpvid
pip install --no-deps -e .

# Create the machine-local configuration (repository-relative resource paths).
cp config/local.example.yaml config/local.yaml

# Download the five archives from the Zenodo resource record, verify and
# extract them under resources/. Archive names, SHA-256 values and steps are
# listed in resources/EXTERNAL_RESOURCES.md:
# sha256sum -c SHA256SUMS.txt
# for a in lpvid-0.2-*.tar.gz; do tar xzf "$a"; done

lpvid doctor        # missing large-resource entries are expected before placement
lpvid prepare       # run this after the adapter package is in place
```

Full requirements, resource deployment and troubleshooting:
[docs/INSTALL.md](docs/INSTALL.md).

## Quick start

```bash
# contigs input, example data
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo

# paired-end reads input, example data
lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo

# proteins input
lpvid run-proteins --input proteins.faa --out runs/protein_demo --sample demo

# inspect and resume
lpvid status runs/demo_contigs
lpvid resume runs/demo_contigs
```

`--gpu-id auto` is the default and selects the GPU with the most free memory.
Use `--gpu-id 0` to pin a device. Run `lpvid` from the repository root so that
relative resource paths resolve; absolute paths in `config/local.yaml` also
work.

A five-minute walkthrough with expected results is in
[docs/QUICKSTART.md](docs/QUICKSTART.md).

## Outputs

Each run directory keeps intermediate stage outputs and a final `results/`
delivery:

- `results/main_results_report_zh.html` and
  `results/main_results_report_en.html`: self-contained bilingual reports.
- `results/phylogenetic_explorer.html`: interactive placement tree browser.
- `results/tables/*.tsv`: per-protein detection, contig-level calls,
  representatives, motifs, supergroup assignments, placements, abundance
  (reads route), risk notices and neighbor evidence, plus `run_summary.json`
  and `delivery_manifest.json`.
- `results/sequences/`: RdRP-positive ORFs, source contigs and representative
  ORFs.
- `results/trees/<supergroup>.nwk`: backbone Newick used by the visualization.
- `results/README.md`: generated bilingual index of the delivered files.

Column definitions are in [docs/RESULTS.md](docs/RESULTS.md).

## Resources

Model weights and large databases are distributed separately from the source
repository. The required layout, inventory and verification rules are in
[resources/README.md](resources/README.md),
[resources/EXTERNAL_RESOURCES.md](resources/EXTERNAL_RESOURCES.md) and
[resources/RESOURCE_MANIFEST.tsv](resources/RESOURCE_MANIFEST.tsv).

## Documentation

| Document | Content |
|---|---|
| [docs/INSTALL.md](docs/INSTALL.md) / [zh-CN](docs/INSTALL.zh-CN.md) | Requirements, environment, resources, deployment validation |
| [docs/QUICKSTART.md](docs/QUICKSTART.md) / [zh-CN](docs/QUICKSTART.zh-CN.md) | Example-data walkthrough with expected results |
| [docs/RUNBOOK.md](docs/RUNBOOK.md) | Routes, defaults, outputs, resume, troubleshooting |
| [docs/RESULTS.md](docs/RESULTS.md) | Result table reference |
| [resources/README.md](resources/README.md) / [zh-CN](resources/README.zh-CN.md) | Resource inventory, layout and verification |
| [examples/README.md](examples/README.md) / [zh-CN](examples/README.zh-CN.md) | Bundled example data |
| [supergroups_readme.md](supergroups_readme.md) | Supergroup scheme and backbone coverage |

## Updating

Releases are tagged. This distribution corresponds to tag `0.2`. To update
a clone:

```bash
git fetch origin
git pull origin main          # or: git checkout <newer-tag>
conda env update -f environment.yml -n lpvid
pip install --no-deps -e .
lpvid doctor
```

Bundled resources are updated with the source repository. Downloadable model,
adapter and index resources are published as a versioned resource record;
update them according to the matching release notes.

## License

lpvid is released under the Apache License 2.0; see [LICENSE](LICENSE).
Bundled components and their license texts are listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Model weights, adapters,
reference data and databases are data assets with terms stated on the
accompanying resource record.

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff). Cite the
repository and the resource record when you use lpvid or its bundled resources.
