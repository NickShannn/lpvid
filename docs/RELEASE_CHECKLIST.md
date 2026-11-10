# lpvid release checklist

Version: 0.2

Status as of 2026-11-11. Checked items are complete for this release; unchecked
items remain open. This checklist covers release governance and external
distribution, not routine test execution.

## License and upstream notices

- [x] Project license selected: Apache License 2.0 (`LICENSE`).
- [x] Bundled components and license texts recorded in
  `THIRD_PARTY_NOTICES.md`.
- [x] Chart.js, Tabulator, jsPDF and svg2pdf.js license texts retained under
  `third_party/`.
- [ ] Confirm redistribution terms for model weights, LoRA adapters, backbone
  trees and reference databases, and record them on the resource record.

## Resource distribution

- [x] Small runtime resources bundled in git with file-level SHA-256 in
  `resources/RESOURCE_MANIFEST.tsv`.
- [x] Downloadable resource archives prepared with archive-level SHA-256 in
  `resources/EXTERNAL_RESOURCES.md`.
- [x] Resource layout matches `config/local.example.yaml` and
  `docs/INSTALL.md`.
- [ ] Publish the resource record and insert its URL/DOI in
  `resources/EXTERNAL_RESOURCES.md` and `resources/EXTERNAL_RESOURCES.zh-CN.md`.
- [ ] Record the resource record version in each run's `run_manifest.json`.
- [x] Large external items marked `verify_after_download` in the manifest.
- [ ] Provide a public download channel for the RefSeq protein DIAMOND database
  (about 214 GB, beyond the resource record limit).

## Repository and automation

- [x] Dedicated public repository created:
  `https://github.com/NickShannn/lpvid`.
- [x] `.gitignore` excludes run data, caches, local configuration, large models
  and databases while retaining source, templates and the fastorf binary.
- [x] No private configuration, run data, caches or model weights in the source
  tree.
- [ ] Configure continuous integration for Python tests, JavaScript core tests
  and documentation link checks.
- [ ] Separate resource-backed integration tests from resource-free unit tests
  in public CI.

## Independent-user acceptance

- [x] Source archive extracted into an empty directory and tested: 182 passed,
  2 explicit resource skips; 59 JavaScript tests passed.
- [x] Example contigs and reads runs completed with the final report template
  and passed the bilingual and tree verifiers.
- [x] Chinese and English reports opened in a browser; desktop and narrow-screen
  layouts checked.
- [ ] Repeat environment creation, `lpvid doctor`, `lpvid prepare` and one
  sample run on a clean machine or container using only public resources.
- [x] Supported platform, Python, CUDA and component versions documented in
  `docs/INSTALL.md` and `docs/INSTALL.zh-CN.md`.

## Optional follow-up

- [ ] Resource downloader or version-management command.
- [ ] Container image or lightweight CI runtime.
- [ ] Self-contained wheel (editable source installation is supported).
