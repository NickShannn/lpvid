# Third-party notices

This file records the provenance and licensing of components bundled with
lpvid. The project source code is distributed under the Apache License 2.0
(see `LICENSE`).

## Project components

The following components are maintained as part of the lpvid project. No
separate upstream license file is included for them; they are distributed
under the project Apache License 2.0.

| Component | Path | Purpose |
|---|---|---|
| fastorf | `third_party/fastorf/` | Multithreaded ORF prediction implemented in Rust |
| placement_view | `third_party/placement_view/` | Offline phylogenetic tree explorer |
| lucaphylo | `third_party/lucaphylo/` | Stage-4 phylogenetic placement inference code |

## Open-source libraries

The following libraries are redistributed unmodified. Their license texts are
included in the repository at the paths shown.

| Library | Version | License | License text |
|---|---|---|---|
| Chart.js | 4.4.7 | MIT | `third_party/report_assets/licenses/CHARTJS-LICENSE` |
| Tabulator | bundled release | MIT | `third_party/report_assets/licenses/TABULATOR-LICENSE` |
| jsPDF | bundled release | MIT | `third_party/placement_view/vendor/LICENSE-jspdf` |
| svg2pdf.js | bundled release | MIT | `third_party/placement_view/vendor/LICENSE-svg2pdf` |

## Models, reference data and databases

Model weights, LoRA adapters, backbone trees, taxonomy tables, reference
sequences and databases are data assets, not source code. The Apache License
2.0 does not automatically apply to them. Their redistribution terms are stated
on the accompanying resource record; see
[`resources/EXTERNAL_RESOURCES.md`](resources/EXTERNAL_RESOURCES.md).

## Chinese summary

- 项目源码采用 Apache License 2.0。
- `fastorf`、`placement_view`、`lucaphylo` 为本项目维护的组件，随项目许可分发。
- Chart.js、Tabulator、jsPDF、svg2pdf.js 为 MIT 许可的第三方库，许可证全文已随仓库提供。
- 模型权重、adapter、骨架树、taxonomy 与参考数据库的再分发条款以资源记录为准，
  不自动继承源码许可。
