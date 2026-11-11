# lpvid

Version: 0.2

lpvid 是一个命令行 RNA 病毒发现与表征工作流，面向宏转录组数据。它支持测序
reads、组装 contigs 或蛋白序列输入，识别含 RdRP 的序列，按 180 超群体系
分类，将代表序列放置到精选系统发育骨架树，对 reads 路线估算丰度，并输出
自包含的中英双语 HTML 报告。

## 主要能力

- 三类输入路线：双端或单端 reads、组装核酸 contigs、氨基酸 FASTA。
- RdRP 识别，输出逐蛋白与逐 contig 判定。
- 180 超群分类。
- 使用 LucaPhylo 在精选骨架树上进行系统发育放置，并提供独立交互式树浏览器。
- reads 路线丰度分析：mapped reads、RPM、平均深度与覆盖广度。
- 中英双语、完全离线的 HTML 报告，支持可排序表格、名录提示，以及
  SVG、PNG、PDF、Newick 导出。
- 可恢复运行：各阶段带指纹签名，`lpvid resume` 只重跑发生变化的阶段。

## 输入

| 路线 | 输入 | 命令 |
|---|---|---|
| reads | 双端或单端 FASTQ | `lpvid run-reads` |
| contigs | 组装核酸 FASTA | `lpvid run-contigs` |
| proteins | 氨基酸 FASTA | `lpvid run-proteins` |

## 安装摘要

```bash
git clone https://github.com/<org>/lpvid.git
cd lpvid

conda env create -f environment.yml -n lpvid
conda activate lpvid
pip install --no-deps -e .

# 生成本机私有配置（默认使用仓库相对资源路径）。
cp config/local.example.yaml config/local.yaml

# 从 Zenodo 资源记录下载 5 个压缩包，校验后解压到 resources/。
# 压缩包名称、SHA-256 与步骤见 resources/EXTERNAL_RESOURCES.zh-CN.md：
# sha256sum -c SHA256SUMS.txt
# for a in lpvid-0.2-*.tar.gz; do tar xzf "$a"; done

lpvid doctor        # 资源放置前，大型资源项显示 missing 属于预期
lpvid prepare       # adapter 资源包就位后执行
```

完整系统要求、资源部署与故障排查见
[docs/INSTALL.zh-CN.md](docs/INSTALL.zh-CN.md)。

## 快速运行

```bash
# contigs 输入，示例数据
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo

# 双端 reads 输入，示例数据
lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo

# 蛋白输入
lpvid run-proteins --input proteins.faa --out runs/protein_demo --sample demo

# 查看与恢复
lpvid status runs/demo_contigs
lpvid resume runs/demo_contigs
```

`--gpu-id auto` 为默认值，会选择空闲显存最多的 GPU；用 `--gpu-id 0` 固定设备。
请在仓库根目录运行 `lpvid`，使相对资源路径正确解析；在 `config/local.yaml`
中使用绝对路径同样可行。

包含预期结果的五分钟演练见 [docs/QUICKSTART.zh-CN.md](docs/QUICKSTART.zh-CN.md)。

## 输出

每个运行目录保留中间阶段产物与最终 `results/` 交付目录：

- `results/main_results_report_zh.html` 与
  `results/main_results_report_en.html`：自包含中英双语报告。
- `results/phylogenetic_explorer.html`：交互式放置树浏览器。
- `results/tables/*.tsv`：逐蛋白识别、contig 级判定、代表序列、motif、
  超群分类、系统发育放置、丰度（reads 路线）、名录提示与近邻证据，以及
  `run_summary.json` 与 `delivery_manifest.json`。
- `results/sequences/`：RdRP 阳性 ORF、源 contig 与代表 ORF。
- `results/trees/<超群>.nwk`：可视化使用的骨架树 Newick。
- `results/README.md`：按本次实际交付文件生成的双语索引。

逐列定义见 [docs/RESULTS.md](docs/RESULTS.md)。

## 资源

模型权重与大型数据库不随源码仓库分发。目录约定、清单与校验规则见
[resources/README.zh-CN.md](resources/README.zh-CN.md)、
[resources/EXTERNAL_RESOURCES.zh-CN.md](resources/EXTERNAL_RESOURCES.zh-CN.md) 与
[resources/RESOURCE_MANIFEST.tsv](resources/RESOURCE_MANIFEST.tsv)。

## 文档导航

| 文档 | 内容 |
|---|---|
| [docs/INSTALL.zh-CN.md](docs/INSTALL.zh-CN.md) / [EN](docs/INSTALL.md) | 系统要求、环境、资源、部署校验 |
| [docs/QUICKSTART.zh-CN.md](docs/QUICKSTART.zh-CN.md) / [EN](docs/QUICKSTART.md) | 示例数据演练与预期结果 |
| [docs/RUNBOOK.md](docs/RUNBOOK.md) | 路线、默认参数、输出、恢复、故障排查 |
| [docs/RESULTS.md](docs/RESULTS.md) | 结果表字段参考 |
| [resources/README.zh-CN.md](resources/README.zh-CN.md) / [EN](resources/README.md) | 资源清单、目录约定与校验 |
| [examples/README.zh-CN.md](examples/README.zh-CN.md) / [EN](examples/README.md) | 内置示例数据 |
| [supergroups_readme.md](supergroups_readme.md) | 超群体系与骨架树覆盖说明 |

## 更新

版本以 tag 发布，当前分发对应 tag `0.2`。更新方法：

```bash
git fetch origin
git pull origin main          # 或 git checkout <更新的 tag>
conda env update -f environment.yml -n lpvid
pip install --no-deps -e .
lpvid doctor
```

仓库内置资源随源码一起更新；可下载的模型、adapter 与索引资源以版本化资源记录
发布，请按对应版本说明更新。

## 许可证

lpvid 采用 Apache License 2.0 许可证，全文见 [LICENSE](LICENSE)。内置组件及
其许可证全文见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。模型权重、
adapter、参考数据与数据库属于数据资产，其条款以随附资源记录为准。

## 引用

引用元数据见 [CITATION.cff](CITATION.cff)。使用 lpvid 或其资源时，请同时引用
代码仓库与资源记录。
