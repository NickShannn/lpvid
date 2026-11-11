# lpvid：安装指南

Version: 0.2

本指南带领一台全新 Linux 机器从零完成 lpvid 的部署校验。源码与大型资源
（模型权重、数据库）分开安装。

> English version: [INSTALL.md](INSTALL.md)

## 1. 系统要求

| 项目 | 要求 |
|---|---|
| 操作系统 | Linux x86_64（glibc；已验证 Ubuntu 20.04+ 与 CentOS 8+） |
| CPU | 建议 8 核以上；组装与数据库工具均为多线程 |
| 内存 | 建议 64 GB；小样本最低 32 GB |
| GPU | NVIDIA GPU，显存至少 16 GB；建议 24 GB |
| NVIDIA 驱动 | 550.x 或更新；用 `nvidia-smi` 确认 |
| CUDA 运行时 | 12.8，随 PyTorch wheel 提供，无需本地 CUDA toolkit |
| 磁盘 | 仓库与环境约 15 GB，另加资源 |
| 软件 | conda 或 mamba（Miniconda 或 Miniforge）、git、pip |

资源需要可观的额外磁盘空间，其中蛋白参考数据库占用最大。发布时的文件级
清单与体量记录在
[resources/RESOURCE_MANIFEST.tsv](../resources/RESOURCE_MANIFEST.tsv)。

Node.js 仅在运行可选的树浏览器 JavaScript 测试时需要，不是运行工作流的
必要依赖。

## 2. 克隆仓库

```bash
git clone https://github.com/NickShannn/lpvid.git
cd lpvid
git status        # 分支或 tag：0.2
```

## 3. 创建 conda 环境

统一环境 `lpvid` 包含 Python 3.10、所需生信工具和 Python 依赖，其中
PyTorch 2.8.0（CUDA 12.8）。

```bash
# conda
conda env create -f environment.yml -n lpvid
conda activate lpvid

# 或 mamba（解算更快）
mamba env create -f environment.yml -n lpvid
conda activate lpvid
```

环境包含 `fastp`、`bowtie2`、`MEGAHIT`、`samtools`、`DIAMOND`、`CD-HIT` 与
`seqkit`。`fastorf` 二进制随仓库提供，位于 `third_party/fastorf/bin/`。

`environment.yml` 已定义 `name: lpvid`；`-n lpvid` 为显式写法。所有直接
依赖均锁定到精确版本；`config/requirements-lpvid.txt` 另行记录完全锁定的
传递依赖，需要逐位复现环境时执行：

```bash
pip install -r config/requirements-lpvid.txt
```

验证 GPU 栈与一个核心工具：

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# 期望：2.8.0 True
diamond --version
```

若 `torch.cuda.is_available()` 为 `False`，用 `nvidia-smi` 检查并升级 NVIDIA
驱动。若镜像源提供 CPU 版 wheel，请从 CUDA 12.8 索引安装 PyTorch：

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128
```

## 4. 安装工作流

```bash
pip install --no-deps -e .
lpvid --version      # 0.2
lpvid --help
```

`pip install -e .` 注册 `lpvid` 命令并保持源码可编辑；仓库根目录的 `./lpvid`
启动器同样可用。不要把模型或数据库安装进 site-packages，它们通过配置引用。

## 5. 放置资源

体积较小的运行资源已随 git 仓库内置，包括 RdRP 参考库、24 个骨架树与
taxonomy 表、分类启动脚本/标签表/运行时模型文件，以及两个 RefSeq taxid
白名单。需要单独获取的仅包括大型模型 checkpoint、24 组 LoRA 放置 adapter、
完整 RefSeq 蛋白 DIAMOND 数据库和 SILVA rRNA 索引。请按以下目录放入
`resources/`；权威清单见
[resources/RESOURCE_MANIFEST.tsv](../resources/RESOURCE_MANIFEST.tsv)，可下载
压缩包、内容与 SHA-256 见
[resources/EXTERNAL_RESOURCES.zh-CN.md](../resources/EXTERNAL_RESOURCES.zh-CN.md)。

```text
resources/
  phylo/
    backbone_dir/<超群>/{*.nwk,*.fas}   # 仓库内置
    adapter_dir/<超群>/...              # Zenodo 资源记录
    taxonomy/<超群>.csv                  # 仓库内置
  rdrp/
    rdrp_reference.dmnd                  # 仓库内置
    rdrp_reference.fas                   # 仓库内置
    rdrp_reference.csv                   # 仓库内置
  models/
    rdrp_detector/                       # Zenodo 资源记录
    lucaphylo/                           # Zenodo 资源记录
    classifier/
      pipeline/                          # 仓库内置
      encoder/                           # Zenodo 资源记录
      model/                             # 仓库内置
  blastdb/refseq_protein/
    refseq_protein.dmnd                  # 单独分发（体积大）
    TAXID_viruses.txt                    # 仓库内置
    TAXID_cellular.txt                   # 仓库内置
  rrna/
    silva_ssu.*                          # Zenodo 资源记录
    silva_lsu.*                          # Zenodo 资源记录
```

下载后先用 `resources/EXTERNAL_RESOURCES.zh-CN.md` 中的 SHA-256 校验压缩包，
再用清单核对文件。仓库内置文件带 SHA-256；压缩包内单个文件在清单中标记为
`verify_after_download`。

## 6. 配置本机路径

```bash
cp config/local.example.yaml config/local.yaml
```

`config/local.yaml` 是本机私有配置，已被 .gitignore 排除。模板默认使用第 3 步
创建的统一环境；若 conda 安装位置不同，请把 `/opt/conda/envs/lpvid` 前缀改为
实际路径（`conda info --base` 可查看根目录）。模板也默认使用第 5 步的
`resources/` 布局。需要调整科学参数或运行参数时，写在该文件中，不要修改
`config/default.yaml`。

请在仓库根目录运行 `lpvid`，使模板中的相对路径正确解析；或在
`config/local.yaml` 中使用绝对路径。

若预编译的 `third_party/fastorf/bin/fastorf` 无法在目标平台运行，用 Rust
重新构建：

```bash
cd third_party/fastorf
cargo build --release
mkdir -p bin && cp target/release/fastorf bin/fastorf
cd ../..
```

## 7. 校验部署

```bash
lpvid doctor      # 环境、工具、模型与数据库逐项检查，全部应为 [OK]
lpvid prepare     # 需骨架树、taxonomy 与 adapter 就位；构建本地参考资产
./scripts/validate.sh   # 单元测试与树浏览器核心测试
```

`prepare` 是确定性流程：输入不变则不会重建。它不会从互联网下载任何资源。

## 8. 首次运行

```bash
# contigs，内置示例数据
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo --gpu-id 0

# 双端 reads，内置示例数据
lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo --gpu-id 0
```

预期结果与五分钟演练见 [QUICKSTART.zh-CN.md](QUICKSTART.zh-CN.md)。

### 输出解读

每个运行目录包含：

```text
run_manifest.json              输入、合并配置与资源指纹
state.json / checkpoint.json    逐阶段状态与签名
step*/                         各阶段中间产物与日志
results/                       交付目录：HTML 报告、tables/、sequences/、trees/
```

`results/` 中的关键文件：

| 路径 | 内容 |
|---|---|
| `main_results_report_zh.html` / `main_results_report_en.html` | 自包含中英双语主报告 |
| `phylogenetic_explorer.html` | 交互式放置树，支持 SVG、PNG、PDF、Newick 导出 |
| `tables/orf_predictions.tsv` | 逐蛋白 RdRP 概率与检出判定 |
| `tables/contig_level_calls.tsv` | 逐 contig 判定；任一 ORF 阳性即 contig 阳性 |
| `tables/supergroup_assignments.tsv` | 180 超群分类及参考匹配或分类证据 |
| `tables/placements.tsv` | LucaPhylo top-5 骨架近邻与距离 |
| `tables/abundance.tsv` | 仅 reads 路线：mapped reads、RPM、深度与覆盖度 |
| `sequences/rdrp_positive_*.fna` / `.faa` | RdRP 阳性源 contig 与 ORF |
| `results/README.md` | 按实际交付文件生成的双语索引 |

逐列字段参考：[RESULTS.md](RESULTS.md)。

## 9. 更新

版本以 tag 发布，当前为 `0.2`。

```bash
git fetch origin
git pull origin main        # 跟随主干更新
# 或切换到更新的 tag：git checkout <tag>
conda env update -f environment.yml -n lpvid
pip install --no-deps -e .
lpvid doctor
```

仓库内置资源随源码一起更新；可下载的模型、adapter 与索引资源以版本化资源记录
发布，请按对应版本说明更新。

## 10. 故障排查

| 现象 | 处理 |
|---|---|
| `lpvid: command not found` | 激活 `lpvid` 环境；重跑 `pip install --no-deps -e .`；或在仓库根目录用 `./lpvid` |
| `doctor` 报 `missing` | 按输出修正 `config/local.yaml` 对应项；不要关闭检查 |
| `torch.cuda.is_available()` 为 `False` | 升级 NVIDIA 驱动到 550 或更新；从 CUDA 12.8 索引重装 PyTorch |
| 找不到 `fastp`、`bowtie2`、`MEGAHIT` 或 `samtools` | 它们位于 `lpvid` 环境内；确认环境已激活且 `config/local.yaml` 的环境路径指向该环境 |
| 参考数据库任务慢 | 大型数据库扫描受 IO 限制；小查询在完整库上也可能需要约 10 分钟。不要并行启动多个全库任务 |
| `fastorf` 二进制无法执行 | 按第 6 步从源码重建，并检查 glibc 版本 |
| GPU 显存不足 | 在 `config/local.yaml` 调低模型批大小，或用 `--gpu-id` 换更大显存的卡 |
| 运行中断 | `lpvid status <run_dir> --json` 定位失败阶段；修复后 `lpvid resume <run_dir>` 只重跑受影响阶段 |

## 11. 卸载

```bash
conda deactivate
conda env remove -n lpvid
pip uninstall lpvid
```

运行结果、资源与缓存独立于安装，保留在各自位置。

## 12. 许可证与引用

lpvid 采用 Apache License 2.0 许可证，全文见仓库根目录 `LICENSE`；内置组件及
其许可证全文见 `THIRD_PARTY_NOTICES.md`。模型权重、adapter、参考数据与数据库
属于数据资产，其条款以随附资源记录为准。

引用元数据见 `CITATION.cff`。
