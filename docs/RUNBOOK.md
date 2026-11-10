# lpvid 运行手册

Version: 0.2

## 1. 运行方式

lpvid 是本地命令行工作流，不启动 Web 服务。每次运行使用独立输出目录，输入、
配置、资源指纹、阶段状态和结果都记录在该目录中。

能力概览：

```text
reads / contigs / proteins
        │
        ▼
   RdRP 识别
        │
        ▼
   超群分类
        │
        ▼
 系统发育放置
        │
        ▼
 丰度分析（仅 reads）
        │
        ▼
 双语报告与结果表
```

请在仓库根目录运行 `lpvid`，或在 `config/local.yaml` 中使用绝对路径。新输入
必须使用新的输出目录；非空目录不会被静默覆盖。

## 2. 预检与资产准备

```bash
cd lpvid   # 仓库根目录
./lpvid doctor
./lpvid prepare
```

`doctor` 检查环境、工具、模型与数据库，所有行应为 `[OK]`。`prepare` 校验
精选骨架树资源并构建本地参考资产；相同输入重复运行不会重建，且不会从互联网
下载资源。

可选地预热系统发育放置缓存：

```bash
./lpvid prepare --prewarm
```

预热会执行真实的放置计算，耗时较长，建议提前完成，不要在业务运行前临时执行。

可用一个命令完成静态预检、单元测试和树浏览器核心测试：

```bash
./scripts/validate.sh
```

大模型、蛋白参考数据库与 rRNA 索引通过 `config/default.yaml` 或
`config/local.yaml` 的路径引用，目录约定见
[资源目录约定](../assets/RESOURCE_LAYOUT.md)。

参考数据库相关步骤受 IO 限制：即使查询集很小，完整库扫描也可能需要约
10 分钟。请为该步骤预留时间，不要并行启动多个全库任务。

## 3. 三种输入

### 蛋白

```bash
./lpvid run-proteins \
  --input /absolute/path/proteins.faa \
  --out /absolute/path/runs/protein_demo \
  --sample demo --gpu-id 1
```

蛋白输入直接进入 RdRP 识别、超群分类与系统发育放置，不做核酸层面的统计；
contig 级表格在该路线为空是预期行为。主分析集要求每条序列至少 150 aa。

### 核酸 contigs

```bash
./lpvid run-contigs \
  --input /absolute/path/contigs.fna \
  --out /absolute/path/runs/contig_demo \
  --sample demo --gpu-id 1
```

### 双端 reads

```bash
./lpvid run-reads \
  --r1 /absolute/path/sample_R1.fastq.gz \
  --r2 /absolute/path/sample_R2.fastq.gz \
  --out /absolute/path/runs/reads_demo \
  --sample demo --gpu-id 1
```

### 单端 reads

```bash
./lpvid run-reads \
  --r1 /absolute/path/sample.fastq.gz \
  --out /absolute/path/runs/reads_se_demo \
  --sample demo --gpu-id 1
```

不要为单端数据伪造 R2。`--gpu-id auto` 是默认值，会选择当前空闲显存最多的
GPU；共享服务器上运行时建议显式指定。

## 4. 默认参数

| 环节 | 默认行为 |
|---|---|
| 输入与主分析长度 | 主分析集要求 ≥150 aa；更短序列不进入主分析 |
| RdRP 识别阈值 | 0.5 |
| 超群分类 | 180 超群体系；结果区分为参考匹配与分类 |
| 代表序列 | 每个非冗余分组保留一条代表序列，选择规则确定 |
| 系统发育放置 | LucaPhylo；最大序列长度 2048；输出 top-5 近邻与距离 |
| 丰度 | `bowtie2 --very-sensitive-local`；RPM 为 primary mapped reads / clean reads × 1e6；双端分母为 R1+R2 reads |
| 线程 | 由 `config/default.yaml` 定义，可在本机 YAML 覆盖 |
| GPU | `auto` 选择空闲显存最多的卡；共享服务器建议显式指定 |

进阶阈值、筛选规则与工具参数在 `config/default.yaml` 中；需要调整时写入
`config/local.yaml`，不要直接修改默认文件。请保留对应运行目录的
`run_manifest.json`。

## 5. 输出与恢复

每次运行有独立目录：

```text
run_manifest.json       输入、配置与资源记录
state.json              各阶段原子更新的状态
checkpoint.json         阶段签名与产物记录
step*/                  各阶段中间输出与日志
results/                最终 FASTA、TSV 与 HTML 交付结果
step*_report/           summary.md、summary.json、图和中间摘要
```

命令示例：

```bash
./lpvid status runs/contig_demo
./lpvid status runs/contig_demo --json
./lpvid resume runs/contig_demo
```

恢复时只有阶段代码、输入、参数、资源签名和产物均一致的阶段才会跳过。上游
变化会将下游标为 stale 并重跑，不会仅因同名文件存在就复用结果。

主要结果：

- `results/tables/orf_predictions.tsv`：逐蛋白 RdRP 识别概率与判定。
- `results/tables/contig_level_calls.tsv`：contig 级检出汇总。
- `results/tables/prefilter_calls.tsv`：蛋白筛选判定与分析集标记。
- `results/sequences/rdrp_positive_orfs.faa`：主分析集阳性 ORF；QA 抽样不混入。
- `results/sequences/representative_orfs.faa` 与
  `results/tables/representative_map.tsv`：代表序列及全部成员映射。
- `results/tables/motifs.tsv`：RdRP motif 位置。
- `results/tables/supergroup_assignments.tsv`：180 超群分类及参考匹配或
  分类证据。
- `step*_classifier/classifier_assignments.tsv`：分类
  的置信度与候选超群（阶段目录，主报告读取后展示）。
- `results/tables/placements.tsv`：LucaPhylo 实际输出的 top-5 近邻与距离。
- `results/tables/abundance.tsv`：仅 reads 路线，contig 级 RPM、深度与覆盖度。
- `results/tables/supergroup_abundance.tsv`：仅 reads 路线，超群级丰度汇总。
- `results/tables/risk_warning.tsv` 与 `results/tables/neighbor_evidence.tsv`：
  类群归属、名录提示与结构化近邻证据。
- `results/main_results_report_zh.html` / `main_results_report_en.html`：
  独立自包含双语主报告。
- `results/phylogenetic_explorer.html`：当前运行专属离线系统发育树浏览器。
- `results/README.md`：按本次实际交付文件生成的双语索引。

所有表格的逐列解释见 [RESULTS.md](RESULTS.md)。

## 6. 科学解释边界

- RdRP 识别概率是模型输出，不等同于实验确认或完整病毒基因组证明。
- 超群分类中的参考序列一致性、比对 E 值和比对分数是证据指标，不是分类概率。
- 系统发育浏览器中的挂枝是最近邻近似，不是枝级似然放置。
- 名录提示只锚定最近邻所属科的类群归属与名录先验，不是风险分级或病原确证。
- RPM 是文库内相对丰度，不是 copies/mL 或绝对载量。
- 内置示例数据与骨架树示例用于连通性检查，不作为独立科学准确率证据。

## 7. 故障排查

1. 先运行 `./lpvid status RUN_DIR --json`，定位 failed 阶段。
2. 查看该阶段目录中的 `run.log`、工具专用日志和 `commands.jsonl`。
3. 修复环境或配置后运行 `resume`；不要删除失败证据。
4. 新输入必须换新输出目录；非空目录不会被静默覆盖。
5. 大库 IO 慢时不要同时启动多个参考数据库任务；GPU 推理与重 IO 步骤尽量错峰。
