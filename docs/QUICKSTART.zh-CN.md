# lpvid：快速上手

Version: 0.2

基于内置示例数据的五分钟演练。默认安装已完成，且 `lpvid doctor` 全部为
`[OK]`；否则请先阅读 [INSTALL.zh-CN.md](INSTALL.zh-CN.md)。

> English version: [QUICKSTART.md](QUICKSTART.md)

示例数据位于 [`examples/`](../examples/README.zh-CN.md)：9 条 contig
（1 条 RdRP 阳性，8 条阴性背景）和模拟的双端 2x150 reads。它用于验证工作流
连通性，不能用于性能评估。

## 1. contigs 路线

```bash
cd lpvid          # 仓库根目录
conda activate lpvid

lpvid run-contigs \
  --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs \
  --sample demo \
  --gpu-id 0
```

结束时输出 `<状态>\t<运行目录>`。可在另一个终端用
`lpvid status runs/demo_contigs` 查看进度。

### 预期结果

| 检查项 | 位置 | 预期 |
|---|---|---|
| RdRP 阳性 ORF | `results/sequences/rdrp_positive_orfs.faa` | 非空；全部来自 `demo_rdrp_positive_1` |
| 阳性源 contig | `results/sequences/rdrp_positive_contigs.fna` | 仅 `demo_rdrp_positive_1` |
| contig 判定 | `results/tables/contig_level_calls.tsv` | `demo_rdrp_positive_1` 阳性；`demo_background_*` 阴性 |
| 超群分类 | `results/tables/supergroup_assignments.tsv` | 阳性 ORF 的代表序列经参考匹配划入 Nido 超群 |
| 系统发育放置 | `results/tables/placements.tsv` | 近邻位于 Nido 骨架树内；最近邻为 alphacoronavirus RdRP 参考序列 |
| 报告 | `results/main_results_report_zh.html` | 阳性 contig 带超群、放置与名录提示出现 |

用浏览器打开 `results/main_results_report_zh.html` 与
`results/phylogenetic_explorer.html`。两者均为单文件离线页面。

## 2. reads 路线

```bash
lpvid run-reads \
  --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads \
  --sample demo \
  --gpu-id 0
```

该路线增加质控、去 rRNA、组装与丰度估算，随后经过相同的识别、分类与放置
步骤。

预期结果：

- `results/tables/abundance.tsv`：与 `demo_rdrp_positive_1` 对应的重新组装
  contig，包含 mapped reads、RPM、平均深度与覆盖广度。示例 reads 在阳性
  基因组上约 25x 覆盖；ID 来自组装结果，与输入 contig 名不同。
- `results/tables/supergroup_abundance.tsv`：检测到的超群各一行。
- `demo_background_*` 可能被部分组装；预期保持 RdRP 阴性。

## 3. 蛋白路线（可选）

蛋白路线期望氨基酸 FASTA，主分析集要求每条序列至少 150 aa。请使用自己的
蛋白输入运行：

```bash
lpvid run-proteins \
  --input proteins.faa \
  --out runs/protein_demo --sample demo --gpu-id 0
```

蛋白输入没有源核酸，contig 级表格为空属于预期行为。

## 4. 恢复与查看

```bash
lpvid status runs/demo_contigs --json   # 逐阶段状态
lpvid resume runs/demo_contigs          # 只重跑受影响阶段
```

## 5. 示例数据的边界

- 示例输入用于验证端到端连通性、配置与报告，不构成灵敏度、特异性或运行
  时间基准；性能评估请使用真实评估队列。
- 模拟 reads 含约 0.8% 替换错误、平均 300 bp 片段；组装质量刻意保持一般。
- 背景 contig 为非病毒组装产物；其筛选偶尔给出 `KEEP_unknown` 属预期。

后续阅读：[RUNBOOK.md](RUNBOOK.md)（生产路线与默认参数）、
[RESULTS.md](RESULTS.md)（结果表定义）。
