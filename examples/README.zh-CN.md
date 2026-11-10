# 内置示例数据

Version: 0.2

用于验证 lpvid 安装的小型测试输入。它以一条已知阳性序列跑通完整工作流，
仅用于连通性检查，不能用于灵敏度、特异性或运行时间评估。

> English version: [README.md](README.md)

## 文件

| 文件 | 内容 |
|---|---|
| `contigs/demo_contigs.fa` | 9 条 contig，合计 39,577 bp：`demo_rdrp_positive_1`（27,452 bp，RdRP 阳性）与 `demo_background_1` 到 `demo_background_8`（703 到 2,238 bp，非病毒） |
| `reads/demo_R1.fastq.gz` | 2,612 对双端 reads，2x150，gzip 后约 114 KB |
| `reads/demo_R2.fastq.gz` | R1 的配对文件 |
| `src/extract_demo_contigs.py` | 生成 contigs 文件的脚本 |
| `src/simulate_demo_reads.py` | reads 模拟器（片段长度 N(300, 30)、0.8% 替换错误、Illumina 风格质量值、种子 20261110） |

阳性 contig 取自经 lpvid 运行验证的真实宏转录组组装：携带 RdRP 阳性 ORF，
被划入 Nido 超群，近邻为 alphacoronavirus 参考序列。背景 contig 为同一组装
中的阴性、非病毒序列。所有序列 ID 已重命名，不含任何样本标识。

## 使用

```bash
lpvid run-contigs --input examples/contigs/demo_contigs.fa \
  --out runs/demo_contigs --sample demo --gpu-id 0

lpvid run-reads --r1 examples/reads/demo_R1.fastq.gz \
  --r2 examples/reads/demo_R2.fastq.gz \
  --out runs/demo_reads --sample demo --gpu-id 0
```

## 预期结果

- `demo_rdrp_positive_1`：至少一个 ORF 判为 RdRP 阳性；contig 出现在
  `results/sequences/rdrp_positive_contigs.fna`；其代表序列经参考匹配划入
  Nido 超群，放置于 alphacoronavirus 骨架近邻
  （`results/tables/placements.tsv`）。
- `demo_background_*`：预期无阳性判定。
- reads 路线：阳性基因组重新组装后约 25x 覆盖，并在
  `results/tables/abundance.tsv` 中得到 RPM、深度与覆盖度行。

完整演练见 [docs/QUICKSTART.zh-CN.md](../docs/QUICKSTART.zh-CN.md)。

## 重新生成

两个脚本仅依赖 Python 标准库：

```bash
python3 src/extract_demo_contigs.py \
    --positive-fna <validation_run>/sequences/rdrp_positive_contigs.fna \
    --calls-tsv <validation_run>/tables/contig_level_calls.tsv \
    --source-contigs <assembly>.contigs.fa --out contigs/demo_contigs.fa

python3 src/simulate_demo_reads.py --contigs contigs/demo_contigs.fa \
    --r1 reads/demo_R1.fastq.gz --r2 reads/demo_R2.fastq.gz
```

`extract_demo_contigs.py` 需要访问原始验证运行与组装文件；随包提供该脚本是
为了来源可溯，不是用于再次分发源数据。
