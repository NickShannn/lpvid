# fastorf 使用说明

**fastorf** 是 NCBI ORFfinder 的多线程重写（Rust），用于在核酸序列中查找开放阅读框（ORF）。
与官方 ORFfinder 0.4.3 的输出坐标与 partial 标注已在 16 组对拍用例上验证完全一致
（见 `xval/xval_all.sh`），在此之上新增：并行加速、ID 规范化、TSV 映射表、GFF3 输出。

- 二进制：`fastorf/target/release/fastorf`（单文件，无运行时依赖，可直接拷贝到其他 x86_64 Linux）
- 项目文档：`PLAN.md`（里程碑）、`DECISIONS.md`（全部决策与对拍结论）、`CODING.md`（编码准则）

---

## 1. 快速上手

```bash
B=~/orf_rework/fastorf/target/release/fastorf

# (1) 最简: 官方默认等价 (gc=1, ATG+替代起始, ml=75, 双链), 蛋白 FASTA 输出到 stdout
$B --in transcripts.fa

# (2) 典型工作流: 蛋白 + CDS + ID映射表 + GFF3, 输入支持 .gz
$B --in transcripts.fa.gz \
    --out orfs.faa \
    --outfmt 1 --out orfs.fnt \
    --tsv orf_map.tsv \
    --gff3 orfs.gff3
# 注: --outfmt 与 --out 是一对; 需要两种 FASTA 时跑两次, 或见下方多格式说明

# (3) 模拟官方某个具体调用
#   官方: ORFfinder -in x.fa -g 11 -s 0 -ml 300 -n true -strand plus -outfmt 0 -out y.faa
$B --in x.fa --gc 11 --start 0 --ml 300 --nested drop --strand plus --out y.faa
```

> **多格式输出说明**：`--outfmt` 只控制 `--out` 主输出的格式（一次一种，与官方一致）；
> 蛋白/CDS 的推荐输出方式是 `--faa`/`--fnt` 独立开关，一轮全出（见上方示例）。

## 2. 参数完整表

### 输入

| 参数 | 默认 | 说明 |
|------|------|------|
| `--in FILE` | 必填 | FASTA；支持 `.gz` 自动解压；`-` 表示 stdin |
| `--begin N` | 1 | 只处理片段起点（1-based，含） |
| `--end N` | 0 | 片段终点（0 = 到序列末尾）；报告坐标仍为全长序列坐标（与官方一致） |
| `--circular` | 关 | 环状序列，允许跨原点 ORF（官方此功能 under development，fastorf 为完整实现） |

### 搜索

| 参数 | 默认 | 说明 |
|------|------|------|
| `--gc N` | 1 | NCBI 遗传密码表：1-6, 9-16, 21-33（数据直接来自 NCBI gc.prt v4.6） |
| `--start 0\|1\|2` | 1 | 0 = 仅 ATG；1 = ATG + 该密码表的官方替代起始密码子；2 = 任意 sense 密码子 |
| `--start-codons A,B` | — | 显式指定起始密码子集合（覆盖 `--start`，超集功能） |
| `--stop-codons A,B` | — | 显式指定终止密码子集合（覆盖表内终止，超集功能） |
| `--ml N` | 75 | 最小 ORF 长度（**不含**终止密码子的 nt 数，与官方口径一致）；< 30 时打印警告但按给定值执行（官方是静默改为 30）。序列（片段）长度小于 `--ml` 时整条跳过（与官方一致） |
| `--nested all\|drop` | all | `all` = 官方默认行为（见下方"语义"）；`drop` = 官方 `-n true`，额外删除完全嵌套于其他 ORF 内的（跨框跨链判定，与官方实测一致） |
| `--strand both\|plus\|minus` | both | 链选择 |
| `--frames 1,3,4-6` | 全部 | 精确到框：1-3 正向、4-6 反向（官方此版无此功能，超集） |

### 输出

| 参数 | 默认 | 说明 |
|------|------|------|
| `--out FILE` | stdout | 主输出文件 |
| `--outfmt 0\|1\|2\|3` | 0 | 0 = ORF 蛋白 FASTA；1 = CDS 核酸 FASTA；2 = Text ASN.1；3 = feature table |
| `--tsv FILE` | — | **ID/坐标映射表**（见第 4 节，本项目核心新增） |
| `--gff3 FILE` | — | GFF3 输出 |
| `--id-mode numeric\|keep` | numeric | numeric = 重编号（`SEQ000001_ORF2`）；keep = 清洗保留原 ID |
| `--prefix STR` | SEQ | numeric 模式的 ID 前缀 |
| `--threads N` | 全部核 | 并行线程数 |

### 其他

`-V/--version` 版本；`-h/--help` 用法。

## 3. ORF 查找语义（全部经官方 0.4.3 对拍验证）

| 语义 | 规则 |
|------|------|
| ORF 定义 | 起始密码子 → 下游第一个终止密码子 |
| 同框嵌套 | **不输出**：同一框同一 stop 区间内只报第一个起始密码子的 ORF（后续 start 嵌套跳过）。这是官方默认行为 |
| `--nested drop` | 在上表基础上，再删除被其他 ORF（任何框、任何链）坐标完全包含的 ORF |
| `-s 1` 替代起始集合 | 按所选 `--gc` 表的 NCBI 官方数据（gc.prt）。例：表 1 = TTG/CTG/ATG（**不含** GTG，与网页版 wprintgc 的展示不同，以工具行为为准） |
| ml 长度口径 | 不含终止密码子的 nt 数 |
| 3' partial | 有起始无终止的 ORF 照常输出（延伸到窗口尾，截到最后一个完整密码子） |
| partial 标注 | 5' 侧 = 起点紧贴搜索窗口头（序列头或 `-b` 片段头）；3' 侧 = 无终止密码子。有 stop 的 ORF 即便终点贴序列尾也标 complete |
| 起点翻译 | 起始密码子译为 M；**但起点贴窗口头（5' partial）时保留原义**（TTG/CTG→L），219 万条蛋白全量对拍确认 |
| 输出顺序 | 序列按输入顺序；序列内按框序 +1,+2,+3,−1,−2,−3，框内按链上位置升序；与线程数无关 |
| 不输出 | 无起始密码子的窗口头读码（即便 `-s 2`） |

## 4. 输出文件格式

### 4.1 TSV 映射表（`--tsv`）——ID 与坐标的权威对照

首行为 `#` 注释（版本与参数），随后表头 + 每 ORF 一行：

```
orf_id  orig_id  desc  strand  frame  start  end  nt_len  aa_len  start_codon  stop_codon  status
SEQ000001_ORF1  transcript_001  hypothetical protein  -  2  423  501  76  21  GTG  TAA  complete
SEQ000001_ORF2  transcript_001  hypothetical protein  -  76  423  348  115  GTG      3partial
```

| 列 | 含义 |
|----|------|
| `orf_id` | 规范化输出 ID（即所有输出文件中 FASTA/GFF3 的 ID） |
| `orig_id` / `desc` | 原始 FASTA ID 与描述；其中 tab/换行替换为空格，**其余字符（含非 ASCII）原样保留** |
| `strand` / `frame` | 链与该链上的框（1-3） |
| `start` / `end` | 1-based 闭区间，**相对原始全长序列**；minus 链同样 start < end，方向看 strand |
| `nt_len` | 链方向 5'→3' 长度，含终止密码子 |
| `aa_len` | 蛋白残基数（起始译 M，不含终止） |
| `start_codon` / `stop_codon` | 链方向的实际密码子；3' partial 的 stop 为空 |
| `status` | `complete` / `5partial` / `3partial` / `partial`（语义见上表） |

### 4.2 ID 规范化

- `numeric`（默认）：`{prefix}{序列号:6位}_ORF{orf序号}`，如 `SEQ000001_ORF2`。任何原始 ID 中的特殊字符只存在于 TSV 的 `orig_id` 列，输出序列文件中绝对干净。
- `keep`：白名单 `[A-Za-z0-9._-]`，其余字符折叠为单个 `_`；首尾 `_` 去除；清洗后为空回退编号；清洗后撞名自动加 `_2/_3` 后缀。
  例：`sp|P12345|XYZ_HUMAN` → `sp_P12345_XYZ_HUMAN`。

### 4.3 GFF3（`--gff3`）

标准 9 列；feature type 为 `ORF`，phase 0；`ID=` 规范化 orf_id，`orig_id=` 为 percent-encoding 后的原始 ID。

### 4.4 outfmt 0/1（蛋白/CDS FASTA）

header 格式：`>{orf_id} [{orig_id}] Nucleotide {start}..{end} strand {s} frame {f} ({status}, {n} aa)`。
CDS 为链方向 5'→3'（minus 链已反向互补），含终止密码子。

### 4.5 outfmt 2/3（ASN.1 / feature table）

近似官方格式（坐标、链、密码表编号一致）；如需与官方逐字节一致请告知，可按样本校准。

## 5. 性能

128 核服务器、5 万条序列 / 150 Mbp / 455 万 ORF / 1 GB 输出：

同一数据（5 万条 / 150 Mbp / 219 万 ORF / ~300 MB 输出）实测：

| 配置 | 耗时 |
|------|------|
| 官方 ORFfinder 0.4.3（单线程） | 339.67 s |
| fastorf 128 线程，**全流程含写盘** | **1.09 s** |
| 加速比 | **≈ 312×** |

（早期 455 万 ORF 版本数据见 DECISIONS.md T13：计算阶段 9.1x 并行扩展、云盘带宽 ~200 MB/s 为输出下限。）

结果跨线程数完全一致（md5 相同）。

## 6. 与官方 ORFfinder 的差异清单

1. **多线程**：官方单线程；fastorf 默认用满所有核，`--threads` 可控
2. **`--ml < 30`**：fastorf 打印警告并按给定值执行；官方静默改为 30
3. **ID 规范化 + TSV 映射表**：fastorf 新增；官方输出的 ID 直接取自输入
4. **`-id`（按 accession 联网取序列）**：fastorf 不实现，请上游 `efetch` 落盘
5. **`-c` 环状**：fastorf 完整实现（跨原点 ORF 用线性化坐标，`end` 可能大于序列长度，见 TSV）；官方标注 under development
6. **`--frames` 单框选择**、**`--start-codons`/`--stop-codons` 显式密码子**：fastorf 超集
7. 输出 header 文案、TSV/GFF3 为 fastorf 自有格式（坐标体系见第 4 节，与官方 0-based 表示法的换算：fastorf 1-based start = 官方 from + 1）

## 7. 常见问题

- **内存**：序列、坐标与格式化输出常驻内存，量级约为 `输入 + 2×输出`；超大输入请分批
- **输出量本身巨大时**（如随机序列测试），耗时下限是磁盘写带宽，与程序无关
- **管道使用**：`--in -` 读 stdin；主输出默认 stdout，stderr 输出统计行（方便 `2>/dev/null` 静默）
- **重跑对拍**：`bash xval/xval_all.sh`（需要官方 ORFfinder 参照二进制与 `official_out/` 内的官方参照输出）
