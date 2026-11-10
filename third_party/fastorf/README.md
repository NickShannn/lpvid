# fastorf

NCBI ORFfinder 的多线程 Rust 重写: 功能对齐 + 并行加速 + ID 规范化 + 映射表输出。

## 构建

```bash
cargo build --release   # 产物: target/release/fastorf, 无运行时依赖
```

## 用法

```bash
# 基本等价 ORFfinder 默认行为 (蛋白 FASTA 到 stdout)
fastorf --in genes.fa

# 常用: 蛋白 + CDS + 映射表 + GFF3, 重编号 ID
fastorf --in genes.fa.gz --out orfs.faa --outfmt 0 \
        --tsv map.tsv --gff3 orfs.gff3

# 模拟官方: -s 0 (仅 ATG), -n true (去嵌套), 密码表 11, 最小 300 nt
fastorf --in genes.fa --start 0 --nested drop --gc 11 --ml 300 --out orfs.faa
```

## 参数 (与官方 ORFfinder 对照)

| 官方 | fastorf | 说明 |
|------|---------|------|
| `-in FILE` | `--in FILE` | FASTA, 支持 .gz 与 `-` (stdin) |
| `-b N` / `-e N` | `--begin N` / `--end N` | 片段 (报告坐标仍是全长序列坐标, 待 M4 校准) |
| `-c` | `--circular` | 真正实现跨原点 ORF (官方 under development) |
| `-g N` | `--gc N` | 遗传密码表 1-6,9-16,21-33 |
| `-s 0\|1\|2` | `--start 0\|1\|2` | 1 = ATG + 该表官方 initiator (gc.prt sncbieaa 'M' 位) |
| (无) | `--start-codons ATG,GTG` | 显式起始密码子 (超集) |
| (无) | `--stop-codons TAA,TAG` | 显式终止密码子 (超集) |
| `-ml N` | `--ml N` | 最小 ORF 长度 nt (含终止密码子); <30 显式警告但照用 (官方静默改 30) |
| `-n t\|f` | `--nested all\|drop` | drop = 删除完全嵌套 ORF (坐标包含判定, 不分框/链, 待 M4 校准) |
| `-strand` | `--strand both\|plus\|minus` | 另有 `--frames 1,3,4-6` 超集 |
| `-outfmt 0-3` | `--outfmt 0-3` | 0=ORF 蛋白 FASTA, 1=CDS FASTA, 2=Text ASN.1*, 3=feature table* |
| (无) | `--tsv FILE` | ID/坐标映射表 |
| (无) | `--gff3 FILE` | GFF3 |
| (无) | `--id-mode numeric\|keep`, `--prefix STR` | numeric 默认: `SEQ000001_ORF2`; keep: 清洗保留 (`[A-Za-z0-9._-]` 白名单, 其余折叠为 `_`, 空名回退编号, 撞名加后缀) |
| (无) | `--threads N` | 默认全部核 |

*ASN.1 与 feature table 为近似格式, M4 对拍后逐字节校准。

## TSV 映射表

```
orf_id  orig_id  desc  strand  frame  start  end  nt_len  aa_len  start_codon  stop_codon  status
```

- `orig_id`/`desc` 中制表符与换行替换为空格 (其余字符原样保留, 包括非 ASCII)
- `status`: `complete` / `3partial` (有起始无终止)
- 坐标: 1-based 闭区间, 相对原始全长序列; minus 链 ORF 的 start < end, 方向看 strand 列
- circular 环绕 ORF 用线性化坐标 (end 可能 > 序列长度)

## 语义口径 (对拍基准, 详见 ../DECISIONS.md)

- ORF = 起始密码子到下游第一个终止密码子; 同框内每个起始密码子各成一条 ORF (嵌套照报, 与官方默认一致)
- 起始密码子在蛋白中译为 M (无论 ATG/GTG/...)
- ml 过滤按链方向全长 (含终止密码子) 计
- 输出顺序 = 输入序列顺序, 序列内按 (框, 链上位置); 与线程数无关

## 性能 (128 核 Xeon, 5 万条 / 150 Mbp / 4.55M ORF / 1 GB 输出)

| 阶段 | 1 线程 | 128 线程 |
|------|--------|----------|
| 扫描+翻译+格式化 | 7.03 s | **0.77 s** |
| 全流程 (输出到云盘) | 11.7 s | 5.3 s (受磁盘带宽 ~200 MB/s 限制) |
| 全流程 (输出到 tmpfs) | — | **0.96 s** |

结果跨线程数完全一致 (md5 相同)。

## 内存

序列、坐标与格式化后的输出驻留内存 (基准集峰值 ~2.3 GB); 大输入请按 输入+约 2x 输出量 估算。
