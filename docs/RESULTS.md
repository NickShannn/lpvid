# lpvid 结果表解读

Version: 0.2

每次运行保留各阶段中间产物，并将主要结果复制到 `results/`。
`results/main_results_report_zh.html` 与 `results/main_results_report_en.html`
用于浏览，`results/tables/` 中的 TSV 用于审核和二次分析，二者数据源一致。
`results/README.md` 是按本次实际交付文件生成的双语索引；序列位于
`results/sequences/`，Newick 位于 `results/trees/`。

## 1. 通用约定

- `query_id`：本次运行内部稳定 ID；原始 ORF、contig 和描述可通过映射表回溯。
- 多样本联合键是 `(sample, query_id)` 或 `(sample, contig_id)`，不能只靠 ORF
  编号。
- 空字段表示“不适用”或“没有该结果”，不自动等价于概率 0。
- `main` 是正式主分析集；`qa_discard` 是筛选丢弃序列的固定种子抽样审计集，
  不能混入主结果计数。
- HTML 报告用于展示，TSV 是精确数值的权威来源；静态图保留在阶段报告目录。

## 2. prefilter_calls.tsv

| 列 | 含义 |
|---|---|
| `sample` | 样本名 |
| `query_id` | 内部蛋白 ID |
| `aa_len` | 蛋白长度，单位 aa |
| `prefilter_call` | `KEEP_viral`、`KEEP_unknown` 或 `DISCARD` |
| `analysis_set` | `main` 或 `qa_discard` |

生产默认只让 ≥150 aa 的蛋白进入主分析。`KEEP_viral` 表示保留为病毒相关候选，
`KEEP_unknown` 表示无明确归属但保守保留，`DISCARD` 表示不进入主分析。被丢弃
序列中有固定比例进入 `qa_discard` 审计集，单独分析，不混入主结果。

## 3. orf_predictions.tsv

每行一个进入分析的蛋白 ORF。

| 列 | 含义 |
|---|---|
| `query_id` | 内部 ORF 蛋白 ID |
| `rdrp_prob` | RdRP 识别概率，范围 0 到 1 |
| `rdrp_positive` | 是否达到默认阈值 0.5 |
| `rdrp_status` | `ok`、`missing` 或 `excluded_input_validation` |
| `call` | `positive` 或 `negative` |
| `analysis_set` | `main` 或 `qa_discard` |

`missing` 表示没有取得识别结果；`excluded_input_validation` 表示该蛋白被输入
校验排除，不等同于生物阴性。概率是模型输出，不是实验确认。不同运行之间可以
比较同一字段的相对高低，但不应把概率直接解释为感染性或病原性。

## 4. contig_level_calls.tsv

| 列 | 含义 |
|---|---|
| `contig_id` | 原始核酸 contig ID |
| `n_orf_model_input` | 该 contig 进入分析的 ORF 数 |
| `rdrp_max_prob` | 该 contig 内最高的 RdRP 识别概率 |
| `call` | contig 级检出结论，`positive` 或 `negative` |

contig 上任一 ORF 为阳性，则该 contig 为阳性。蛋白输入没有 contig 映射，因此
该表为空是预期行为。

## 5. representative_map.tsv

| 列 | 含义 |
|---|---|
| `cluster_id` | 非冗余分组 ID |
| `query_id` | 组内成员原始内部 ID |
| `representative_id` | 该成员最终归属的代表 ID |
| `aa_len` | 成员长度，aa |
| `motif_class_count` | F/A/B/C 中至少命中一次的类别数，0 到 4 |
| `motif_classes` | 实际命中的 motif 类别 |
| `is_representative` | 1 表示最终代表，0 表示其他成员 |
| `exact_duplicate_of` | 精确重复对应的 canonical ID；非精确重复为空 |

代表选择是确定性的。只有代表进入超群分类和系统发育放置，所有成员仍保留在
该表中。

## 6. motifs.tsv

| 列 | 含义 |
|---|---|
| `query_id` | 阳性 ORF 内部 ID |
| `aa_len` | 蛋白长度，aa |
| `F` / `A` / `B` / `C` | 各 motif 的命中位置；多个位置用分号分隔，空为未命中 |
| `motif_class_count` | 至少命中一次的 motif 类别数 |
| `motif_classes` | 命中类别集合 |

位置为 0-based。该表是 motif 正则扫描结果，不应单独当作 RdRP 身份或完整性的
证明。

## 7. supergroup_assignments.tsv

| 列 | 含义 |
|---|---|
| `query_id` | 代表序列 ID |
| `supergroup` | 180 超群体系中的名称（如 Picorna、Bunya-Arena、Supergroup012）；未命中为空 |
| `status` | 分类与骨架树状态，见下表 |
| `classification_method` | `reference match` 或 `classifier classification` |
| `reference_id` | 参考序列 ID |
| `reference_identity` | 参考序列一致性，百分比 |
| `alignment_e_value` | 比对 E 值 |
| `alignment_score` | 比对分数 |
| `qcovhsp` / `scovhsp` | query / subject 的 HSP 覆盖度，百分比 |

`status` 取值：

| 值 | 含义 |
|---|---|
| `assigned` | 已分类，且有可用骨架树 |
| `classified_no_tree` | 已分类，但当前没有骨架树 |
| `classified_pending_merge` | 已分类，等待合并建树 |
| `classified_pending_rename` | 已分类，等待名称规范化 |
| `unassigned` | 未获得超群分类 |

参考序列一致性、比对 E 值和比对分数不是分类概率。`classified_no_tree` 序列已分类，
但因没有骨架树，暂不进行系统发育放置。辅助分类结果保留在阶段目录，不属于本交付
目录的核心结果表。

## 8. placements.tsv

每行是一条成功放置的 query。近邻列按 `nn1` 到 `nn5` 重复：

| 列 | 含义 |
|---|---|
| `nn1` 到 `nn5` | 最近的骨架树序列 ID |
| `nn1_dist` 到 `nn5_dist` | LucaPhylo 预测距离 |
| `method` | `lucaphylo_placement` |

这里只保留 LucaPhylo 实际输出的近邻和距离；没有额外的枝级置信度或衍生权重
列。可视化中的挂枝是最近邻近似。

## 8b. placement_annotations.tsv

`results/tables/placement_annotations.tsv` 是系统发育浏览器同一份数据的表格
形式，不进浏览器也能直接取用：

| 列 | 含义 |
|---|---|
| `query_id` | 本次运行内部序列 ID |
| `original_id` | fastorf 原始 ORF ID，可在原始 FASTA 中检索 |
| `supergroup` | 所属超群 |
| `nn1` 到 `nn5` | 最近 5 个骨架树序列 ID |
| `nn1_dist` 到 `nn5_dist` | LucaPhylo 预测距离 |

## 8c. trees/*.nwk

`results/trees/<超群>.nwk` 是本次可视化使用的骨架树 Newick，与浏览器完全一致，
可直接导入 iTOL、FigTree 等工具；浏览器内的“导出”也可输出 Newick 与注释表。

## 8d. risk_warning.tsv 与 neighbor_evidence.tsv

`results/tables/risk_warning.tsv` 每条放置代表一行，展示类群归属、原始 ORF 与
contig 映射、以及名录提示。它不构成新序列的风险分级或病原确证。缺失距离与
未知宿主标记导出为空，不填 0。

`results/tables/neighbor_evidence.tsv` 是长表，每个 query 的每个真实
`nn1` 到 `nn5` 近邻一行，空白近邻列不重排 rank。列为 `query_id`、`rank`、
`reference`、`distance`、`family`、`genus`、`human`、`vertebrate`、`plant`；
taxonomy 缺失时分类和宿主字段为空。

## 9. abundance.tsv

| 列 | 含义 |
|---|---|
| `contig_id` | RdRP 阳性源 contig |
| `length` | contig 长度，bp |
| `mapped_reads` | primary mapped reads 数 |
| `RPM` | `mapped_reads / total_clean_reads × 1e6` |
| `mean_depth` | 平均测序深度 |
| `coverage_breadth_pct` | ≥1× 覆盖的 contig 比例，百分比 |
| `total_clean_reads` | 完成已启用清理步骤后的分母 |
| `count_unit` | 当前为 `primary_mapped_reads` |

双端分母是 R1+R2 reads，分子也按 reads 计数；单端只计算 R1。RPM 是文库内
相对丰度，不是 copies/mL 或绝对载量。

## 10. supergroup_abundance.tsv

| 列 | 含义 |
|---|---|
| `supergroup` | 唯一归属超群 |
| `n_contigs` | 纳入汇总的 contig 数 |
| `mapped_reads` | 这些 contig 的 primary mapped reads 合计 |
| `RPM` | 合计 mapped reads 按 clean reads 重新计算的 RPM |
| `mean_depth_weighted` | 按 contig 长度加权的平均深度 |
| `coverage_breadth_pct` | 按 contig 长度加权的覆盖广度 |
| `total_clean_reads` | 与 contig 表相同的分母 |

多重超群 contig 不重复计入超群总量，在阶段目录的
`abundance_contig_assignments.tsv` 中标记为 `ambiguous`。

## 10b. 补充明细表

`results/tables/` 还可能包含以下运行明细：

| 文件 | 内容 |
|---|---|
| `supergroup_details.tsv` | 超群分类明细，每条代表一行 |
| `placement_details.tsv` | 系统发育放置明细，每条代表一行 |
| `stage_timings.tsv` | 各阶段状态与耗时 |
| `sequence_funnel.tsv` | 各步骤序列数量与类型 |
| `stage_metrics.tsv` | 各阶段输入、输出与耗时 |

这些表用于运行审计与复现，不是主结果的替代品。它们的实际列集合以本次交付的
`results/README.md` 和文件本身为准。

## 11. summary 与 HTML

- `results/tables/run_summary.json`：交付的结构化运行汇总；`risk_rows` 为对象
  列表，每个 `neighbors` 为嵌套证据对象列表，缺失距离和宿主标记为 `null`。
- `results/tables/delivery_manifest.json`：本轮交付成功哨兵。只有报告、JSON
  和 README 全部成功后才原子写入 `status: "complete"`；`files` 为相对
  `results/` 的文件清单。任一中途失败都不会保留完成标记。
- `step*_report/summary.md`：便于阅读的数量、链接、警告和阶段状态中间稿。
- `step*_report/summary.json`：与阶段汇总对应的机器可读中间产物。
- `results/main_results_report_zh.html` / `main_results_report_en.html`：独立、
  自包含的中英文主报告。
- `results/phylogenetic_explorer.html`：独立的系统发育放置浏览器。
- `results/README.md`：按本次实际交付文件生成的双语索引。

交付根目录不复制 `summary.md` 或静态图。三个 HTML 文件均不需要启动 Web 服务。
