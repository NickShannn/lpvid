# 超群体系说明

Version: 0.2

## 1. 超群命名

当前分类使用 180 超群体系。早期体系使用 24 个超群编号；二者的对应关系记录在
仓库根目录的 `supergroups_mapping.csv`，以其中 `LucaPhylo_SGs` 列的名称为准。
结果表和报告中的 `supergroup` 字段统一使用当前名称，例如 Picorna、
Bunya-Arena、Supergroup012。

映射表保留历史编号是为了兼容既有运行和资源目录，不代表当前分类仍使用旧体系。

## 2. 分类结果与状态

每个代表序列在 `results/tables/supergroup_assignments.tsv` 中有一行，主要字段
包括 `supergroup`、`status`、`classification_method`、`reference_id`、
`reference_identity`、`alignment_e_value`、`alignment_score`、`qcovhsp` 与
`scovhsp`。

`classification_method` 取值：

- `reference match`：代表序列与参考库序列匹配后获得超群名称；
- `classifier`：由分类模型给出超群预测。

`status` 取值：

| 值 | 含义 |
|---|---|
| `assigned` | 已分类，且有可用骨架树，可进入系统发育放置 |
| `classified_no_tree` | 已分类，但当前没有骨架树，暂不放置 |
| `classified_pending_merge` | 已分类，等待合并建树 |
| `classified_pending_rename` | 已分类，等待名称规范化 |
| `unassigned` | 未获得超群分类 |

参考序列一致性、比对 E 值和比对分数是证据指标，不是分类概率。低置信度分类结果
建议结合其他证据人工复核。

## 3. 骨架树覆盖范围

当前带注释的骨架树覆盖 180 超群体系中的一个子集，共 24 棵。能够被分类但暂无
骨架树的超群会在结果中标记为 `classified_no_tree`，不进行系统发育放置。

部分超群需要合并建树，例如将相近的类群合并为单棵骨架树。合并建树后需人工
检查原有 family 与 genus 是否仍为单系。

## 4. 合并建树参考命令

以下工具仅在建树时使用，不随统一环境 `environment.yml` 安装。

```bash
input_file=$1
base_name=$(basename "$input_file" .fas)

# 多序列比对：MAFFT
mafft --maxiterate 1000 --localpair --thread -1 \
    "$input_file" > "${base_name}_mafft.fas"

# 序列修剪：trimAl
trimal -in "${base_name}_mafft.fas" \
    -out "${base_name}_mafft_trimal_cons5_gt0.9.fas" \
    -gt 0.9 -cons 5 -keepheader

# 建树：IQ-TREE
iqtree2 -s "${base_name}_mafft_trimal_cons5_gt0.9.fas" \
    -m LG --alrt 1000 -B 1000 -T AUTO \
    --prefix "${base_name}_mafft_trimal_cons5_gt0.9_LG"
```

资源目录约定见 [assets/RESOURCE_LAYOUT.md](assets/RESOURCE_LAYOUT.md)。
