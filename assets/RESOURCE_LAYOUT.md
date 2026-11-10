# lpvid 外部资源目录约定

Version: 0.2

源码与大型资源分离发布。推荐将资源放在仓库内 `resources/` 目录，
`config/local.example.yaml` 默认采用该相对布局；也可以放在任意位置，然后在
`config/local.yaml` 中指向真实路径。文件级清单见
[`resources/RESOURCE_MANIFEST.tsv`](../resources/RESOURCE_MANIFEST.tsv)，双语
说明见 [`resources/README.md`](../resources/README.md)。

推荐结构：

```text
resources/
  phylo/
    backbone_dir/<超群>/{*.nwk,*.fas}
    adapter_dir/<超群>/{adapter_config.json,adapter_model.safetensors,...}
    taxonomy/<超群>.csv
  rdrp/
    rdrp_reference.dmnd
    rdrp_reference.fas
    rdrp_reference.csv
  models/
    rdrp_detector/
    lucaphylo/
    classifier/
      pipeline/
      encoder/
      model/
  blastdb/refseq_protein/
    refseq_protein.dmnd
    TAXID_viruses.txt
    TAXID_cellular.txt
  rrna/
    silva_ssu.*
    silva_lsu.*
```

## 资源组与用途

| 资源组 | 内容与用途 |
|---|---|
| `phylo/` | 精选系统发育骨架树、adapter 与 taxonomy；用于超群分类后的放置与可视化 |
| `rdrp/` | RdRP 参考库（`.dmnd`、`.fas`、`.csv`）；用于超群分类与参考匹配 |
| `models/rdrp_detector/` | RdRP 检测模型权重与运行文件 |
| `models/lucaphylo/` | LucaPhylo 系统发育模型文件 |
| `models/classifier/` | 分类模型的 `pipeline`、`encoder` 与 `model` 文件 |
| `blastdb/refseq_protein/` | RefSeq 蛋白参考数据库与 taxid 白名单；用于蛋白筛选 |
| `rrna/` | SILVA SSU 与 LSU bowtie2 索引；用于 reads 路线的 rRNA 去除 |

报告资产、树浏览器资产、名录资产与超群映射随 git 仓库提供，不需要单独下载。

## 配置与本机路径

- 清单中的安装路径均相对仓库根目录。
- `config/local.example.yaml` 使用相同的相对路径；复制为
  `config/local.yaml` 后按本机情况调整。
- 请在仓库根目录运行 `lpvid`，使相对路径解析到仓库内资源；也可以在
  `config/local.yaml` 中使用绝对路径。
- `config/local.yaml` 是本机私有配置，不纳入 git。

## 不随源码发布的内容

- 模型权重与大型参考数据库；
- 系统发育骨架树资源包；
- rRNA 索引；
- 运行输出、缓存和本机私有配置；
- 原始验证数据与历史运行目录。

代码许可不等同于模型或数据许可。再分发前请确认各资源组的使用条款。

## 校验与准备

```bash
lpvid doctor
lpvid prepare
```

`doctor` 逐项检查环境、工具、模型与数据库路径；`prepare` 只校验本地骨架树
资产并构建 `assets/` 下的预构建参考资产，不会从互联网下载任何内容。

大型外部下载项在清单中标记为 `verify_after_download`。下载完成后请计算并记录
SHA-256，再投入业务使用。
