# lpvid 资源说明

Version: 0.2

本目录记录 lpvid 运行所需的外部资源。源码与大型资源分开分发：git 仓库保持
轻量，模型与数据库放置在本目录（或任意位置），由 `config/local.yaml` 引用。

完整文件级清单见 [`RESOURCE_MANIFEST.tsv`](RESOURCE_MANIFEST.tsv)，各列含义：

| 列 | 含义 |
|---|---|
| `category` | 资源类别 |
| `item` | 文件名或目录条目 |
| `required` | `required`（必需）、`optional`（可选）或 `generated`（由流程生成） |
| `install_path` | 相对仓库根目录的安装路径 |
| `purpose` | 工作流中的用途 |
| `size_bytes` | 字节数；目录为内部文件合计 |
| `sha256` | 打包时已知的 SHA-256；超大型外部条目留空 |
| `source_hint` | 来源：仓库内置、资源包分发，或发布渠道 |
| `notes` | 获取或校验说明，例如 `verify_after_download` |

最终发布清单是文件数量、体量与校验和的权威来源。本文档不重复这些统计数字。

## 目录约定

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

清单中的安装路径一律相对仓库根目录，`config/local.example.yaml` 使用完全
一致的相对路径。若资源放在其他位置，请在 `config/local.yaml` 中使用绝对路径。

## 资源组

| 资源组 | 内容 | 分发方式 |
|---|---|---|
| 报告资产 | HTML 报告的离线图表与表格库 | git 仓库内置 |
| 浏览器资产 | 树浏览器的 PDF 与矢量导出库 | git 仓库内置 |
| 名录资产 | 报告提示使用的参考名录 | git 仓库内置 |
| 超群映射 | 历史编号到当前超群名称的映射 | git 仓库内置 |
| 系统发育骨架树 | 精选骨架比对、树文件与 taxonomy 表 | git 仓库内置 |
| 系统发育 adapter | 24 个超群的 LoRA 放置 adapter | Zenodo 资源记录 |
| RdRP 参考库 | 参考数据库、序列与标签 | git 仓库内置 |
| RdRP 检测模型 | 检测模型权重与运行文件 | Zenodo 资源记录 |
| 系统发育模型 | LucaPhylo 基础模型与分类编码器 | Zenodo 资源记录 |
| 分类模型 | 启动脚本、标签表与运行时模型文件 | git 仓库内置 |
| 分类编码器 | 分类流程使用的编码器 checkpoint | Zenodo 资源记录 |
| 蛋白参考数据库 | RefSeq 蛋白 DIAMOND 数据库 | 单独分发（体积大） |
| RefSeq taxid 白名单 | 粗筛使用的病毒与细胞 taxid 名单 | git 仓库内置 |
| rRNA 索引 | SILVA SSU 与 LSU bowtie2 索引 | Zenodo 资源记录 |
| 预构建资产 | `assets/` 下由流程生成的本地参考资产 | `lpvid prepare` 生成 |

## 校验原则

- git 仓库内置文件的 SHA-256 已写入清单。
- 可下载资源压缩包及其 SHA-256 见
  [`EXTERNAL_RESOURCES.zh-CN.md`](EXTERNAL_RESOURCES.zh-CN.md)。
- 压缩包内单个文件在清单中标记为 `verify_after_download`；以压缩包校验和作为
  完整性判据。
- 工作流通过 `lpvid doctor` 校验资源就绪状态；doctor 只检查配置路径是否存在，
  不会对大型数据库重新哈希。
- `lpvid prepare` 校验本地骨架树、taxonomy 表与 adapter，并构建 `assets/`
  下的预构建参考资产；它不会从互联网下载任何内容。

## 获取说明

- **仓库内置资源**：`resources/rdrp/`、`resources/phylo/backbone_dir/`、
  `resources/phylo/taxonomy/`、分类启动脚本与标签表、分类运行时模型文件，
  以及两个 RefSeq taxid 白名单均已随 git 仓库提供，SHA-256 写入清单。
- **系统发育 adapter**：`resources/phylo/adapter_dir/` 包含 24 个 LoRA 放置
  adapter，整组约 545 MiB。
- **外部模型资源**：`models/rdrp_detector`、`models/lucaphylo` 与
  `models/classifier/encoder` 通过 Zenodo 资源记录分发；压缩包名称、大小与
  校验和见
  [`EXTERNAL_RESOURCES.zh-CN.md`](EXTERNAL_RESOURCES.zh-CN.md)。
- **蛋白参考数据库**：由 NCBI RefSeq 蛋白序列构建；两个 `TAXID_*.txt`
  白名单已内置，DIAMOND 数据库需单独下载。
- **SILVA SSU 与 LSU**：由 SILVA 核糖体 RNA 数据库构建的 bowtie2 索引，索引
  文件共享前缀，例如 `resources/rrna/silva_ssu`。
- **预构建资产**：无需下载；`lpvid prepare` 会基于系统发育骨架包重新生成。

资源就位后，在仓库根目录执行：

```bash
lpvid doctor
lpvid prepare
```

## English version

See [README.md](README.md).
