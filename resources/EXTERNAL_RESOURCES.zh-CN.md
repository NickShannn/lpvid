# lpvid 0.2 外部资源

资源记录：https://doi.org/10.5281/zenodo.24122793（概念 DOI，始终指向最新版本）
DOI：`10.5281/zenodo.24122793`

资源记录包含五个压缩包。请全部下载到临时目录，先校验 SHA-256，再在仓库根目录
解压。每个压缩包都解压到 `resources/`，文件会落到 `config/default.yaml` 与
`config/local.example.yaml` 预期的位置。

## 压缩包

| 压缩包 | 内容 | 安装路径 | 字节数 | SHA-256 |
|---|---|---|---|---|
| `lpvid-0.2-phylo-adapters.tar.gz` | 24 个 LoRA 放置 adapter（168 个文件） | `resources/phylo/adapter_dir/` | 503,943,636 | `0f63127818bf11ecc566badd679d7d38bc63ec63ae770b13dfa3622e0d9f9535` |
| `lpvid-0.2-lucaphylo-base.tar.gz` | LucaPhylo 基础 checkpoint，仅推理文件 | `resources/models/lucaphylo/checkpoint-81000/` | 2,423,201,239 | `d596d006f6ca005af70476186c832ba3c6d6c44eb6c7abab0c751676aa2d0aa9` |
| `lpvid-0.2-classifier-encoder.tar.gz` | 分类编码器 checkpoint，仅推理文件 | `resources/models/classifier/encoder/checkpoint-18000/` | 2,432,855,340 | `4f43cf2cc5ed94c51bf8d1adc129b8a0ea5a728bc4d4a7d267e2279b4c9b4158` |
| `lpvid-0.2-rdrp-detector.tar.gz` | 检测代码、运行配置与一个 checkpoint | `resources/models/rdrp_detector/` | 75,154,141 | `a5e45d09e8275bf79500b93a31bc913209d30d5e82765b275efaae59d435f1ed` |
| `lpvid-0.2-silva-rrna.tar.gz` | 12 个 SILVA SSU/LSU bowtie2 索引文件 | `resources/rrna/` | 1,430,723,048 | `3b05c54a664c4ad9b6a52ecf2f9e1c48c8e3bcb5815d9e01a839e1f1e1fc2e0f` |

模型压缩包只包含推理所需文件（`config.json`、`model.safetensors`、
`tokenizer_config.json`、`special_tokens_map.json`、`vocab.txt`），不包含训练
优化器状态、随机数状态或训练日志。

adapter 配置中的基础模型路径已改为仓库相对路径
`resources/models/lucaphylo/checkpoint-81000`，可在任意安装目录使用；
`lpvid prepare` 会解析该路径并拒绝真实不一致。

## RefSeq 蛋白数据库

RefSeq 蛋白 DIAMOND 数据库
（`resources/blastdb/refseq_protein/refseq_protein.dmnd`）约 214 GB，超过资源
记录体积上限，需单独分发；两个 `TAXID_*.txt` 白名单已随仓库内置。当前下载
渠道见发布说明。

## 安装

```bash
cd /path/to/lpvid
sha256sum -c SHA256SUMS.txt
for archive in lpvid-0.2-*.tar.gz; do tar xzf "$archive"; done
lpvid doctor
lpvid prepare
```

`SHA256SUMS.txt` 与资源记录一同发布。仓库内置资源的逐文件校验仍在
[`RESOURCE_MANIFEST.tsv`](RESOURCE_MANIFEST.tsv)。

## 条款

源码采用 Apache License 2.0。模型权重、adapter、骨架树、taxonomy 表与数据库
属于数据资产，其再分发条款以资源记录为准。详见
[`../THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md)。
