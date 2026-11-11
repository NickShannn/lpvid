# External resources for lpvid 0.2

Resource record: https://doi.org/10.5281/zenodo.24122793 (concept DOI; always
resolves to the latest version)
DOI: `10.5281/zenodo.24122793`

The resource record contains five archives. Download all five into a temporary
directory, verify their SHA-256 values, and extract them from the repository
root. Every archive extracts into `resources/`, so the files land at the paths
expected by `config/default.yaml` and `config/local.example.yaml`.

## Archives

| Archive | Contents | Installed path | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| `lpvid-0.2-phylo-adapters.tar.gz` | 24 LoRA placement adapters (168 files) | `resources/phylo/adapter_dir/` | 503,943,636 | `0f63127818bf11ecc566badd679d7d38bc63ec63ae770b13dfa3622e0d9f9535` |
| `lpvid-0.2-lucaphylo-base.tar.gz` | LucaPhylo base checkpoint, inference files only | `resources/models/lucaphylo/checkpoint-81000/` | 2,423,201,239 | `d596d006f6ca005af70476186c832ba3c6d6c44eb6c7abab0c751676aa2d0aa9` |
| `lpvid-0.2-classifier-encoder.tar.gz` | Classifier encoder checkpoint, inference files only | `resources/models/classifier/encoder/checkpoint-18000/` | 2,432,855,340 | `4f43cf2cc5ed94c51bf8d1adc129b8a0ea5a728bc4d4a7d267e2279b4c9b4158` |
| `lpvid-0.2-rdrp-detector.tar.gz` | Detection code, runtime configuration and one checkpoint | `resources/models/rdrp_detector/` | 75,154,141 | `a5e45d09e8275bf79500b93a31bc913209d30d5e82765b275efaae59d435f1ed` |
| `lpvid-0.2-silva-rrna.tar.gz` | 12 SILVA SSU/LSU bowtie2 index files | `resources/rrna/` | 1,430,723,048 | `3b05c54a664c4ad9b6a52ecf2f9e1c48c8e3bcb5815d9e01a839e1f1e1fc2e0f` |

The model archives contain the files required for inference (`config.json`,
`model.safetensors`, `tokenizer_config.json`, `special_tokens_map.json`,
`vocab.txt`). Training optimizer states, RNG states and trainer logs are not
included.

The adapter configurations use the repository-relative base-model path
`resources/models/lucaphylo/checkpoint-81000`, so they work from any installation
directory. `lpvid prepare` resolves the path and rejects real mismatches.

## RefSeq protein database

The RefSeq protein DIAMOND database (`resources/blastdb/refseq_protein/refseq_protein.dmnd`)
is about 214 GB and exceeds the resource-record size limit. It is distributed
separately; the two bundled `TAXID_*.txt` whitelists accompany it. See the
release announcement for the current download channel.

## Install

```bash
cd /path/to/lpvid
sha256sum -c SHA256SUMS.txt
for archive in lpvid-0.2-*.tar.gz; do tar xzf "$archive"; done
lpvid doctor
lpvid prepare
```

`SHA256SUMS.txt` is published in the same resource record. File-level checks for
the bundled resources remain in [`RESOURCE_MANIFEST.tsv`](RESOURCE_MANIFEST.tsv).

## Terms

Source code is distributed under the Apache License 2.0. Model weights,
adapters, backbone trees, taxonomy tables and databases are data assets; their
redistribution terms are stated on the resource record. See
[`../THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md).
