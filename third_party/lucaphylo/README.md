# LucaPhylo Training Release

This directory contains a GitHub-ready release for the four-stage LucaPhylo pipeline.

The goal of this refactor is the minimum open-source version of the training code:

- keep the core training logic used in the paper,
- remove lab-specific absolute paths and backup logic,
- unify argument names across stages,
- leave only the components needed for reproduction and extension.

The original research scripts are preserved in [`source_raw/`](../source_raw). The public release lives in [`release/`](.).

This release is intended for basic reproduction of the stage-specific checkpoints described in the manuscript. It preserves the main objectives, data interfaces, and training flow for stages 1 to 4, but it is not meant to be a byte-for-byte replay of every internal experiment script.

## Directory Layout

```text
release/
|-- README.md
|-- infer_stage_4.py
|-- requirements.txt
|-- train_stage_1.py
|-- train_stage_2_3.py
|-- train_stage_4.py
|-- scripts/
|   |-- run_stage_1.sh
|   |-- run_stage_2.sh
|   |-- run_stage_3.sh
|   |-- run_stage_4_infer.sh
|   `-- run_stage_4_train.sh
`-- src/
    |-- __init__.py
    |-- dataset.py
    |-- model.py
    |-- stage4_data.py
    |-- stage4_model.py
    `-- utils.py
```

## Stage Overview

1. `train_stage_1.py`
   Continued pretraining with masked language modeling on protein sequences. Stage 1 keeps the base ESM positional embedding configuration unchanged. When `max_seq_length` exceeds the checkpoint-supported context window, inputs are left-truncated to the supported length and keep the right-hand sequence region.
2. `train_stage_2_3.py --training_stage 2`
   Label-aware contrastive learning. By default, the public launcher now follows the original paper-training implementation more closely: each stage-2 training item dynamically samples multiple labels, draws one positive pair per label, and optimizes a class-balanced contrastive batch. You can disable this behavior with `--stage2_legacy_mode false` to recover the simplified row-wise release implementation.
3. `train_stage_2_3.py --training_stage 3`
   Paired contrastive learning between anchor and paired domain sequences, with orthogonality and separation regularization.
4. `train_stage_4.py`
   Distance regression on phylogenetic pairs using the stage-2 or stage-3 checkpoint as the base encoder. Stage 4 applies LoRA to the ESM attention projections and keeps the hyperbolic distance head trainable. Public checkpoints save only the LoRA adapter plus the hyperbolic `scale` module.

## Expected Input Formats

### Stage 1

- `.txt`: one protein sequence per line
- `.csv` or `.json`: one column containing sequences
- if the sequence column is not named `text`, `sequence`, or `seq`, pass `--text_column`

### Stage 2

CSV with at least:

```text
seq,label
SEQUENCE_A,Label_1
SEQUENCE_B,Label_1
SEQUENCE_C,Label_2
```

Default columns:

- `seq`: anchor sequence
- `label`: class / supergroup label

Default stage-2 launcher behavior:

- `release/scripts/run_stage_2.sh` uses the legacy class-balanced stage-2 sampler by default
- each sampled stage-2 training item contains 4 positive pairs drawn from 4 distinct labels
- the launcher keeps `per_device_train_batch_size=1` and `gradient_accumulation_steps=2`, which matches the original effective batch construction more closely
- training logs keep the original Trainer `epoch` and also record `conventional_epoch_estimate`, an expected sequence-coverage epoch for easier comparison with standard training runs

### Stage 3

CSV with at least:

```text
seq,domain
POLYPROTEIN_SEQ_1,DOMAIN_SEQ_1
POLYPROTEIN_SEQ_2,DOMAIN_SEQ_2
```

Default columns:

- `seq`: anchor sequence
- `domain`: paired positive sequence

### Stage 4 Train

- `--tree_path`: backbone tree in Newick format
- `--fasta_path`: FASTA whose record IDs match the tree leaves
- optional `--split_file`:
  - either columns `train,valid`, `train,test`, or `train,valid,test`
  - or columns `id,split`

Stage-4 training expands the tree into supervised sequence pairs:

- train: `train x train`
- valid: `valid x train` plus `valid x valid`, if provided
- test: `test x train`, `test x valid`, and `test x test`, if provided

For same-split pairs, the public release keeps one copy per unordered pair and defaults to excluding diagonal self-pairs, which matches the original stage-4 training script more closely.

The default `release/scripts/run_stage_4_train.sh` launcher enables training and validation only. Add `--do_predict` manually when you also provide a test split and want test-set metrics during the same run. If `--do_predict` is not enabled, any test IDs present in the split file are ignored.

#### Fast pipeline and training controls (2026-10)

- `--data_pipeline indexed|legacy` (default `indexed`): `indexed` pre-tokenizes each unique sequence once and feeds pairs by ID reference — removes the O(n^2) per-epoch tokenization and the O(n^2) string duplication in the pair DataFrames. `legacy` reproduces the historical on-the-fly behaviour byte-for-byte (A/B baseline).
- `--pair_budget_ratio` (default `1.0`): fraction of train pairs sampled per epoch, reshuffled every epoch with length bucketing (`src/stage4_sampler.py`). Validated as loss-neutral down to ~0.25 on pair-subsampling experiments.
- `--eval_pair_budget` (default `0` = full): cap validation pairs per evaluation via distance-stratified subsampling — full valid pairs cost ~0.13*n^2, rivaling a training epoch.
- `--split_seed` (default `-1`): decouples the taxa-split seed from the training seed (`-1` keeps the legacy single-seed behaviour). Multi-split runs use one seed per run.
- `--split_mode random|clade` (default `random`): `clade` splits whole clades instead of individual leaves — evaluation then measures generalization to unseen lineages.
- `--early_stopping_patience` (default `0` = off; launcher default `3`): requires `--load_best_model_at_end` and `--metric_for_best_model`.
- `--attn_implementation sdpa|eager` (default `sdpa`): ESM attention kernel; falls back to `eager` automatically if unavailable.

Recommended single-supergroup run (launcher env): `DATA_PIPELINE=indexed PAIR_BUDGET_RATIO=0.25 EVAL_PAIR_BUDGET=30000 EARLY_STOP_PATIENCE=3`.

Multi-split protocol (k split seeds on one GPU each, then mean ± std):

```bash
setsid nohup env SPLIT_SEEDS="42 43 44 45 46" GPU_IDS=0,1,2,3,4,5,6,7 \
  TREE_PATH=... FASTA_PATH=... MODEL_NAME_OR_PATH=... \
  bash scripts/run_stage_4_multisplit.sh > multisplit.log 2>&1 &
```

### Stage 4 Infer

- `--model_name_or_path`: stage-2 or stage-3 base checkpoint
- `--stage4_adapter_path`: stage-4 adapter directory saved by `train_stage_4.py`; a comma-separated list enables fold-ensemble inference (per-adapter distances are averaged)
- `--backbone_tree`: backbone tree in Newick format
- `--backbone_seqs`: backbone FASTA
- `--query_seqs`: query FASTA

## Quick Start

Install dependencies:

```bash
pip install -r release/requirements.txt
```

If `ete3` is installed, stage-4 uses it by default for faster backbone-tree distance parsing and falls back to Biopython otherwise.

Run the four stages:

```bash
bash release/scripts/run_stage_1.sh
bash release/scripts/run_stage_2.sh
bash release/scripts/run_stage_3.sh
bash release/scripts/run_stage_4_train.sh
bash release/scripts/run_stage_4_infer.sh
```

For stage 2, stage 3, and stage 4, the recommended starting point is the checkpoint produced by the previous stage. Stages 1 to 3 keep the positional embedding size from the loaded checkpoint and left-truncate overlength inputs to the supported context window.

For stage 2 specifically, the default launcher is tuned to be closer to the original paper-training behavior than the earlier simplified release refactor. If you prefer the simplified release behavior for ablations, pass `--stage2_legacy_mode false` and adjust the batch-size settings accordingly.

Each stage saves both model weights and tokenizer files into its output directory, so the next stage can load directly from the previous stage checkpoint folder.

Stage 4 is slightly different: when LoRA is enabled, the saved checkpoint contains the adapter weights and the hyperbolic `scale` module, while inference still loads the stage-2 or stage-3 checkpoint as the base model.

The shell launchers are templates. Override paths and hardware settings with environment variables:

```bash
GPU_IDS=0,1,2,3 \
NUM_GPUS=4 \
MODEL_NAME_OR_PATH=/path/to/stage1_checkpoint \
TRAIN_FILE=/path/to/stage2.csv \
OUTPUT_DIR=/path/to/output \
bash release/scripts/run_stage_2.sh

GPU_IDS=4,5,6,7 \
MODEL_NAME_OR_PATH=/path/to/stage3_checkpoint \
TREE_PATH=/path/to/backbone.nwk \
FASTA_PATH=/path/to/backbone.fasta \
OUTPUT_DIR=/path/to/stage4_adapter \
bash release/scripts/run_stage_4_train.sh

GPU_IDS=7 \
MODEL_NAME_OR_PATH=/path/to/stage3_checkpoint \
STAGE4_ADAPTER_PATH=/path/to/stage4_adapter \
BACKBONE_TREE=/path/to/backbone.nwk \
BACKBONE_SEQS=/path/to/backbone.fasta \
QUERY_SEQS=/path/to/query.fasta \
OUTPUT_DIR=/path/to/infer_output \
bash release/scripts/run_stage_4_infer.sh
```

For the stage-4 train and inference launchers, `GPU_IDS` sets `CUDA_VISIBLE_DEVICES`. When `GPU_IDS` is provided for training and `NUM_GPUS` is omitted, the launcher automatically uses the number of listed GPUs. For inference, keep `DEVICE=cuda:0` unless you intentionally want a different index within the visible devices.

Extra CLI flags can also be appended after the script call and will be forwarded to the Python entrypoint.

## What Was Intentionally Removed

To follow a minimum-release principle, this public version does not include:

- machine-specific absolute paths,
- one-off experiment backup scripts,
- checkpoint schedules tied to a private server workflow,
- commented-out branches that were not required for the main stage pipeline.

If you want, you can keep `source_raw/` and `stage4_train_infer/` private for archival purposes and publish only `release/` as the public training and inference code.
