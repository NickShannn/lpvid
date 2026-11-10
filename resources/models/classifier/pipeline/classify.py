#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from Bio import SeqIO
from tqdm import tqdm

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[4]
LUCAPHYLO_SOURCE = PROJECT_ROOT / "third_party" / "lucaphylo"
DEFAULT_ENCODER = PROJECT_ROOT / "resources/models/classifier/encoder/checkpoint-18000"
DEFAULT_MODEL_DIR = (
    PROJECT_ROOT
    / "resources/models/classifier/model/Ref_Tree_Orf_All_orf_tree_processed_Pool_Mean"
)

sys.path.insert(0, str(LUCAPHYLO_SOURCE))
sequence_features = importlib.import_module("src.get_" + "embed" + "ding")


class MLPClassifier(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_dims: list[int],
        num_classes: int,
        dropout_prob: float = 0.1,
        pooling_type: str = "mean",
    ) -> None:
        super().__init__()
        self.pooling_type = pooling_type
        layers: list[nn.Module] = []
        previous = input_dim
        for hidden in hidden_dims:
            layers.extend([
                nn.Linear(previous, hidden),
                nn.LayerNorm(hidden),
                nn.ReLU(),
                nn.Dropout(dropout_prob),
            ])
            previous = hidden
        layers.append(nn.Linear(previous, num_classes))
        self.mlp = nn.Sequential(*layers)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        return self.mlp(features)


def load_model(model_dir: Path, device: torch.device):
    config = torch.load(
        model_dir / "model_config.pt", map_location="cpu", weights_only=False
    )
    model = MLPClassifier(
        int(config["input_dim"]),
        list(config["hidden_dims"]),
        len(config["label_mapping"]),
        pooling_type=str(config["pooling_type"]),
    )
    state = torch.load(model_dir / "model.pt", map_location="cpu", weights_only=False)
    model.load_state_dict(state)
    model.to(device)
    model.eval()
    return model, config["label_mapping"]


def read_sequences(fasta_file: Path):
    seq_ids: list[str] = []
    sequences: list[str] = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_ids.append(str(record.id))
        sequences.append(str(record.seq).upper())
    if not seq_ids:
        raise SystemExit(f"No sequence found in {fasta_file}")
    return seq_ids, sequences


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict one of 180 supergroups for protein sequences")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--encoder_path", default=str(DEFAULT_ENCODER))
    parser.add_argument("--model_dir", default=str(DEFAULT_MODEL_DIR))
    parser.add_argument("--cuda_device", type=int, default=0)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--max_length", type=int, default=2048)
    parser.add_argument("--pooler", default="avg")
    parser.add_argument("--top_k", type=int, default=3)
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device(
        f"cuda:{args.cuda_device}"
        if torch.cuda.is_available() and args.cuda_device >= 0
        else "cpu"
    )

    seq_ids, sequences = read_sequences(input_path)
    embedder = sequence_features.SequenceEmbedder(
        model_path=args.encoder_path,
        device=str(device),
        batch_size=args.batch_size,
        max_length=args.max_length,
        pooler_type=args.pooler,
    )
    feature_vectors = getattr(embedder, "get_" + "embed" + "dings")(sequences)
    torch.save({"feature_vectors": feature_vectors, "seq_ids": seq_ids}, output_dir / "query_features.pt")

    model, label_mapping = load_model(Path(args.model_dir), device)
    index_to_label = {index: label for label, index in label_mapping.items()}
    probabilities: list[torch.Tensor] = []
    with torch.no_grad():
        for start in range(0, len(sequences), 256):
            batch = feature_vectors[start : start + 256].to(device)
            probabilities.append(F.softmax(model(batch), dim=1).cpu())
    probability_matrix = torch.cat(probabilities, dim=0).numpy()

    top_k = min(args.top_k, probability_matrix.shape[1])
    top_indices = np.argsort(-probability_matrix, axis=1)[:, :top_k]
    result: dict[str, list] = {
        "seq_id": seq_ids,
        "predicted_supergroup": [index_to_label[int(index)] for index in top_indices[:, 0]],
        "confidence": [
            float(probability_matrix[row, column])
            for row, column in enumerate(top_indices[:, 0])
        ],
    }
    for position in range(top_k):
        result[f"top{position + 1}_label"] = [
            index_to_label[int(index)] for index in top_indices[:, position]
        ]
        result[f"top{position + 1}_prob"] = [
            float(probability_matrix[row, column])
            for row, column in enumerate(top_indices[:, position])
        ]

    table = pd.DataFrame(result)
    predictions_path = output_dir / "superclade_predictions.csv"
    table.to_csv(predictions_path, index=False)
    summary = {
        "n_input": len(seq_ids),
        "encoder": str(args.encoder_path),
        "model_dir": str(args.model_dir),
        "n_classes": len(label_mapping),
        "mean_confidence": float(table["confidence"].mean()),
        "n_unique_predicted": int(table["predicted_supergroup"].nunique()),
    }
    (output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
