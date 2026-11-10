import json
import os
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
import torch
from Bio import Phylo
from torch.utils.data import Dataset


def parse_fasta(path: str) -> "OrderedDict[str, str]":
    sequences: "OrderedDict[str, str]" = OrderedDict()
    current_id = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(chunks).upper()
                current_id = line[1:].strip().split()[0]
                chunks = []
                continue

            chunks.append("".join(line.split()))

    if current_id is not None:
        sequences[current_id] = "".join(chunks).upper()

    if not sequences:
        raise ValueError(f"No FASTA records found in `{path}`.")

    return sequences


def _load_tree_distance_matrix_biopython(tree_path: str):
    tree = Phylo.read(tree_path, "newick")
    taxa = [clade.name for clade in tree.get_terminals()]
    if any(name is None for name in taxa):
        raise ValueError("All leaf nodes in the backbone tree must have names.")

    distance_matrix = np.zeros((len(taxa), len(taxa)), dtype=np.float32)
    for i, left_name in enumerate(taxa):
        for j in range(i + 1, len(taxa)):
            distance = float(tree.distance(left_name, taxa[j]))
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    id2index = {taxon_id: index for index, taxon_id in enumerate(taxa)}
    return distance_matrix, taxa, id2index


def _ete3_row_chunk(start_index, end_index, taxa, node_map):
    chunk_results = []
    for i in range(start_index, end_index):
        left_node = node_map[taxa[i]]
        row_values = []
        for j in range(i + 1, len(taxa)):
            row_values.append((j, float(left_node.get_distance(node_map[taxa[j]]))))
        chunk_results.append((i, row_values))
    return chunk_results


def _load_tree_distance_matrix_ete3(tree_path: str, num_workers: int = 0):
    from ete3 import Tree

    tree = Tree(tree_path, format=1)
    taxa = [leaf.name for leaf in tree.get_leaves()]
    if any(name is None for name in taxa):
        raise ValueError("All leaf nodes in the backbone tree must have names.")

    node_map = {leaf.name: leaf for leaf in tree.get_leaves()}
    distance_matrix = np.zeros((len(taxa), len(taxa)), dtype=np.float32)
    worker_count = max(1, num_workers or min(os.cpu_count() or 1, 8))

    if len(taxa) <= 100 or worker_count == 1:
        row_chunks = [_ete3_row_chunk(0, len(taxa), taxa, node_map)]
    else:
        chunk_size = max(1, len(taxa) // worker_count)
        ranges = [(start, min(start + chunk_size, len(taxa))) for start in range(0, len(taxa), chunk_size)]
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            row_chunks = list(
                executor.map(lambda item: _ete3_row_chunk(item[0], item[1], taxa, node_map), ranges)
            )

    for chunk in row_chunks:
        for row_index, row_values in chunk:
            for col_index, distance in row_values:
                distance_matrix[row_index, col_index] = distance
                distance_matrix[col_index, row_index] = distance

    id2index = {taxon_id: index for index, taxon_id in enumerate(taxa)}
    return distance_matrix, taxa, id2index


def load_tree_distance_matrix(tree_path: str, backend: str = "auto", num_workers: int = 0):
    backend = backend.lower()
    if backend not in {"auto", "ete3", "biopython"}:
        raise ValueError("`backend` must be one of: auto, ete3, biopython.")

    if backend in {"auto", "ete3"}:
        try:
            return _load_tree_distance_matrix_ete3(tree_path, num_workers=num_workers)
        except ImportError:
            if backend == "ete3":
                raise

    return _load_tree_distance_matrix_biopython(tree_path)


def normalize_split_ids(
    split_file: str,
    valid_aliases: Iterable[str] = ("valid", "val", "validation"),
) -> Dict[str, List[str]]:
    dataframe = pd.read_csv(split_file)
    columns = {column.lower(): column for column in dataframe.columns}

    valid_column = None
    for alias in valid_aliases:
        if alias in columns:
            valid_column = columns[alias]
            break

    if "train" in columns and (valid_column or "test" in columns):

        return {
            "train_id": dataframe[columns["train"]].dropna().astype(str).tolist(),
            "valid_id": dataframe[valid_column].dropna().astype(str).tolist() if valid_column else [],
            "test_id": dataframe[columns["test"]].dropna().astype(str).tolist() if "test" in columns else [],
        }

    if {"id", "split"}.issubset(columns):
        split_ids = {"train_id": [], "valid_id": [], "test_id": []}
        valid_names = set(valid_aliases)
        for _, row in dataframe.iterrows():
            sample_id = str(row[columns["id"]])
            split_name = str(row[columns["split"]]).strip().lower()
            if split_name == "train":
                split_ids["train_id"].append(sample_id)
            elif split_name in valid_names:
                split_ids["valid_id"].append(sample_id)
            elif split_name == "test":
                split_ids["test_id"].append(sample_id)
            else:
                raise ValueError(f"Unsupported split label `{split_name}` in `{split_file}`.")
        return split_ids

    raise ValueError(
        "Split file must contain a `train` column plus optional `valid`/`test` columns, "
        "or an `id/split` table."
    )


def create_random_split(
    all_ids: Iterable[str],
    train_ratio: float,
    valid_ratio: float,
    test_ratio: float,
    seed: int,
) -> Dict[str, List[str]]:
    ratio_sum = train_ratio + valid_ratio + test_ratio
    if not np.isclose(ratio_sum, 1.0, atol=1e-6):
        raise ValueError(
            f"`train_ratio + valid_ratio + test_ratio` must equal 1.0, got {ratio_sum:.6f}."
        )

    ids = list(all_ids)
    rng = np.random.default_rng(seed)
    rng.shuffle(ids)

    train_end = int(len(ids) * train_ratio)
    valid_end = train_end + int(len(ids) * valid_ratio)
    return {
        "train_id": ids[:train_end],
        "valid_id": ids[train_end:valid_end],
        "test_id": ids[valid_end:],
    }


def validate_split_ids(split_ids: Dict[str, List[str]], available_ids: Iterable[str]) -> None:
    available_set = set(available_ids)
    train_set = set(split_ids["train_id"])
    valid_set = set(split_ids["valid_id"])
    test_set = set(split_ids["test_id"])

    missing = (train_set | valid_set | test_set).difference(available_set)
    if missing:
        preview = sorted(missing)[:10]
        raise ValueError(f"Split file contains IDs missing from the FASTA/tree inputs: {preview}")

    overlap = {
        "train/valid": train_set & valid_set,
        "train/test": train_set & test_set,
        "valid/test": valid_set & test_set,
    }
    duplicated = {key: value for key, value in overlap.items() if value}
    if duplicated:
        summary = {key: sorted(list(value))[:10] for key, value in duplicated.items()}
        raise ValueError(f"Split IDs overlap across splits: {summary}")


def construct_pair_dataframe(
    anchor_ids: Iterable[str],
    sequence_by_id: Dict[str, str],
    positive_ids: Iterable[str],
    distance_matrix: np.ndarray,
    id2index: Dict[str, int],
    include_self_pairs: bool = False,
) -> pd.DataFrame:
    anchor_ids = list(anchor_ids)
    positive_ids = list(positive_ids)
    same_id_order = anchor_ids == positive_ids
    rows = []

    if same_id_order:
        for anchor_index, left_id in enumerate(anchor_ids):
            end_index = anchor_index + 1 if include_self_pairs else anchor_index
            for positive_index in range(0, end_index):
                right_id = positive_ids[positive_index]
                rows.append(
                    {
                        "id1": left_id,
                        "id2": right_id,
                        "seq1": sequence_by_id[left_id],
                        "seq2": sequence_by_id[right_id],
                        "distance": float(distance_matrix[id2index[left_id], id2index[right_id]]),
                    }
                )
    else:
        for left_id in anchor_ids:
            for right_id in positive_ids:
                rows.append(
                    {
                        "id1": left_id,
                        "id2": right_id,
                        "seq1": sequence_by_id[left_id],
                        "seq2": sequence_by_id[right_id],
                        "distance": float(distance_matrix[id2index[left_id], id2index[right_id]]),
                    }
                )

    return pd.DataFrame(rows, columns=["id1", "id2", "seq1", "seq2", "distance"])


def save_stage4_dataset_artifacts(
    output_dir: str,
    split_ids: Dict[str, List[str]],
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    test_df: pd.DataFrame,
    metadata: Dict[str, object],
) -> None:
    dataset_dir = os.path.join(output_dir, "dataset")
    os.makedirs(dataset_dir, exist_ok=True)

    for split_name, ids in split_ids.items():
        file_name = split_name.replace("_id", "_ids.txt")
        with open(os.path.join(dataset_dir, file_name), "w", encoding="utf-8") as handle:
            handle.write("\n".join(ids))

    train_df.to_csv(os.path.join(dataset_dir, "train_pairs.csv"), index=False)
    valid_df.to_csv(os.path.join(dataset_dir, "valid_pairs.csv"), index=False)
    test_df.to_csv(os.path.join(dataset_dir, "test_pairs.csv"), index=False)

    with open(os.path.join(dataset_dir, "split_metadata.json"), "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def matrix_to_lower_triangle(distance_matrix: np.ndarray) -> List[List[float]]:
    matrix = np.asarray(distance_matrix)
    return [[float(matrix[i, j]) for j in range(i + 1)] for i in range(matrix.shape[0])]


def write_distance_csv(distance_matrix: np.ndarray, labels: List[str], output_path: str) -> None:
    dataframe = pd.DataFrame(distance_matrix, index=labels, columns=labels)
    dataframe.to_csv(output_path)


class Stage4DistanceDataset(Dataset):
    def __init__(self, dataframe, tokenizer, max_length: int = 2048):
        self.data = dataframe.reset_index(drop=True)
        self.tokenizer = tokenizer
        self.max_length = max_length

        required_columns = {"seq1", "seq2", "distance"}
        missing = required_columns.difference(self.data.columns)
        if missing:
            raise ValueError(f"Missing required stage-4 columns: {sorted(missing)}")

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        row = self.data.iloc[index]
        return {
            "seq1": str(row["seq1"]),
            "seq2": str(row["seq2"]),
            "labels": np.float32(row["distance"]),
        }


class Stage4DataCollator:
    def __init__(self, tokenizer, max_length: int = 2048):
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __call__(self, features):
        seq1 = [feature["seq1"] for feature in features]
        seq2 = [feature["seq2"] for feature in features]
        labels = torch.tensor([feature["labels"] for feature in features], dtype=torch.float32)
        batch_size = len(features)

        # Match the legacy stage-4 trainer more closely by padding only to the
        # longest sequence in the current batch instead of always padding to max_length.
        # Tokenize the two sequence columns together so both branches share the
        # same padded length and can be stacked into [batch, 2, seq_len].
        encoding = self.tokenizer(
            seq1 + seq2,
            truncation=True,
            padding=True,
            max_length=self.max_length,
            return_tensors="pt",
        )

        return {
            "input_ids": torch.stack(
                (encoding["input_ids"][:batch_size], encoding["input_ids"][batch_size:]),
                dim=1,
            ),
            "attention_mask": torch.stack(
                (encoding["attention_mask"][:batch_size], encoding["attention_mask"][batch_size:]),
                dim=1,
            ),
            "labels": labels,
        }


def pretokenize_sequences(sequences: Dict[str, str], tokenizer, max_length: int) -> Dict[str, List[int]]:
    """Tokenize every unique sequence exactly once using per-record calls,
    because batched tokenizer input can drop per-record state."""
    tokenized: Dict[str, List[int]] = {}
    for record_id, sequence in sequences.items():
        encoded = tokenizer(sequence, truncation=True, max_length=max_length)["input_ids"]
        ids = encoded[0] if (isinstance(encoded, (list, tuple)) and len(encoded) == 1 and isinstance(encoded[0], (list, tuple))) else encoded
        tokenized[record_id] = [int(token_id) for token_id in ids]
    return tokenized


class Stage4IndexedDataset(Dataset):
    """Pair dataset over pre-tokenized sequences; pairs reference IDs, not strings."""

    def __init__(self, pair_dataframe: pd.DataFrame, tokenized_sequences: Dict[str, List[int]]):
        required_columns = {"id1", "id2", "distance"}
        missing = required_columns.difference(pair_dataframe.columns)
        if missing:
            raise ValueError(
                f"Stage4IndexedDataset expects columns {sorted(required_columns)}, missing: {sorted(missing)}"
            )
        missing_ids = set(pair_dataframe["id1"]).union(pair_dataframe["id2"]).difference(tokenized_sequences)
        if missing_ids:
            raise ValueError(
                f"Pairs reference IDs missing from tokenized sequences: {sorted(missing_ids)[:10]}"
            )
        self.id1 = pair_dataframe["id1"].astype(str).to_numpy()
        self.id2 = pair_dataframe["id2"].astype(str).to_numpy()
        self.distance = pair_dataframe["distance"].to_numpy(dtype=np.float32)
        self.tokenized = tokenized_sequences
        self.lengths = np.maximum(
            [len(tokenized_sequences[record_id]) for record_id in self.id1],
            [len(tokenized_sequences[record_id]) for record_id in self.id2],
        ).astype(np.int64)

    def __len__(self):
        return len(self.id1)

    def __getitem__(self, index):
        left_id = self.id1[index]
        right_id = self.id2[index]
        return {
            "ids1": self.tokenized[left_id],
            "ids2": self.tokenized[right_id],
            "labels": self.distance[index],
            "length": int(self.lengths[index]),
            "anchor": left_id,
        }


class Stage4IndexedCollator:
    """Pads pre-tokenized pairs into [batch, 2, seq_len] batches.

    Both branches share one padded length (the batch-wide max) so the model's
    merged single-pass forward stays available. Pad positions use the tokenizer
    pad id and a zero attention mask, matching the legacy collator byte-for-byte
    on untruncated inputs (covered by tests).
    """

    def __init__(self, pad_token_id: int):
        self.pad_token_id = int(pad_token_id)

    def __call__(self, features):
        seqs1 = [list(feature["ids1"]) for feature in features]
        seqs2 = [list(feature["ids2"]) for feature in features]
        labels = torch.tensor([float(feature["labels"]) for feature in features], dtype=torch.float32)
        max_len = max(max(len(s) for s in seqs1), max(len(s) for s in seqs2))

        def pad_branch(seqs):
            ids = torch.full((len(seqs), max_len), self.pad_token_id, dtype=torch.long)
            mask = torch.zeros((len(seqs), max_len), dtype=torch.long)
            for row, seq in enumerate(seqs):
                ids[row, : len(seq)] = torch.tensor(seq, dtype=torch.long)
                mask[row, : len(seq)] = 1
            return ids, mask

        ids1, mask1 = pad_branch(seqs1)
        ids2, mask2 = pad_branch(seqs2)
        anchor_groups, group_map = [], {}
        for feature in features:
            anchor = str(feature.get("anchor", ""))
            if anchor not in group_map:
                group_map[anchor] = len(group_map)
            anchor_groups.append(group_map[anchor])
        return {
            "input_ids": torch.stack((ids1, ids2), dim=1),
            "attention_mask": torch.stack((mask1, mask2), dim=1),
            "labels": labels,
            "anchor_groups": torch.tensor(anchor_groups, dtype=torch.long),
        }


def subsample_pairs_by_distance(pair_dataframe: pd.DataFrame, budget: int, seed: int = 0) -> pd.DataFrame:
    """Distance-stratified subsampling that keeps the label distribution.

    Used to cap per-epoch validation cost: full valid pairs cost ~0.13*n^2
    which rivals the training epoch itself. Draws exactly `budget` rows across
    distance-decile strata, proportionally to stratum size.
    """
    if budget <= 0 or len(pair_dataframe) <= budget:
        return pair_dataframe.reset_index(drop=True)

    rng = np.random.default_rng(seed)
    distances = pair_dataframe["distance"].to_numpy(dtype=np.float64)
    edges = np.quantile(distances, np.linspace(0.0, 1.0, 11)[1:-1])
    strata = np.digitize(distances, edges)

    rows_by_stratum = {}
    for stratum in np.unique(strata):
        rows_by_stratum[stratum] = np.where(strata == stratum)[0]

    shares = {
        stratum: budget * len(rows) / len(pair_dataframe)
        for stratum, rows in rows_by_stratum.items()
    }
    take = {stratum: int(np.floor(share)) for stratum, share in shares.items()}
    remainder = budget - sum(take.values())
    for stratum in sorted(shares, key=lambda s: shares[s] - take[s], reverse=True):
        if remainder <= 0:
            break
        if take[stratum] < len(rows_by_stratum[stratum]):
            take[stratum] += 1
            remainder -= 1

    keep = []
    for stratum, count in take.items():
        if count > 0:
            keep.extend(rng.choice(rows_by_stratum[stratum], size=count, replace=False).tolist())
    return pair_dataframe.iloc[np.sort(keep)].reset_index(drop=True)


def create_clade_split(
    tree_path: str,
    train_ratio: float,
    valid_ratio: float,
    test_ratio: float,
    seed: int,
) -> Dict[str, List[str]]:
    """Split taxa by whole clades instead of individual leaves.

    Leaves that share a recent ancestor stay in the same split, so evaluation
    measures generalization to unseen lineages rather than to near-duplicate
    neighbours of training sequences. Greedy bin-packing: repeatedly place the
    largest clade that fits a remaining quota; descend (replace a node by its
    children) when nothing fits.
    """
    from ete3 import Tree

    tree = Tree(tree_path, format=1)
    all_leaf_names = [leaf.name for leaf in tree.get_leaves()]
    total = len(all_leaf_names)
    if total == 0:
        raise ValueError("Backbone tree has no leaves.")

    rng = np.random.default_rng(seed)
    quotas = {
        "train_id": int(round(total * train_ratio)),
        "valid_id": int(round(total * valid_ratio)),
    }
    quotas["test_id"] = max(0, total - quotas["train_id"] - quotas["valid_id"])

    size_cache = {}

    def clade_size(node) -> int:
        if node.is_leaf():
            return 1
        key = id(node)
        if key not in size_cache:
            size_cache[key] = sum(clade_size(child) for child in node.children)
        return size_cache[key]

    def leaf_names(node) -> List[str]:
        return [leaf.name for leaf in node.get_leaves()]

    split_ids: Dict[str, List[str]] = {"train_id": [], "valid_id": [], "test_id": []}
    pool = list(tree.children) or [tree]

    def largest_open_quota():
        return max(quotas, key=lambda name: quotas[name])

    while pool:
        if max(quotas.values()) <= 0:
            break
        target = largest_open_quota()
        candidates = [node for node in pool if clade_size(node) <= quotas[target]]
        if candidates:
            node = max(candidates, key=clade_size)
            split_ids[target].extend(leaf_names(node))
            quotas[target] -= clade_size(node)
            pool.remove(node)
        else:
            # Nothing fits anywhere: descend the largest pooled node one level.
            biggest = max(pool, key=clade_size)
            pool.remove(biggest)
            if biggest.is_leaf():
                # A single leaf never fits only when every quota is exhausted.
                split_ids["train_id"].append(biggest.name)
            else:
                pool.extend(biggest.children)

    # Sweep leftovers into the largest remaining quota (shuffled for determinism-per-seed).
    assigned = set(split_ids["train_id"]) | set(split_ids["valid_id"]) | set(split_ids["test_id"])
    leftovers = [name for name in all_leaf_names if name not in assigned]
    if leftovers:
        rng.shuffle(leftovers)
        fallback = largest_open_quota() if max(quotas.values()) > 0 else "train_id"
        split_ids[fallback].extend(leftovers)

    for key in split_ids:
        split_ids[key] = list(dict.fromkeys(split_ids[key]))
    return split_ids


class Stage4AnchorDataset(Dataset):
    """Sequence-level dataset for anchor (within-batch all-pairs) training.

    HDEPP-style pairing: a batch of sequences is encoded once each, then every
    within-batch pair contributes a distance-regression term. Logical length
    is `repeats * n_taxa` so one epoch averages `repeats` appearances per
    sequence; a plain RandomSampler over the tiled indices provides the
    epoch-to-epoch reshuffling.
    """

    def __init__(self, taxon_ids, tokenized_sequences: Dict[str, List[int]], repeats: int = 8):
        if repeats < 1:
            raise ValueError("anchor repeats must be >= 1")
        missing = set(taxon_ids).difference(tokenized_sequences)
        if missing:
            raise ValueError(f"Taxa missing from tokenized sequences: {sorted(missing)[:10]}")
        self.taxon_ids = list(taxon_ids)
        self.tokenized = tokenized_sequences
        self.repeats = int(repeats)

    def __len__(self):
        return len(self.taxon_ids) * self.repeats

    def __getitem__(self, index):
        taxon_id = self.taxon_ids[index % len(self.taxon_ids)]
        return {"ids": self.tokenized[taxon_id], "taxon_id": taxon_id}


class AnchorPairCollator:
    """Pads an anchor batch and attaches the [B, B] ground-truth distance block."""

    def __init__(self, distance_matrix: np.ndarray, id2index: Dict[str, int], pad_token_id: int):
        self.distance_matrix = distance_matrix
        self.id2index = id2index
        self.pad_token_id = int(pad_token_id)

    def __call__(self, features):
        sequences = [list(feature["ids"]) for feature in features]
        taxon_ids = [str(feature["taxon_id"]) for feature in features]
        max_len = max(len(seq) for seq in sequences)
        input_ids = torch.full((len(sequences), max_len), self.pad_token_id, dtype=torch.long)
        attention_mask = torch.zeros((len(sequences), max_len), dtype=torch.long)
        for row, seq in enumerate(sequences):
            input_ids[row, : len(seq)] = torch.tensor(seq, dtype=torch.long)
            attention_mask[row, : len(seq)] = 1

        indices = [self.id2index[taxon_id] for taxon_id in taxon_ids]
        block = self.distance_matrix[np.ix_(indices, indices)]
        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "labels": torch.as_tensor(block, dtype=torch.float32),
        }


class DualModeCollator:
    """Dispatches between anchor-batch and pair-batch features in one collator.

    Trainer installs a single collator for both loops; anchor training steps
    send {"ids", "taxon_id"} features while (pair-level) validation steps send
    {"ids1", "ids2", "labels"} features, so evaluation keeps the exact pair
    protocol used by the pair mode for cross-mode comparability.
    """

    def __init__(self, anchor_collator: AnchorPairCollator, pair_collator: Stage4IndexedCollator):
        self.anchor_collator = anchor_collator
        self.pair_collator = pair_collator

    def __call__(self, features):
        if features and "taxon_id" in features[0]:
            return self.anchor_collator(features)
        return self.pair_collator(features)
