from typing import Optional, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers.models.esm.modeling_esm import EsmModel, EsmPreTrainedModel

from .utils import Pooler


def minkowski_inner_product(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    x = x.to(torch.float32)
    y = y.to(torch.float32)
    return -(x[..., 0] * y[..., 0]) + (x[..., 1:] * y[..., 1:]).sum(dim=-1)


def hyperboloid_distance(
    x: torch.Tensor,
    y: torch.Tensor,
    curvature: Optional[Union[torch.Tensor, float]] = None,
    epsilon: float = 1e-7,
) -> torch.Tensor:
    inner = torch.clamp(-minkowski_inner_product(x, y), min=1.0 + epsilon)
    distance = torch.acosh(inner)

    if curvature is None:
        return distance

    curvature_tensor = torch.as_tensor(curvature, device=x.device, dtype=torch.float32)
    scale_factor = 1.0 / torch.sqrt(torch.abs(curvature_tensor).clamp_min(epsilon))
    return distance * scale_factor


class HyperboloidEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.projection = nn.Linear(dim, dim + 1)

    def project(self, x: torch.Tensor) -> torch.Tensor:
        x_proj = self.projection(x.to(torch.float32))
        spatial = x_proj[..., 1:]
        time_component = torch.sqrt(1.0 + spatial.pow(2).sum(dim=-1, keepdim=True) + 1e-6)
        return torch.cat([time_component, spatial], dim=-1)


class PairEmbeddingDistance(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.curvature = nn.Parameter(torch.tensor([-1.0], dtype=torch.float32))
        self.embedding = HyperboloidEmbedding(dim)
        nn.init.uniform_(self.curvature, -2.0, -0.5)

    def _apply_constraints(self) -> None:
        with torch.no_grad():
            self.curvature.clamp_(max=-1e-6)

    def forward(self, emb1: torch.Tensor, emb2: torch.Tensor) -> torch.Tensor:
        self._apply_constraints()
        emb1_h = self.embedding.project(emb1)
        emb2_h = self.embedding.project(emb2)
        return hyperboloid_distance(emb1_h, emb2_h, curvature=self.curvature)


class LucaPhyloHyper(EsmPreTrainedModel):
    _keys_to_ignore_on_load_missing = [r"position_ids"]

    def __init__(self, config, *model_args, **model_kwargs):
        super().__init__(config)
        self.model_args = model_kwargs.get("model_args")
        pooler_type = getattr(self.model_args, "pooler_type", "avg")
        self.loss_weighted_method = getattr(self.model_args, "loss_weighted_method", "ols")
        self.loss_weight_eps = getattr(self.model_args, "loss_weight_eps", 1e-2)
        self.loss_target_transform = getattr(self.model_args, "loss_target_transform", "none")
        if self.loss_target_transform not in ("none", "log", "sqrt"):
            raise ValueError(
                f"Unsupported loss_target_transform `{self.loss_target_transform}`; "
                "expected one of: none, log, sqrt."
            )
        if self.loss_weighted_method not in ("ols", "fm", "be"):
            raise ValueError(
                f"Unsupported loss_weighted_method `{self.loss_weighted_method}`; "
                "expected one of: ols, fm, be."
            )
        self.rank_lambda = float(getattr(self.model_args, "rank_lambda", 0.0))
        self.rank_margin = float(getattr(self.model_args, "rank_margin", 0.05))

        self.esm = EsmModel(config, add_pooling_layer=False)
        self.pooler = Pooler(pooler_type)
        self.scale = PairEmbeddingDistance(config.hidden_size)
        self.init_weights()

    def _encode(self, input_ids: torch.Tensor, attention_mask: torch.Tensor):
        return self.esm(
            input_ids=input_ids,
            attention_mask=attention_mask,
            output_hidden_states=self.pooler.requires_hidden_states,
            return_dict=True,
        )

    def encode_sequences(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        outputs = self._encode(input_ids=input_ids, attention_mask=attention_mask)
        return self.pooler(attention_mask, outputs).to(torch.float32)

    def predict_distance_from_embeddings(
        self,
        anchor_embeddings: torch.Tensor,
        positive_embeddings: torch.Tensor,
    ) -> torch.Tensor:
        return self.scale(anchor_embeddings.to(torch.float32), positive_embeddings.to(torch.float32))

    def _transform_labels(self, labels: torch.Tensor) -> torch.Tensor:
        if self.loss_target_transform == "log":
            return torch.log(labels.clamp_min(self.loss_weight_eps))
        if self.loss_target_transform == "sqrt":
            return torch.sqrt(labels)
        return labels

    def _inverse_transform(self, raw: torch.Tensor) -> torch.Tensor:
        """Map model outputs back to the original distance space so that
        logits/metrics/inference stay comparable across transforms."""
        if self.loss_target_transform == "log":
            return torch.exp(raw.clamp_max(15.0))
        if self.loss_target_transform == "sqrt":
            return raw ** 2
        return raw

    def forward(self, input_ids=None, attention_mask=None, labels=None, anchor_groups=None, **kwargs):
        # Anchor (HDEPP-style) mode: 2-D inputs encode a single set of
        # sequences; all within-batch pairs regress at once against a [B, B]
        # distance block, and each sequence is encoded exactly once per step.
        if input_ids is not None and input_ids.dim() == 2:
            embeddings = self.encode_sequences(input_ids, attention_mask)
            raw = self.predict_distance_from_embeddings(
                embeddings.unsqueeze(1), embeddings.unsqueeze(0)
            )
            loss = None
            if labels is not None:
                labels = labels.to(raw.device, dtype=raw.dtype)
                t_labels = self._transform_labels(labels)
                off_diagonal = ~torch.eye(raw.size(0), dtype=torch.bool, device=raw.device)
                loss = ((raw[off_diagonal] - t_labels[off_diagonal]) ** 2).mean()
                if self.rank_lambda > 0 and self.training:
                    # Pairwise hinge over every within-batch ordering: for each
                    # anchor i, force raw[i, nearer] + margin <= raw[i, farther].
                    # Optimizes the quantity placement is actually scored on
                    # (nearest-neighbour ordering), which MSE alone underfits.
                    # Gated on self.training so eval_loss stays pure MSE.
                    diff_labels = t_labels.unsqueeze(2) - t_labels.unsqueeze(1)  # [B, nearer, farther]
                    order_mask = diff_labels < 0
                    hinge = F.relu(
                        raw.unsqueeze(2) - raw.unsqueeze(1) + self.rank_margin
                    )
                    rank_loss = hinge[order_mask].mean()
                    loss = loss + self.rank_lambda * rank_loss
            return SequenceClassifierOutput(
                loss=loss,
                logits=self._inverse_transform(raw),
                hidden_states=None,
                attentions=None,
            )

        if isinstance(input_ids, (list, tuple)):
            input_ids_1, input_ids_2 = input_ids
            attention_mask_1, attention_mask_2 = attention_mask
        else:
            if input_ids is None or attention_mask is None or input_ids.dim() != 3:
                raise ValueError(
                    "Stage-4 expects either a list of two tensors or tensors with shape "
                    "[batch, 2, sequence_length]."
                )
            input_ids_1 = input_ids[:, 0]
            input_ids_2 = input_ids[:, 1]
            attention_mask_1 = attention_mask[:, 0]
            attention_mask_2 = attention_mask[:, 1]

        # Merge the two branches into a single ESM forward pass. ESM layers are
        # per-token / per-sample (LayerNorm + self-attention over the batch dim
        # independently), so concatenating along the batch dimension is
        # numerically equivalent to two separate forward passes -- only the
        # launch/sync overhead is saved. Peak activation memory is unchanged,
        # since both branches' activations are needed for backward anyway.
        # The standard collator stacks the two branches with a shared padded
        # length, so they share seq_len and can be concatenated. When the two
        # branches come in with different lengths (e.g. a list/tuple call with
        # independently padded tensors), fall back to two passes.
        batch_size = input_ids_1.size(0)
        if input_ids_1.size(1) == input_ids_2.size(1):
            merged_ids = torch.cat((input_ids_1, input_ids_2), dim=0)
            merged_mask = torch.cat((attention_mask_1, attention_mask_2), dim=0)
            merged = self.encode_sequences(merged_ids, merged_mask)
            z1, z2 = merged[:batch_size], merged[batch_size:]
        else:
            z1 = self.encode_sequences(input_ids_1, attention_mask_1)
            z2 = self.encode_sequences(input_ids_2, attention_mask_2)
        raw = self.predict_distance_from_embeddings(z1, z2)

        loss = None
        if labels is not None:
            labels = labels.to(raw.device, dtype=raw.dtype).view_as(raw)
            t_labels = self._transform_labels(labels)
            if self.loss_weighted_method == "ols":
                loss = F.mse_loss(raw, t_labels)
            else:
                # HDEPP-style distance weighting (Biology 2022): scale each
                # squared error with a power of the true distance so close
                # pairs dominate (fm = Fitch-Margoliash 1/D^2, be = Beyer
                # 1/D). Weights are mean-normalised to keep the loss on the
                # same scale as plain MSE, and labels are clamped so
                # near-zero patristic distances cannot blow the weight up
                # under fp16. Weighting uses the ORIGINAL distance scale.
                ref = labels.clamp_min(self.loss_weight_eps)
                if self.loss_weighted_method == "fm":
                    weight = 1.0 / ref.pow(2)
                else:
                    weight = 1.0 / ref
                weight = weight / weight.mean()
                loss = ((raw - t_labels) ** 2 * weight).mean()
            if anchor_groups is not None and self.rank_lambda > 0 and self.training:
                # Group-wise ranking hinge: within each anchor group, force
                # raw[nearer] + margin <= raw[farther]. Eval stays pure MSE so
                # early stopping remains comparable across rank on/off arms.
                same_group = anchor_groups.unsqueeze(1) == anchor_groups.unsqueeze(0)
                closer = t_labels.unsqueeze(1) < t_labels.unsqueeze(0)
                not_self = ~torch.eye(raw.size(0), dtype=torch.bool, device=raw.device)
                mask = same_group & closer & not_self
                if mask.any():
                    hinge = F.relu(
                        raw.unsqueeze(1) - raw.unsqueeze(0) + self.rank_margin
                    )
                    loss = loss + self.rank_lambda * hinge[mask].mean()

        return SequenceClassifierOutput(
            loss=loss,
            logits=self._inverse_transform(raw),
            hidden_states=None,
            attentions=None,
        )


def read_target_transform_from_run_config(adapter_dir: str) -> str:
    """Read loss_target_transform from a trained run's stage4_run_config.json.

    Defaults to "none" when the file/key is absent (legacy runs predate the flag).
    """
    import json
    import os

    path = os.path.join(adapter_dir, "stage4_run_config.json")
    if not os.path.exists(path):
        return "none"
    with open(path, "r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    return cfg.get("model_args", {}).get("loss_target_transform", "none")
