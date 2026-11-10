"""跨模型统一的原始预测记录。"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelPrediction:
    query_id: str
    prob: float
    label: str = ""
    status: str = "ok"
