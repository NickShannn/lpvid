import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.distributed as dist

import transformers
from transformers.models.esm.modeling_esm import EsmPreTrainedModel, EsmModel, EsmLMHead

from transformers.activations import gelu
from transformers.file_utils import (
    add_code_sample_docstrings,
    add_start_docstrings,
    add_start_docstrings_to_model_forward,
    replace_return_docstrings,
)
from transformers.modeling_outputs import SequenceClassifierOutput, BaseModelOutputWithPoolingAndCrossAttentions



class MLPLayer(nn.Module):
    """
    Head for getting sentence representations over LucaPhylo's CLS representation.
    """

    def __init__(self, config):
        super().__init__()
        self.dense = nn.Linear(config.hidden_size, config.hidden_size)
        self.activation = nn.Tanh()

    def forward(self, features, **kwargs):
        x = self.dense(features)
        x = self.activation(x)

        return x

def square_distance(t1_emb, t2_emb):
    D = t1_emb - t2_emb
    d = torch.sum(D * D, dim=-1)
    return d
def hyperbolic_distance(u, v, epsilon=1e-7):
    # 确保输入向量范数小于1
    u_norm = torch.sum(u * u, dim=-1)
    v_norm = torch.sum(v * v, dim=-1)
    
    # 检查是否满足双曲空间约束
    if torch.any(u_norm >= 1) or torch.any(v_norm >= 1):
        raise ValueError("Input vectors must have norm < 1 for hyperbolic distance")
    
    # 计算双曲距离
    diff_squared = torch.sum((u - v) ** 2, dim=-1)
    alpha = (1 - u_norm)
    beta = (1 - v_norm)
    
    # 增加数值稳定性
    gamma = 1 + 2 * diff_squared / (alpha * beta) + epsilon
    delta = torch.sqrt(gamma * gamma - 1 + epsilon)
    
    return torch.log(gamma + delta)
def hyperbolic_distance_old(u, v, epsilon=1e-7):  # changed from epsilon=1e-7 to reduce error
    sqdist = torch.sum((u - v) ** 2, dim=-1)
    squnorm = torch.sum(u ** 2, dim=-1)
    sqvnorm = torch.sum(v ** 2, dim=-1)
    x = 1 + 2 * sqdist / ((1 - squnorm) * (1 - sqvnorm)) + epsilon
    z = torch.sqrt(x ** 2 - 1)
    return torch.log(x + z)

class Pooler(nn.Module):
    """
    Parameter-free poolers to get the sentence embedding
    'cls': [CLS] representation with ESM's MLP pooler.
    'cls_before_pooler': [CLS] representation without the original MLP pooler.
    'avg': average of the last layers' hidden states at each token.
    'avg_top2': average of the last two layers.
    'avg_first_last': average of the first and the last layers.
    """
    def __init__(self, pooler_type):
        super().__init__()
        self.pooler_type = pooler_type
        assert self.pooler_type in ["cls", "cls_before_pooler", "avg", "avg_top2", "avg_first_last"], "unrecognized pooling type %s" % self.pooler_type

    def forward(self, attention_mask, outputs):
        last_hidden = outputs.last_hidden_state
        pooler_output = outputs.pooler_output
        hidden_states = outputs.hidden_states

        if self.pooler_type in ['cls_before_pooler', 'cls']:
            return last_hidden[:, 0]
        elif self.pooler_type == "avg":




            # 下面这个方面在pooling算平均时，可能包括了[CLS] 和[EOS] token的embedding
            # return ((last_hidden * attention_mask.unsqueeze(-1)).sum(1) / attention_mask.sum(-1).unsqueeze(-1))
        

        # 复制一份 mask 以免影响其他计算
            # 注意：这里假设 batch 中所有序列都有 EOS，且 attention_mask 能够正确反映长度
            # 但由于 padding 的存在，EOS 的位置在每个序列中是不固定的，处理起来比较麻烦。
            
            # 更稳妥的方法是利用 input_ids (如果能传进来) 识别特殊 token，但这里只传了 mask。
            # 对于 ESM，通常 CLS 是第 0 位，EOS 是最后一个非 padding 位。
            
            # --- 简化且通用的修正方案 ---
            # 既然难以精确定位变长的 EOS，我们可以直接切片去掉第 0 位（CLS）
            # 对于 EOS，因为它是随长度变化的，如果不想改动太多数据加载逻辑，
            # 很多研究实际上只去掉了 CLS，或者保留现状。
            
            # 但如果我们追求严谨，可以利用 attention_mask 计算长度来定位 EOS
            
            seq_lengths = attention_mask.sum(1).long() # 获取每个序列的真实长度
            
            # 创建一个新的 mask，排除首尾
            # 1. 把所有位置为 1 的 mask 拷贝
            mask_no_special = attention_mask.clone()
            
            # 2. 去掉 CLS (第 0 列)
            mask_no_special[:, 0] = 0
            
            # 3. 去掉 EOS (利用 range 索引)
            # EOS 的索引是 seq_lengths - 1
            batch_indices = torch.arange(attention_mask.size(0), device=attention_mask.device)
            # 确保长度大于0，避免索引错误
            mask_no_special[batch_indices, seq_lengths - 1] = 0
            
            # 4. 重新计算平均值
            # 注意分母也要重新计算，防止除以 0 (虽然对于正常序列不可能)
            sum_embeddings = (last_hidden * mask_no_special.unsqueeze(-1)).sum(1)
            valid_lengths = mask_no_special.sum(-1).unsqueeze(-1)
            
            # 加上一个极小值防止除零错误 (epsilon)
            return sum_embeddings / (valid_lengths + 1e-9)



        elif self.pooler_type == "avg_first_last":

            first_hidden = hidden_states[1]
            last_hidden = hidden_states[-1]
            pooled_result = ((first_hidden + last_hidden) / 2.0 * attention_mask.unsqueeze(-1)).sum(1) / attention_mask.sum(-1).unsqueeze(-1)
            return pooled_result
        elif self.pooler_type == "avg_top2":
            second_last_hidden = hidden_states[-2]
            last_hidden = hidden_states[-1]
            pooled_result = ((last_hidden + second_last_hidden) / 2.0 * attention_mask.unsqueeze(-1)).sum(1) / attention_mask.sum(-1).unsqueeze(-1)
            return pooled_result
        else:
            raise NotImplementedError


# def cl_init(cls, config):
#     """
#     Contrastive learning class init function.
#     """
#     cls.pooler_type = cls.model_args.pooler_type
#     cls.pooler = Pooler(cls.model_args.pooler_type)
#     if cls.model_args.pooler_type == "cls":
#         cls.mlp = MLPLayer(config)
#     # cls.sim = Similarity(temp=cls.model_args.temp)
#     cls.init_weights()


def distance_forward(cls,
    encoder,
    scale,
    input_ids=None,
    attention_mask=None,
    position_ids=None,
    head_mask=None,
    inputs_embeds=None,
    labels=None,
    output_attentions=None,
    output_hidden_states=None,
    return_dict=None,
):
    return_dict = return_dict if return_dict is not None else cls.config.use_return_dict
    # input_ids = torch.stack(input_ids,dim=1)
    # 将每个子列表堆叠成一个 tensor
    tensor_list = [torch.stack(sub_list) for sub_list in input_ids]

    # 将所有 tensor 堆叠成一个 tensor
    input_ids = torch.stack(tensor_list)

    tensor_list = [torch.stack(sub_list) for sub_list in attention_mask]
    attention_mask = torch.stack(tensor_list)



    batch_size = input_ids.size(0)
    # Number of sentences in one instance
    # 2: pair instance; 3: pair instance with a hard negative
    num_sent = input_ids.size(1)

    mlm_outputs = None
    # Flatten input for encoding
    input_ids = input_ids.view((-1, input_ids.size(-1))) # (bs * num_sent, len)
    attention_mask = attention_mask.view((-1, attention_mask.size(-1))) # (bs * num_sent len)
    # if token_type_ids is not None:
    #     token_type_ids = token_type_ids.view((-1, token_type_ids.size(-1))) # (bs * num_sent, len)

    # Get raw embeddings
    outputs = encoder(
        input_ids,
        attention_mask=attention_mask,
        # token_type_ids=token_type_ids,
        position_ids=position_ids,
        head_mask=head_mask,
        inputs_embeds=inputs_embeds,
        output_attentions=output_attentions,
        output_hidden_states=False,
        return_dict=True,
    )

    # Pooling
    pooler_output = cls.pooler(attention_mask, outputs)
    pooler_output = pooler_output.view((batch_size, num_sent, pooler_output.size(-1))) # (bs, num_sent, hidden)

    # If using "cls", we add an extra MLP layer
    # (same as BERT's original implementation) over the representation.
    # if cls.pooler_type == "cls":
    #     pooler_output = cls.mlp(pooler_output)

    # Separate representation
    z1, z2 = pooler_output[:,0], pooler_output[:,1]

    # Gather all embeddings if using distributed training
    if dist.is_initialized() and cls.training:
        # Gather hard negative
        if num_sent >= 3:
            z3_list = [torch.zeros_like(z3) for _ in range(dist.get_world_size())]
            dist.all_gather(tensor_list=z3_list, tensor=z3.contiguous())
            z3_list[dist.get_rank()] = z3
            z3 = torch.cat(z3_list, 0)

        # Dummy vectors for allgather
        z1_list = [torch.zeros_like(z1) for _ in range(dist.get_world_size())]
        z2_list = [torch.zeros_like(z2) for _ in range(dist.get_world_size())]
        # Allgather
        dist.all_gather(tensor_list=z1_list, tensor=z1.contiguous())
        dist.all_gather(tensor_list=z2_list, tensor=z2.contiguous())

        # Since allgather results do not have gradients, we replace the
        # current process's corresponding embeddings with original tensors
        z1_list[dist.get_rank()] = z1
        z2_list[dist.get_rank()] = z2
        # Get full batch embeddings: (bs x N, hidden)
        z1 = torch.cat(z1_list, 0)
        z2 = torch.cat(z2_list, 0)

    # cos_sim = cls.sim(z1.unsqueeze(1), z2.unsqueeze(0))
    # add distance compute
    distances = scale(z1,z2)
    # labels = labels.to(device)
    if labels!=None:
        loss_fct = nn.MSELoss()
        loss = loss_fct(distances, labels)
    else:
        loss = None


    # if not return_dict:
    #     output = (cos_sim,) + outputs[2:]
    #     return ((loss,) + output) if loss is not None else output

    return SequenceClassifierOutput(
        loss=loss,
        logits=distances,
        hidden_states=outputs.hidden_states,
        attentions=outputs.attentions,
    )

class PairEmbeddingDistance(nn.Module):
    """ Wrapper model for a general encoder, computes pairwise distances and applies projections """
    
    def __init__(self, dim=None, distance_type='hyperbolic'):
        super(PairEmbeddingDistance, self).__init__()
        self.distance_type = distance_type
        self.scaling = nn.Parameter(torch.Tensor([1.]), requires_grad=True)
        
        if distance_type == 'hyperboloid':
            if dim is None:
                raise ValueError("dim parameter is required for hyperboloid distance")
            self.embedding = HyperboloidEmbedding(dim)
        else:
            self.radius = nn.Parameter(torch.Tensor([1e-2]), requires_grad=True)

    def normalize_embeddings(self, embeddings):
        """ Project embeddings to an hypersphere of a certain radius """
        if self.distance_type == 'hyperbolic':
            # 确保半径严格小于1
            radius = self.radius.clamp(1e-7, 0.9)
            normalized = F.normalize(embeddings, p=2, dim=1) * radius
            return normalized
        elif self.distance_type == 'square':
            return F.normalize(embeddings, p=2, dim=1) * self.radius
        return embeddings  # hyperboloid 模型不需要这种归一化

        
    def normalize_embeddings_old(self, embeddings):
        """ Project embeddings to an hypersphere of a certain radius """
        min_scale = 1e-7

        if self.distance_type == 'hyperbolic':
            max_scale = 1 - 1e-3
        else:
            max_scale = 1e10

        return F.normalize(embeddings, p=2, dim=1) * self.radius.clamp_min(min_scale).clamp_max(max_scale)

    def forward(self, emb1, emb2):
        if self.distance_type == 'hyperbolic':
            distance = hyperbolic_distance(
                self.normalize_embeddings(emb1), 
                self.normalize_embeddings(emb2)
            )
        elif self.distance_type == 'hyperboloid':
            emb1_h = self.embedding.project(emb1)
            emb2_h = self.embedding.project(emb2)
            distance = hyperboloid_distance(emb1_h, emb2_h)
        elif self.distance_type == 'square':
            distance = square_distance(
                self.normalize_embeddings(emb1), 
                self.normalize_embeddings(emb2)
            )
        else:
            raise ValueError(f"Unknown distance type: {self.distance_type}")
            
        return distance * self.scaling



class LucaPhylo(EsmPreTrainedModel):
    _keys_to_ignore_on_load_missing = [r"position_ids"]

    def __init__(self, config):
        super().__init__(config)
        # self.model_args = config["model_args"]
        self.esm = EsmModel(config, add_pooling_layer=False)
        self.scale = PairEmbeddingDistance()
        self.pooler = Pooler('avg')
        self.init_weights()
        # cl_init(self, config)

    def forward(self,
        input_ids=None,
        attention_mask=None,
        # token_type_ids=None,
        position_ids=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
        sent_emb=False,
    ):

        return distance_forward(self, self.esm, self.scale,
            input_ids=input_ids,
            attention_mask=attention_mask,
            # token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            labels=labels,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
            # mlm_input_ids=mlm_input_ids,
            # mlm_labels=mlm_labels,
        )
    




# 'Loid （Hyperboloid Model） 双曲面模型


def minkowski_inner_product(x, y):
    """
    闵可夫斯基内积: -x₀y₀ + x₁y₁ + ... + xₙyₙ
    """
    # 分离时间like分量和空间like分量
    t = x[..., 0] * y[..., 0]
    spatial = torch.sum(x[..., 1:] * y[..., 1:], dim=-1)
    return -t + spatial

def hyperboloid_distance(x, y, epsilon=1e-7):
    """
    双曲面模型（Hyperboloid/Minkowski模型）的距离计算
    在 n+1 维闵可夫斯基空间中，满足 x₀² - x₁² - ... - xₙ² = 1 的上片双曲面
    """
    # 检查点是否在双曲面上 (x₀² - x₁² - ... - xₙ² = 1)
    x_norm = minkowski_inner_product(x, x)
    y_norm = minkowski_inner_product(y, y)
    
    if not (torch.allclose(x_norm, torch.ones_like(x_norm), atol=1e-4) and 
            torch.allclose(y_norm, torch.ones_like(y_norm), atol=1e-4)):
        raise ValueError("Points must lie on the hyperboloid")
    
    # 计算双曲距离
    inner = -minkowski_inner_product(x, y)
    # 确保数值稳定性
    inner = torch.clamp(inner, min=1.0 + epsilon)
    
    return torch.acosh(inner)

class HyperboloidEmbedding(nn.Module):
    """
    双曲面模型的嵌入层
    """
    def __init__(self, dim):
        super().__init__()
        self.dim = dim  # 空间维度
        
    def expmap0(self, v):
        """
        从切空间映射到双曲面（以原点为参考点）
        """
        v_norm = torch.norm(v, p=2, dim=-1, keepdim=True)
        # 避免除零
        v_norm = torch.clamp(v_norm, min=1e-7)
        
        # 计算第一个分量（时间like分量）
        t = torch.cosh(v_norm)
        # 计算其余分量（空间like分量）
        spatial = torch.sinh(v_norm) * (v / v_norm)
        
        # 组合成双曲面上的点
        return torch.cat([t, spatial], dim=-1)

    def project(self, x):
        """
        将点投影到双曲面上
        """
        # 确保时间like分量大于1
        t = torch.sqrt(1.0 + torch.sum(x[..., 1:] ** 2, dim=-1, keepdim=True))
        return torch.cat([t, x[..., 1:]], dim=-1)