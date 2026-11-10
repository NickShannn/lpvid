import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import os
from typing import List, Union, Dict, Optional, Tuple
from .src_models import Pooler


class SeqDataset(Dataset):
    """序列数据集类"""
    def __init__(self, sequences: List[str], tokenizer, max_length: int = 1024):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.max_length = max_length
        self.max_seq_len = min(
            max(len(seq) for seq in sequences),
            self.max_length
        )

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx]
        seq = seq[-self.max_length:] # 截断序列,保留最后右侧max_length个氨基酸
        
        encoding = self.tokenizer(
            seq,
            truncation=True,
            max_length=self.max_seq_len + 2,
            padding='max_length',
            return_tensors='pt'
        )
        return {
            'input_ids': encoding['input_ids'].squeeze(),
            'attention_mask': encoding['attention_mask'].squeeze()
        }

class SequenceEmbedder:
    """
    序列嵌入计算器
    用于获取蛋白质序列的嵌入向量表示
    
    Attributes:
        model_path (str): 预训练模型的路径
        device (str): 计算设备
        batch_size (int): 批处理大小
        max_length (int): 最大序列长度
        pooler_type (str): 池化方式 ('avg', 'cls', 'cls_before_pooler', 'avg_top2', 'avg_first_last')
    """
    
    def __init__(self, 
                 model_path: str, 
                 device: str = 'cuda:0',
                 batch_size: int = 32,
                 max_length: int = 1-24,
                 pooler_type: str = 'avg'):
        """
        初始化序列嵌入计算器
        
        Args:
            model_path: 预训练模型的路径
            device: 计算设备 ('cuda:0', 'cuda:1', ..., 'cpu')
            batch_size: 批处理大小
            max_length: 最大序列长度
            pooler_type: 池化方式
        """
        self.model_path = model_path
        self.device = device
        self.batch_size = batch_size
        self.max_length = max_length
        self.pooler_type = pooler_type
        
        # 初始化分词器和模型
        self.tokenizer = AutoTokenizer.from_pretrained(model_path)
        self.model = AutoModel.from_pretrained(model_path)
        # self.tokenizer = AutoTokenizer.from_pretrained(model_path,ignore_mismatched_sizes=True)
        # self.model = AutoModel.from_pretrained(model_path,ignore_mismatched_sizes=True)

        self.model.to(device)
        
        # 初始化池化层
        self.pooler = Pooler(pooler_type)
        
    def get_embeddings(self, 
                      sequences: Union[List[str], str],
                      return_numpy: bool = False) -> Union[torch.Tensor, np.ndarray]:
        """
        获取序列的嵌入向量
        
        Args:
            sequences: 单个序列或序列列表
            return_numpy: 是否返回numpy数组
            
        Returns:
            torch.Tensor or np.ndarray: 序列的嵌入向量
        """
        # 处理单个序列输入
        if isinstance(sequences, str):
            sequences = [sequences]
            
        # 创建数据加载器
        dataset = SeqDataset(sequences, self.tokenizer, self.max_length)
        dataloader = DataLoader(dataset, 
                              batch_size=self.batch_size, 
                              shuffle=False)
        
        embeddings = []
        self.model.eval()
        
        with torch.no_grad():
            for batch in tqdm(dataloader, desc="计算序列嵌入向量"):
                input_ids = batch['input_ids'].to(self.device)
                attention_mask = batch['attention_mask'].to(self.device)
                
                outputs = self.model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                    output_hidden_states=True,
                    return_dict=True
                )
                
                # 使用Pooler获取序列表示
                batch_embeddings = self.pooler(attention_mask, outputs)
                embeddings.append(batch_embeddings.cpu())
        
        embeddings = torch.cat(embeddings, dim=0)
        
        return embeddings.numpy() if return_numpy else embeddings
    
    def get_embeddings_from_fasta(self, 
                                fasta_file: str,
                                return_dict: bool = False
                                ) -> Union[Dict[str, torch.Tensor], 
                                         Dict[str, np.ndarray],
                                         torch.Tensor,
                                         np.ndarray]:
        """
        从FASTA文件获取序列的嵌入向量
        
        Args:
            fasta_file: FASTA文件路径
            return_dict: 是否返回字典形式 {序列ID: 嵌入向量}
            
        Returns:
            embeddings: 嵌入向量或字典
        """
        sequences = []
        seq_ids = []
        
        # 读取FASTA文件
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_ids.append(record.id)
            sequences.append(str(record.seq))
            
        # 获取嵌入向量
        embeddings = self.get_embeddings(sequences)
        
        if return_dict:
            return {seq_id: emb for seq_id, emb in zip(seq_ids, embeddings)}
        return embeddings
    
    def save_embeddings(self, 
                       embeddings: Union[torch.Tensor, np.ndarray],
                       output_file: str,
                       metadata: Optional[pd.DataFrame] = None):
        """
        保存嵌入向量到文件
        
        Args:
            embeddings: 嵌入向量 (torch.Tensor 或 numpy.ndarray)
            output_file: 输出文件路径 (.pt文件)
            metadata: 可选的元数据DataFrame，包含seq_ids等信息
        """
        # 确保embeddings是torch.Tensor类型
        if isinstance(embeddings, np.ndarray):
            embeddings = torch.from_numpy(embeddings)
            
        # 创建输出目录
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # 保存嵌入向量
        torch.save(embeddings, output_file)
        
        # 如果提供了元数据，保存到同名的csv文件
        if metadata is not None:
            metadata_file = output_file.rsplit('.', 1)[0] + '_metadata.csv'
            metadata.to_csv(metadata_file, index=False)
    
    @staticmethod
    def load_embeddings(input_file: str) -> Union[torch.Tensor, Tuple[torch.Tensor, pd.DataFrame]]:
        """
        从文件加载嵌入向量
        
        Args:
            input_file: 输入文件路径 (.pt文件)
            
        Returns:
            如果存在元数据文件，返回(embeddings, metadata)元组
            否则只返回embeddings
        """
        # 加载嵌入向量
        embeddings = torch.load(input_file)
        
        # 检查是否存在对应的元数据文件
        metadata_file = input_file.rsplit('.', 1)[0] + '_metadata.csv'
        if os.path.exists(metadata_file):
            metadata = pd.read_csv(metadata_file)
            return embeddings, metadata
        
        return embeddings

def main():
    """使用示例"""
    # 初始化嵌入计算器
    model_path = "path/to/your/model"
    embedder = SequenceEmbedder(
        model_path=model_path,
        device='cuda:0',
        batch_size=32,
        # max_length=2048,
        max_length=1024,
        pooler_type='avg'
    )
    
    # 示例1：从序列列表获取嵌入向量
    sequences = [
        "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG",
        "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    ]
    embeddings = embedder.get_embeddings(sequences)
    print("Embeddings shape:", embeddings.shape)
    
    # 示例2：从FASTA文件获取嵌入向量
    fasta_file = "path/to/your/sequences.fasta"
    embeddings_dict = embedder.get_embeddings_from_fasta(fasta_file, return_dict=True)
    
    # 示例3：保存和加载嵌入向量
    output_file = "embeddings.pt"
    metadata = pd.DataFrame({
        'seq_id': ['seq1', 'seq2'],
        'description': ['sequence 1', 'sequence 2']
    })
    
    # 保存嵌入向量和元数据
    embedder.save_embeddings(embeddings, output_file, metadata)
    
    # 加载嵌入向量和元数据
    loaded_embeddings, loaded_metadata = SequenceEmbedder.load_embeddings(output_file)
    print("Loaded embeddings shape:", loaded_embeddings.shape)
    print("Loaded metadata:\n", loaded_metadata)

if __name__ == "__main__":
    main()