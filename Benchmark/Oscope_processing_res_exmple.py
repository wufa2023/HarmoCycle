
import os
os.chdir('/root/Cycle/lattest_20250512/')
import scanpy as sc

from HarmoCycle_v1 import * 

# 1. 并行数据加载和预处理函数


import torch
import numpy as np
import random
import os

def set_seed(seed=3407):
    """固定所有随机种子以确保可复现性"""
    # 1. PyTorch
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # 多GPU情况
        torch.backends.cudnn.deterministic = True  # 确定性卷积算法
        torch.backends.cudnn.benchmark = False     # 关闭自动优化
    
    # 2. Numpy
    np.random.seed(seed)
    
    # 3. Python内置随机
    random.seed(seed)
    
    # 4. 环境变量（某些库会依赖这个）
    os.environ['PYTHONHASHSEED'] = str(seed)

# 使用示例
set_seed(222)  # 设置随机种子为42


import scanpy as sc
ref_dataset1 = sc.read_h5ad('/remote-home/share/data2/wbl/wbl_cellcycle/ground-truth/ref_dataset_Buettner.h5ad')  
ref_dataset1 = ref_dataset1[np.random.permutation(ref_dataset1.n_obs), :]

sc.pp.log1p(ref_dataset1)

import pandas as pd
gene_oscope1 = pd.read_table('/root/Cycle/HarmoCycle_v1/Result/Oscope/Buettner_gene_cluster3.txt').x.values
temp_adata1 = ref_dataset1[:, gene_oscope1].copy()
len(gene_oscope1)
sc.pp.scale(temp_adata1, max_value=10)
sc.tl.pca(temp_adata1, n_comps=2)

# 计算轮廓系数（基于 'stage' 分组）
sil_score1 = calculate_silhouette_score(temp_adata1, cluster_key='stage')
print(f"Silhouette Score for 'stage': {sil_score1:.3f}")

pd.DataFrame({'X': temp_adata1.obs.index.tolist(),
            'PC1': temp_adata1.obsm['X_pca'][:, 0], 
              'PC2': temp_adata1.obsm['X_pca'][:, 1], 
              'true': temp_adata1.obs['cell_phase'],
              'stage': temp_adata1.obs['stage'],
              'sil_score': [sil_score1] * len(temp_adata1)}).to_csv('/root/Cycle/HarmoCycle_v1/Result/Oscope/Buettner_projection_gene_cluster3.csv')
