import pandas as pd
import anndata
import scanpy as sc
import scipy  # 导入 scipy 库
from scipy.sparse import issparse, hstack
import scib
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy import sparse


adata = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/output/big-concat-ensembl-nodub.h5ad")

##查看细胞名称是否重复
unique_obs, counts = np.unique(adata.obs_names, return_counts=True)
# Filter the observation names with more than one count
duplicate_obs = unique_obs[counts > 1]
print(duplicate_obs)
# 修改细胞名，确保唯一 
adata.obs_names_make_unique()

adata.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/output/big-concat-ensembl-nodub.h5ad")

def aggregate_genes_sparse(adata):
    print(f"确保数据为稀疏矩阵格式")
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    
    print(f"获取基因名称并生成列映射")
    genes = adata.var_names.to_numpy()
    unique_genes, column_mapping = pd.factorize(genes)
    
    print(f"转换为COO格式以便处理坐标")
    X_coo = adata.X.tocoo()
    
    print(f"生成新的列索引")
    new_columns = unique_genes[X_coo.col] 
    
    print(f"创建新的COO矩阵")
    aggregated_coo = sparse.coo_matrix(
        (X_coo.data, (X_coo.row, new_columns)),
        shape=(adata.n_obs, len(column_mapping))  # 使用唯一基因数量作为列数
    )
    
    print(f"合并重复项并求和")
    aggregated_coo.sum_duplicates()
    
    print(f"转换回高效存储格式")
    aggregated_csr = aggregated_coo.tocsr()
    
    print(f"构建新的var元数据_保留首次出现的基因信息")
    unique_var = adata.var[~adata.var.index.duplicated(keep='first')]
    
    print(f"创建聚合后的AnnData对象")
    adata_agg = anndata.AnnData(
        X=aggregated_csr,
        obs=adata.obs,
        var=unique_var,
        dtype=aggregated_csr.dtype
    )
    
    return adata_agg

adata_agg = aggregate_genes_sparse(adata)
print(f"保存文件")
adata_agg.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/output/big-concat-ensembl-nodub-aggr-X.h5ad")

