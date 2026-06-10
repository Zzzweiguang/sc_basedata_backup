import pandas as pd
import anndata
import scanpy as sc
import scipy
from scipy.sparse import issparse, coo_matrix, csr_matrix
import os
import numpy as np
import warnings
from scipy import sparse


# 读取已处理好的adata_agg（包含处理后的X矩阵）
print("读取已处理好的adata_agg数据...")
# adata_agg = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-ensembl-aggr-X-py.h5ad")
adata_agg = sc.read_h5ad("./output_260420/small-concat-ensembl-aggr-X-py-noCxG_2.h5ad")


def aggregate_layer(layer_data, var_names):
    """对单个layer进行基因聚合处理"""
    print(f"确保数据为稀疏矩阵格式")
    if not issparse(layer_data):
        layer_data = csr_matrix(layer_data)
    
    print(f"获取基因名称并生成列映射")
    genes = var_names.to_numpy()
    unique_genes, column_mapping = pd.factorize(genes)
    
    print(f"转换为COO格式以便处理坐标")
    layer_coo = layer_data.tocoo()
    
    print(f"生成新的列索引")
    new_columns = unique_genes[layer_coo.col] 
    
    print(f"创建新的COO矩阵")
    aggregated_coo = coo_matrix(
        (layer_coo.data, (layer_coo.row, new_columns)),
        shape=(layer_data.shape[0], len(column_mapping))
    )
    
    print(f"合并重复项并求和")
    aggregated_coo.sum_duplicates()
    
    print(f"转换回高效存储格式")
    aggregated_csr = aggregated_coo.tocsr()
    
    return aggregated_csr


# 定义需要处理的三个layer名称（请替换为实际的层名称）
layers_to_process = ['log1p_scran_samplewise', 'raw_decontXcounts', 'rounded_corrected_counts','uncorrected_counts']  # 这里替换为你的层名称

# 从原始adata中读取需要处理的层数据（假设原始数据路径不变）
print("读取原始数据以获取需要处理的层...")
# adata_original = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-ensembl.h5ad")
adata_original = sc.read_h5ad("./output_260420/small-concat-ensembl-noCxG_2.h5ad")

# 对每个指定的层进行处理并添加到adata_agg
for layer in layers_to_process:
    if layer in adata_original.layers:
        print(f"\n开始处理层: {layer}")
        # 从原始数据中获取该层数据
        original_layer = adata_original.layers[layer]
        # 对层数据进行聚合处理
        aggregated_layer = aggregate_layer(original_layer, adata_original.var_names)
        # 将处理后的层添加到adata_agg
        adata_agg.layers[layer] = aggregated_layer
        print(f"已将处理后的 {layer} 添加到adata_agg")
    else:
        print(f"警告: 层 {layer} 不存在于原始数据中，已跳过")


# 保存包含处理后层数据的结果
print("\n保存处理后的adata_agg(包含处理后的层)...")
# output_path = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-ensembl-aggr.h5ad"
output_path = "./output_260420/small-concat-ensembl-aggr-noCxG_2.h5ad"
adata_agg.write_h5ad(output_path)
print(f"处理完成，文件已保存至: {output_path}")
    