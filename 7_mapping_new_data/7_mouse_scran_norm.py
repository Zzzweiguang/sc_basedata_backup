import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import logging
import os
from scipy.sparse import issparse

# =========================
# rpy2 setup
# =========================
import rpy2.robjects as ro
import rpy2.rinterface_lib.callbacks
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri

# silence R warnings
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

pandas2ri.activate()
anndata2ri.activate()

# =========================
# Load data
# =========================
print("Loading AnnData...")
adata_agg = sc.read_h5ad("./output_all_mouse1226/Mouse_aggregated-gene(human).h5ad")  # ← 改成你的路径

# =========================
# Check zero-count cells
# =========================
print("Check how many cells have zero counts for all genes...")
cellwise_sum = np.array(adata_agg.X.sum(axis=1)).flatten()
num_cells_zero_counts = (cellwise_sum == 0).sum()

if num_cells_zero_counts > 0:
    print(f"{num_cells_zero_counts} cells have 0 counts. Removing...")
    adata_agg = adata_agg[cellwise_sum > 0, :]

adata_final = adata_agg.copy()
print("Remaining cells:", adata_final.n_obs)

# =========================
# Pre-clustering for scran
# =========================
print("Preprocessing for scran clustering...")
adata_pp = adata_final.copy()

sc.pp.normalize_total(adata_pp, target_sum=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, svd_solver="arpack")
sc.pp.neighbors(adata_pp, n_pcs=30)
sc.tl.leiden(adata_pp, key_added="groups", resolution=0.22)

input_groups = adata_pp.obs["groups"]

# =========================
# Prepare matrix for scran
# =========================
print("Preparing matrix for scran...")
data_mat = adata_final.X.T   # genes x cells → cells x genes (scran format)

if not issparse(data_mat):
    raise ValueError("adata_final.X must be a sparse matrix for scran!")

# =========================
# Load R packages
# =========================
print("Loading R packages...")
scran = importr("scran")
Matrix = importr("Matrix")

# =========================
# Pass data to R
# =========================
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = input_groups

# =========================
# Run scran via rpy2
# =========================
print("Running scran::calculateSumFactors ...")

ro.r("""
library(scran)
library(Matrix)

message("Converting to dgCMatrix...")
data_mat <- as(data_mat, "dgCMatrix")

message("Calculating size factors...")
size_factors <- calculateSumFactors(
    data_mat,
    clusters = input_groups,
    min.mean = 0.1
)

message("scran finished.")
""")

# =========================
# Retrieve result
# =========================
size_factors = np.array(ro.r("size_factors")).flatten()
print("Size factors computed:", size_factors.shape)

# =========================
# Save back to AnnData
# =========================
adata_final.obs['size_factors'] = size_factors

print("prepare normalization")
import scipy.sparse as sp
# 检查X是否为COO格式
if isinstance(adata_final.X, sp.coo_matrix):
    # 转换为CSR格式（更适合行访问）
    adata_final.X = adata_final.X.tocsr()
adata_final.write_h5ad("./output_all_mouse1226/Mouse-aggred-prenor.h5ad")

##---------------------------------------------------
print("normalization start")
#Normalize adata 
# adata_final.X /= adata_final.obs['size_factors'].values[:,None]
# sc.pp.log1p(adata_final)
min_threshold = 1e-4  # 根据数据调整，避免除以过小的值
adata_final.obs['size_factors'] = np.clip(
    adata_final.obs['size_factors'],
    a_min=min_threshold,
    a_max=None
)
# 2. 应用归一化并处理异常值
adata_final.X /= adata_final.obs['size_factors'].values[:, None]
# 检查并处理无穷大或NaN值
if np.isinf(adata_final.X.data).any() or np.isnan(adata_final.X.data).any():
    print("检测到异常值，正在处理...")
    # 替换无穷大为极大值，NaN为0
    max_val = np.finfo(np.float64).max
    adata_final.X.data = np.nan_to_num(
        adata_final.X.data,
        nan=0.0,
        posinf=max_val,
        neginf=-max_val
    )
sc.pp.log1p(adata_final)
print("SAVE FILE")
adata_final.write_h5ad("./output_all_mouse1226/Mouse-aggred-normalized.h5ad")
print("ALL DONE.")
