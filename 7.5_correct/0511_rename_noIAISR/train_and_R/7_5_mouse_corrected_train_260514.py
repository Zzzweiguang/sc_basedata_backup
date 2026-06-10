import scanpy as sc
import scarches
from scarches.models.scpoli import scPoli
import scib
import os
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import re
from scipy.sparse import issparse
import scipy.sparse as sparse
import warnings
import torch


print("CUDA available:", torch.cuda.is_available())
print("CUDA device count:", torch.cuda.device_count())

if torch.cuda.is_available():
    print("Using GPU:", torch.cuda.get_device_name(0))
else:
    print("Using CPU")
print(os.getcwd())

adata_final = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/mouse_scpoli_concat_corrected_all.h5ad")
adata_final
#不转密集
def preprocess_adata(adata_final, mode, batch_key):

    def _unique_row_indices(mat):
        if sparse.issparse(mat):
            mat = mat.tocsr(copy=True)
            mat.sum_duplicates()
            mat.sort_indices()
            seen = set()
            unique_indices = []
            for i in range(mat.shape[0]):
                start, end = mat.indptr[i], mat.indptr[i + 1]
                key = (tuple(mat.indices[start:end]), tuple(mat.data[start:end]))
                if key not in seen:
                    seen.add(key)
                    unique_indices.append(i)
            return np.array(unique_indices, dtype=int)
        else:
            mat = np.asarray(mat)
            _, unique_indices = np.unique(mat, axis=0, return_index=True)
            return np.sort(unique_indices)

    adata_final = adata_final.copy()

    adata_final.obs[batch_key] = adata_final.obs[batch_key].astype("category")
    adata_final.obs[batch_key] = adata_final.obs[batch_key].cat.remove_unused_categories()

    if "counts" not in adata_final.layers:
        if "rounded_corrected_counts" in adata_final.layers:
            adata_final.layers["counts"] = adata_final.layers["rounded_corrected_counts"].copy()
        else:
            raise ValueError("No counts layer found.")

    # 关键：从 counts 重新构建干净的 log-normalized .X
    print("Rebuilding adata.X from counts...")
    adata_final.X = adata_final.layers["counts"].copy()
    sc.pp.normalize_total(adata_final, target_sum=1e4)
    sc.pp.log1p(adata_final)

    if mode == "auto":
        print("Using auto mode to do hvg and pca...")

        hvg_list = scib.preprocessing.hvg_batch(
            adata_final,
            batch_key=batch_key,
            target_genes=2000,
            flavor="seurat"
        )

        adata_final.var["highly_variable"] = adata_final.var_names.isin(hvg_list)

        sc.pp.pca(
            adata_final,
            use_highly_variable=True,
            n_comps=50
        )

    elif mode == "manual":
        print("Using manual mode to do hvg and pca...")

        sc.pp.highly_variable_genes(
            adata_final,
            n_top_genes=2000,
            batch_key=batch_key,
            flavor="seurat"
        )

        sc.pp.pca(
            adata_final,
            use_highly_variable=True,
            n_comps=50
        )

    else:
        raise ValueError("Mode must be 'auto' or 'manual'")

    print("Subset for HVGs...")
    hvg = adata_final.var[adata_final.var["highly_variable"]].index.tolist()
    adata_final = adata_final[:, hvg].copy()

    print("Check how many cells have zero counts for all genes...")
    cellwise_sum = np.asarray(adata_final.X.sum(axis=1)).ravel()
    num_cells_zero_counts = int((cellwise_sum == 0).sum())

    if num_cells_zero_counts > 0:
        print(num_cells_zero_counts, "cells were found with 0 expression across all HVGs! Removing these cells now...")
        adata_final = adata_final[cellwise_sum > 0, :].copy()

    print("Check for duplicate gene expressions")

    counts = adata_final.layers["counts"]
    unique_indices_sorted = _unique_row_indices(counts)
    diff = counts.shape[0] - len(unique_indices_sorted)

    if diff > 0:
        print(diff, "non-unique cell expression profiles found! Removing them...")
        adata_final = adata_final[unique_indices_sorted, :].copy()
    else:
        print("No non-unique cells found.")

    print("Adding atlas key")
    adata_final.obs["atlas_key"] = [
        "query" if cell_type == "unknown" else "ref"
        for cell_type in adata_final.obs["cell_type_level1"]
    ]

    print("Preprocessing finished!")

    return adata_final
adata_pp = preprocess_adata(adata_final, mode = "auto", batch_key="sample")
adata_pp.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/mouse_scpoli_concat_corrected_all_pp.h5ad")
adata_pp

adata_pp.X = adata_pp.X.astype("float32")
print(adata_pp.obs["cell_type_level1"].value_counts())

print("Train reference")
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
    
scpoli_model = scPoli(
adata=adata_pp,
condition_keys="sample",
cell_type_keys="cell_type_level1_corrected",
embedding_dims=15,
recon_loss='mse',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=10,
    use_gpu=True
)
print("save model")
scpoli_model.save("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/all_mouse_model_corrected/", overwrite=True)

print("get latent")
scpoli_model.model.eval()
data_latent = scpoli_model.get_latent(
    adata_pp,
    mean=True
)

adata_latent_full = sc.AnnData(data_latent)
adata_latent_full.obs = adata_pp.obs.copy()

adata_latent_full.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/mouse_scpoli_concat_corrected_all_preumap.h5ad")
print("umap!!!")
# adata_latent_full = adata_latent_full[adata_latent_full.obs["dataset"] != "IAISR"].copy()
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
adata_latent_full.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/mouse_scpoli_concat_corrected_all_no_gene.h5ad")

print("finished!!!")