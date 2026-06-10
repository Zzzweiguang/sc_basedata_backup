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

adata_final = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_allhuman/work_0513_2/scPoli_concat_corrected_all.h5ad")
adata_final
def preprocess_adata(adata_final, mode, batch_key, copy=False):
    import gc
    import hashlib
    import numpy as np
    from scipy import sparse

    def _unique_row_indices(mat):
        if sparse.issparse(mat):
            mat = mat.tocsr(copy=False)
            mat.sum_duplicates()
            mat.sort_indices()
            seen = set()
            unique_indices = []
            for i in range(mat.shape[0]):
                start, end = mat.indptr[i], mat.indptr[i + 1]
                h = hashlib.blake2b(digest_size=16)
                h.update(mat.indices[start:end].tobytes())
                h.update(mat.data[start:end].tobytes())
                key = h.digest()
                if key not in seen:
                    seen.add(key)
                    unique_indices.append(i)
            return np.asarray(unique_indices, dtype=np.int64)
        else:
            mat = np.asarray(mat)
            _, unique_indices = np.unique(mat, axis=0, return_index=True)
            return np.sort(unique_indices)

    if copy:
        adata_final = adata_final.copy()

    adata_final.obs[batch_key] = adata_final.obs[batch_key].astype("category").cat.remove_unused_categories()

    if "counts" not in adata_final.layers:
        if "rounded_corrected_counts" in adata_final.layers:
            adata_final.layers["counts"] = adata_final.layers["rounded_corrected_counts"]
        else:
            raise ValueError("No counts layer found.")

    print("Rebuilding adata.X from counts...")
    X = adata_final.layers["counts"]
    if sparse.issparse(X):
        adata_final.X = X.copy().astype(np.float32)
    else:
        adata_final.X = np.asarray(X, dtype=np.float32).copy()

    sc.pp.normalize_total(adata_final, target_sum=1e4)
    sc.pp.log1p(adata_final)

    if mode == "auto":
        print("Using auto mode to do hvg...")
        hvg_list = scib.preprocessing.hvg_batch(adata_final, batch_key=batch_key, target_genes=2000, flavor="seurat")
        adata_final.var["highly_variable"] = adata_final.var_names.isin(hvg_list)
    elif mode == "manual":
        print("Using manual mode to do hvg...")
        sc.pp.highly_variable_genes(adata_final, n_top_genes=2000, batch_key=batch_key, flavor="seurat")
    else:
        raise ValueError("Mode must be 'auto' or 'manual'")

    print("Subset for HVGs before PCA...")
    hvg = adata_final.var_names[adata_final.var["highly_variable"]].tolist()
    adata_final = adata_final[:, hvg].copy()
    gc.collect()

    print("Running PCA...")
    sc.pp.pca(adata_final, n_comps=50)

    print("Check how many cells have zero counts for all genes...")
    cellwise_sum = np.asarray(adata_final.X.sum(axis=1)).ravel()
    keep = cellwise_sum > 0
    num_cells_zero_counts = int((~keep).sum())

    if num_cells_zero_counts > 0:
        print(num_cells_zero_counts, "cells were found with 0 expression across all HVGs! Removing these cells now...")
        adata_final = adata_final[keep, :].copy()
        gc.collect()

    print("Check for duplicate gene expressions")
    counts = adata_final.layers["counts"]
    unique_indices_sorted = _unique_row_indices(counts)
    diff = counts.shape[0] - len(unique_indices_sorted)

    if diff > 0:
        print(diff, "non-unique cell expression profiles found! Removing them...")
        adata_final = adata_final[unique_indices_sorted, :].copy()
        gc.collect()
    else:
        print("No non-unique cells found.")

    print("Adding atlas key")
    adata_final.obs["atlas_key"] = np.where(
        adata_final.obs["cell_type_level1"].astype(str).values == "unknown",
        "query",
        "ref"
    )

    print("Preprocessing finished!")
    return adata_final
adata_pp = preprocess_adata(adata_final, mode = "auto", batch_key="sample")
adata_pp.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_allhuman/work_0513_2/scPoli_concat_corrected_all_pp.h5ad")###合并的时候改的dc名字
adata_pp
# adata_pp = sc.read_h5ad("./output_query_human/all_human_corrected_counts_no_IAISR_pp_rename.h5ad")
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
scpoli_model.save("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_mouse_scpoli/work_0513/all_human_corected", overwrite=True)##在小鼠

print("get latent")
scpoli_model.model.eval()
data_latent = scpoli_model.get_latent(
    adata_pp,
    mean=True
)

adata_latent_full = sc.AnnData(data_latent)
adata_latent_full.obs = adata_pp.obs.copy()
adata_latent_full = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_allhuman/work_0513_2/scPoli_concat_corrected_all_preumap.h5ad")
print("umap!!!")
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
adata_latent_full.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/0511_rename_noIAISR/output_allhuman/work_0513_2/scPoli_concat_corrected_all_no_gene.h5ad")

print("finished!!!")