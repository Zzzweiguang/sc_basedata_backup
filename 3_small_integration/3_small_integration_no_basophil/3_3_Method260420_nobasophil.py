# %%
import warnings
# 忽略所有 DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning, module="scib|scgen|typing_extensions")

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import os
from scipy.sparse import issparse
import scvi
import scib
from scarches.models.scpoli import scPoli
import scarches as sca
from harmony import harmonize
import pyliger
from scarches.dataset.trvae.data_handling import remove_sparsity
import dask.dataframe as dd
import time
import torch
import gc

# %%
# GPU check
USE_GPU = torch.cuda.is_available()
print("CUDA available:", USE_GPU)

if USE_GPU:
    print("GPU:", torch.cuda.get_device_name(0))

# %% [markdown]
# In this notebook the annotated datasets are concatinated, a ensembl <-> gene name mapping dictionary is created and different integration methods are implemented/executed on the dataset

# #### integration

# %%
# adata_final = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-aggr-sparse-scranlog1p.h5ad")
adata_final = sc.read_h5ad(
    "/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_260420/small-concat-aggr-sparse-scranlog1p-noCxG_3.h5ad"
)
adata_final = adata_final[adata_final.obs['cell_type_level1'] != 'Basophil'].copy()
print(adata_final.obs['cell_type_level1'].value_counts())

# %%
def prepare_adata(adata_final):
    
    # delete all layers that are not needed to reduce memory
    del adata_final.layers['log1p_scran_samplewise']  # only for small integrations
    del adata_final.layers['raw_decontXcounts']
    del adata_final.layers['uncorrected_counts']

    # rename layer to counts
    adata_final.layers['counts'] = adata_final.layers['rounded_corrected_counts']
    del adata_final.layers['rounded_corrected_counts']

    # add "_{dataset}" to the sample name to make them unique
    adata_final.obs["sample"] = (
        adata_final.obs["sample"].astype(str)
        + "_"
        + adata_final.obs["dataset"].astype(str)
    )
    adata_final.obs['sample'] = adata_final.obs['sample'].astype('category')   

# %%
def preprocess_adata(adata_final, mode, batch_key):

    # remove unknown and doublets
    print("Removing unknowns and doublets...")
    adata_final = adata_final[
        ~adata_final.obs['cell_type_level1'].isin(['unknown', 'doublets'])
    ].copy()

    if mode == "auto":
        print("Using auto mode to do hvg and pca...")
        scib.preprocessing.reduce_data(adata_final, batch_key=batch_key)
        
    elif mode == "manual":
        print("Using manual mode to do hvg and pca...")
        sc.pp.highly_variable_genes(
            adata_final,
            n_top_genes=2000,
            batch_key=batch_key,
        )
        sc.pp.pca(adata_final, use_highly_variable=True)

    else:
        raise ValueError("Mode must be 'auto' or 'manual'")

    # subset only hvgs
    print("Subset for HVGs...")
    hvg = adata_final.var[adata_final.var['highly_variable']].index.tolist()
    adata_final = adata_final[:, hvg].copy()

    # check how many cells have zero counts for all genes
    print("Check how many cells have zero counts for all genes...")
    cellwise_sum = adata_final.X.sum(axis=1)
    num_cells_zero_counts = (cellwise_sum == 0).sum()
    
    if num_cells_zero_counts > 0:
        print(
            num_cells_zero_counts,
            " cells were found with 0 counts across all genes! Removing these cells now...",
        )
        adata_final = adata_final[cellwise_sum > 0, :]

    # check for duplicate gene expressions across cells after HVG selection
    print("Check for dublicate gene expressions")
    before = adata_final.layers["counts"].toarray().shape[0]
    after = np.unique(adata_final.layers["counts"].toarray(), axis=0).shape[0]
    diff = before - after

    if diff > 0:
        print(diff, " Non unique cell expression profiles found! Removing them...")
        
        counts_array = adata_final.layers["counts"].toarray()

        # Find the indices of unique rows
        _, unique_indices = np.unique(counts_array, axis=0, return_index=True)

        # Sort the indices to maintain the original order
        unique_indices_sorted = np.sort(unique_indices)

        # Filter the AnnData object to keep only unique rows
        adata_final = adata_final[unique_indices_sorted]
    else:
        print("No non unique cells found.")

    print("Preprocessing finished!")  

    return adata_final

# %%
def integrate_fixed_methods(adata, batch_key, cell_type):
    
    adata_final = adata

    adata_counts = adata_final.copy()
    adata_counts.X = adata_counts.layers["counts"].copy()

    # scVI
    print("Integrating with scVI...")

    sca.models.SCVI.setup_anndata(adata_counts, batch_key=batch_key)

    vae = sca.models.SCVI(
        adata_counts,
        n_layers=2,
        encode_covariates=True,
        deeply_inject_covariates=False,
        use_layer_norm="both",
        use_batch_norm="none",
    )

    vae.train(
        accelerator="gpu" if USE_GPU else "cpu",
        devices=1 if USE_GPU else "auto",
    )

    adata_final.obsm["scVI"] = vae.get_latent_representation()

    # scANVI
    print("Integrating with scANVI...")

    scanvae = sca.models.SCANVI.from_scvi_model(
        vae,
        unlabeled_category="none",
        labels_key=cell_type,
    )

    scanvae.train(
        max_epochs=40,
        accelerator="gpu" if USE_GPU else "cpu",
        devices=1 if USE_GPU else "auto",
    )

    adata_final.obsm["scANVI"] = scanvae.get_latent_representation()

    # clean GPU cache after scVI/scANVI
    if USE_GPU:
        torch.cuda.empty_cache()

    # Harmony
    print("Integrating with Harmony...")
    adata_final.obsm["Harmony"] = harmonize(
        adata_final.obsm["X_pca"],
        adata_final.obs,
        batch_key=batch_key,
    )

    # LIGER
    print("Integrating with LIGER...")

    batch_cats = adata_final.obs["sample"].cat.categories
    bdata = adata_final.copy()
    bdata.X = bdata.layers["counts"]

    adata_list = [
        bdata[bdata.obs[batch_key] == b].copy()
        for b in batch_cats
    ]

    for i, ad in enumerate(adata_list):
        ad.uns["sample_name"] = batch_cats[i]
        ad.uns["var_gene_idx"] = np.arange(bdata.n_vars)
    
    liger_data = pyliger.create_liger(
        adata_list,
        remove_missing=False,
        make_sparse=False,
    )

    liger_data.var_genes = bdata.var_names

    pyliger.normalize(liger_data)
    pyliger.scale_not_center(liger_data)
    pyliger.optimize_ALS(liger_data, k=30)
    pyliger.quantile_norm(liger_data)
    
    adata_final.obsm["LIGER"] = np.zeros(
        (
            adata_final.shape[0],
            liger_data.adata_list[0].obsm["H_norm"].shape[1],
        )
    )

    for i, b in enumerate(batch_cats):
        adata_final.obsm["LIGER"][adata_final.obs[batch_key] == b] = (
            liger_data.adata_list[i].obsm["H_norm"]
        )

    print("Finished fixed integrations")
    
    return adata_final

# %%
def integrate_scpoli_only(
    adata,
    batch_key,
    cell_type,
    scpoli_embedding_dims=10,
    scpoli_eta=5,
    scpoli_n_epochs=50,
    scpoli_pretraining_epochs=40,
    scpoli_recon_loss="mse",
):
    
    print(
        "Integrating with scPoli..."
        f" dim={scpoli_embedding_dims}, eta={scpoli_eta}, "
        f"n_epochs={scpoli_n_epochs}, pretraining_epochs={scpoli_pretraining_epochs}"
    )

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
        adata=adata,
        condition_keys=batch_key,
        cell_type_keys=cell_type,
        embedding_dims=scpoli_embedding_dims,
        recon_loss=scpoli_recon_loss,
    )
    
    # Different scarches/scPoli versions handle GPU arguments differently.
    # First try explicit GPU usage; if that version does not support use_gpu, fall back safely.
    try:
        scpoli_model.train(
            n_epochs=scpoli_n_epochs,
            pretraining_epochs=scpoli_pretraining_epochs,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=scpoli_eta,
            use_gpu=True if USE_GPU else False,
        )
    except TypeError:
        print(
            "scPoli train() does not support use_gpu argument in this version. "
            "Training with default device."
        )
        scpoli_model.train(
            n_epochs=scpoli_n_epochs,
            pretraining_epochs=scpoli_pretraining_epochs,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=scpoli_eta,
        )

    scpoli_latent = scpoli_model.get_latent(adata)
    adata.obsm["scPoli"] = scpoli_latent

    if USE_GPU:
        torch.cuda.empty_cache()
    
    return adata

# %%
start_time = time.time()

prepare_adata(adata_final)

adata_final = preprocess_adata(
    adata_final,
    mode="auto",
    batch_key="sample",
)

# 固定方法只跑一次
adata_final = integrate_fixed_methods(
    adata_final,
    batch_key="sample",
    cell_type="cell_type_level1",
)

# rename obsm X_pca to PCA in adata_final
adata_final.obsm["PCA"] = adata_final.obsm["X_pca"]
del adata_final.obsm["X_pca"]

scpoli_grid = [
    {"embedding_dims": 5,  "eta": 2},
    {"embedding_dims": 5,  "eta": 5},
    {"embedding_dims": 5,  "eta": 10},

    {"embedding_dims": 10, "eta": 2},
    {"embedding_dims": 10, "eta": 5},
    {"embedding_dims": 10, "eta": 10},
    {"embedding_dims": 10, "eta": 15},

    {"embedding_dims": 15, "eta": 2},
    {"embedding_dims": 15, "eta": 5},
    {"embedding_dims": 15, "eta": 10},
    {"embedding_dims": 15, "eta": 15},

    {"embedding_dims": 20, "eta": 5},
    {"embedding_dims": 20, "eta": 10},
]

output_dir = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/3_small_integration_no_basophil/output_nobasophil"

os.makedirs(output_dir, exist_ok=True)

for params in scpoli_grid:
    dim = params["embedding_dims"]
    eta = params["eta"]

    output_file = os.path.join(
        output_dir,
        f"small-concat-ensembl-aggr-sparse-method-noIAISR-mse-dim{dim}-eta{eta}.h5ad",
    )

    if os.path.exists(output_file):
        print("File already exists, skipping:", output_file)
        continue

    print("=" * 80)
    print(f"Running scPoli parameters: dim={dim}, eta={eta}")
    print("=" * 80)

    adata_run = adata_final.copy()

    adata_run = integrate_scpoli_only(
        adata_run,
        batch_key="sample",
        cell_type="cell_type_level1",
        scpoli_embedding_dims=dim,
        scpoli_eta=eta,
        scpoli_n_epochs=50,
        scpoli_pretraining_epochs=40,
        scpoli_recon_loss="mse",
    )

    for embed in ["Harmony", "LIGER", "PCA", "scANVI", "scPoli", "scVI"]:
        diff = adata_run.obsm[embed].shape[0] - np.unique(
            adata_run.obsm[embed],
            axis=0,
        ).shape[0]
        print(diff, "for", embed)

    adata_run.write_h5ad(output_file)
    print("Saved:", output_file)

    del adata_run
    gc.collect()

    if USE_GPU:
        torch.cuda.empty_cache()

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")