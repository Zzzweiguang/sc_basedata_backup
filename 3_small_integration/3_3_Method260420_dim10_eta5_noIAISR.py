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

print("torch version:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("GPU name:", torch.cuda.get_device_name(0))

# %% [markdown]
# In this notebook the annotated datasets are concatinated, a ensembl <-> gene name mapping dictionary is created and different integration methods are implemented/executed on the dataset

# #### integration

# %%
# adata_final = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-aggr-sparse-scranlog1p.h5ad")
adata_final = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_260420/small-concat-aggr-sparse-scranlog1p-noCxG.h5ad")
adata_final = adata_final[adata_final.obs['dataset'] != 'IAISR', :].copy()
# %%
def prepare_adata(adata_final):
    
    #delete all layers that are not needed to reduce memory
    del adata_final.layers['log1p_scran_samplewise'] #only for small integrations
    del adata_final.layers['raw_decontXcounts']
    del adata_final.layers['uncorrected_counts']

    #rename layer to counts
    adata_final.layers['counts'] = adata_final.layers['rounded_corrected_counts']
    del adata_final.layers['rounded_corrected_counts']

    # add "_{dataset}" to the sample name to make them unique
    adata_final.obs["sample"] = adata_final.obs["sample"].astype(str) + "_" + adata_final.obs["dataset"].astype(str)
    adata_final.obs['sample'] = adata_final.obs['sample'].astype('category')   

# %%
def preprocess_adata(adata_final, mode, batch_key):

    # remove unknown and doublets
    print("Removing unknowns and doublets...")
    adata_final = adata_final[~adata_final.obs['cell_type_level1'].isin(['unknown', 'doublets'])].copy()

    
    if mode=="auto":
        print("Using auto mode to do hvg and pca...")
        scib.preprocessing.reduce_data(adata_final, batch_key=batch_key)
        
    elif mode=="manual":
        print("Using manual mode to do hvg and pca...")
        sc.pp.highly_variable_genes(adata_final,  n_top_genes = 2000, batch_key = batch_key)
        sc.pp.pca(adata_final, use_highly_variable=True)

    else:
        raise ValueError("Mode must be 'auto' or 'manual'")

    
    # subset only hvgs
    print("Subset for HVGs...")
    hvg = adata_final.var[adata_final.var['highly_variable']].index.tolist()
    adata_final = adata_final[:, hvg].copy()

    
    #check how many cells have zero counts for all genes
    print("Check how many cells have zero counts for all genes...")
    cellwise_sum = adata_final.X.sum(axis=1)
    num_cells_zero_counts = (cellwise_sum == 0).sum()
    
    if num_cells_zero_counts>0:
        print(num_cells_zero_counts, " cells were found with 0 counts across all genes! Removing these cells now...")
        adata_final = adata_final[cellwise_sum > 0, :]

    
    # check for dublicate gene expressions across cells after HVG selection
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
        "No non unique cells found."

    print("Preprocessing finished!")  

    return adata_final


# %%
def integrate_sca(adata, batch_key, cell_type,):
    
    # adata has to be filtered, hvg filtered and pca applied in obsm[X_pca]
    # log normalized data in .X and counts in .layers["counts"] 

    
    adata_final = adata

    adata_counts = adata_final.copy()
    adata_counts.X = adata_counts.layers["counts"].copy()

    # print("scGen scib implementation")
    # adata_after = scib.ig.scgen(adata_final, batch=batch_key, cell_type=cell_type)
    # adata_final.obsm["scGen"] = adata_after.obsm["corrected_latent"].copy()
    
    '''
    # scGen too slow
    print("Integrating with scGen...")

    adata_final = remove_sparsity(adata_final) # remove sparsity
    adata_final.X = adata_final.X.astype("float32")
    
    epoch = 100
    early_stopping_kwargs = {
        "early_stopping_metric": "val_loss",
        "patience": 25,
        "threshold": 0,
        "reduce_lr": True,
        "lr_patience": 13,
        "lr_factor": 0.1,
    }

    network = sca.models.scgen(adata = adata_final,
                       hidden_layer_sizes=[256,128])

    network.train(n_epochs=epoch, early_stopping_kwargs = early_stopping_kwargs, batch_size = 32)

    adata_after = network.batch_removal(adata_final, batch_key=batch_key, cell_label_key=cell_type,return_latent=True)
    adata_final.obsm["scGen"] = adata_after.obsm["corrected_latent"].copy()
    '''

    
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

    try:
        vae.train(accelerator="gpu", devices=1)
    except TypeError:
        vae.train(use_gpu=True)
    adata_final.obsm["scVI"] = vae.get_latent_representation()

    
    # scANVI
    print("Integrating with scANVI...")

    scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "none", labels_key=cell_type)
    # scanvae.train(max_epochs=40)
    try:
        scanvae.train(max_epochs=40, accelerator="gpu", devices=1)
    except TypeError:
        scanvae.train(max_epochs=40, use_gpu=True)
    adata_final.obsm["scANVI"] = scanvae.get_latent_representation()

    
    
    # scPoli
    print("Integrating with scPoli...")
    early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 10,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
    }
    
    scpoli_model = scPoli(
    adata=adata_final,
    condition_keys=batch_key,
    cell_type_keys=cell_type,
    embedding_dims=10,
    recon_loss='mse',
    )
    
    # scpoli_model.train(
    #     n_epochs=50,
    #     pretraining_epochs=40,
    #     early_stopping_kwargs=early_stopping_kwargs,
    #     eta=6,
    # )
    try:
        scpoli_model.train(
            n_epochs=50,
            pretraining_epochs=40,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=5,
            accelerator="gpu",
            devices=1,
        )
    except TypeError:
        scpoli_model.train(
            n_epochs=50,
            pretraining_epochs=40,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=5,
            use_gpu=True,
        )
    scpoli_latent = scpoli_model.get_latent(adata_final)
    adata_final.obsm["scPoli"] = scpoli_latent


    '''
    ####################################################################################
    # scPoli counts
    print("Integrating with scPoli using raw counts...")
    early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
    }
    
    scpoli_model2 = scPoli(
    adata=adata_counts,
    condition_keys=batch_key,
    cell_type_keys=cell_type,
    embedding_dims=5,
    recon_loss='nb',
    )
    
    scpoli_model2.train(
        n_epochs=50,
        pretraining_epochs=40,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=5,
    )

    scpoli_latent2 = scpoli_model2.get_latent(adata_counts)
    adata_final.obsm["scPoli_counts"] = scpoli_latent2
    ####################################################################################
    '''
    
    # Harmony
    print("Integrating with Harmony...")
    adata_final.obsm["Harmony"] = harmonize(adata_final.obsm["X_pca"], adata_final.obs, batch_key=batch_key)

    
    # LIGER
    print("Integrating with LIGER...")

    batch_cats = adata_final.obs["sample"].cat.categories
    bdata = adata_final.copy()
    # Pyliger normalizes by library size with a size factor of 1
    # So here we give it the count data
    bdata.X = bdata.layers["counts"]
    # List of adata per dataset
    adata_list = [bdata[bdata.obs[batch_key] == b].copy() for b in batch_cats]
    for i, ad in enumerate(adata_list):
        ad.uns["sample_name"] = batch_cats[i]
        # Hack to make sure each method uses the same genes
        ad.uns["var_gene_idx"] = np.arange(bdata.n_vars)
    
    
    liger_data = pyliger.create_liger(adata_list, remove_missing=False, make_sparse=False)
    # Hack to make sure each method uses the same genes
    liger_data.var_genes = bdata.var_names
    pyliger.normalize(liger_data)
    pyliger.scale_not_center(liger_data)
    pyliger.optimize_ALS(liger_data, k=30)
    pyliger.quantile_norm(liger_data)
    
    
    adata_final.obsm["LIGER"] = np.zeros((adata_final.shape[0], liger_data.adata_list[0].obsm["H_norm"].shape[1]))
    for i, b in enumerate(batch_cats):
        adata_final.obsm["LIGER"][adata_final.obs[batch_key] == b] = liger_data.adata_list[i].obsm["H_norm"]


    
    print("Finished integrating the data")
    
    return adata_final

# %%
import time
start_time = time.time()

# %%
prepare_adata(adata_final)

# %%
adata_final

# %%
adata_final = preprocess_adata(adata_final, mode = "auto", batch_key="sample")

# %%
adata_final = integrate_sca(adata_final, batch_key = "sample", cell_type = "cell_type_level1")##no-scgen(版本问题)

# %%
end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")

# %%
# check for dublicates in embeddings

# %%
##no-scgen计算重复行的数量
for embed in ['Harmony', 'LIGER', 'X_pca','scANVI', 'scPoli', 'scVI']:
    diff = adata_final.obsm[embed].shape[0] - np.unique(adata_final.obsm[embed], axis=0).shape[0]
    print(diff," for ", embed )

# %%
#rename obsm X_pca to PCA in adata_final
adata_final.obsm["PCA"] = adata_final.obsm["X_pca"]
del adata_final.obsm["X_pca"]

# %%
# adata_final.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-ensembl-aggr-sparse-method-mse.h5ad") # used in the 3-1_evaluation-HPC.py script
# adata_final.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_1226/small-concat-ensembl-aggr-sparse-method-noCxG-mse.h5ad") # used in the 3-1_evaluation-HPC.py script
adata_final.write_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output_260420/small-concat-ensembl-aggr-sparse-method-noCxG-mse-dim10-eta5-noIAISR.h5ad") # used in the 3-1_evaluation-HPC.py script


