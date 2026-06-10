import scanpy as sc
import scarches
from scarches.models.scpoli import scPoli
import scib
import os
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import scipy.sparse as sparse
import warnings
print(os.getcwd())

# Filter out DeprecationWarnings
#warnings.filterwarnings("ignore", category=DeprecationWarning)


adata_final = sc.read_h5ad("./output_1226/big-concat-ensembl-nodub-aggr-scranlog1p-annot-nodoublets.h5ad")
print(adata_final)
def preprocess_adata(adata_final, mode, batch_key):

    #rename layer to counts
    adata_final.layers['counts'] = adata_final.layers['rounded_corrected_counts']
    del adata_final.layers['rounded_corrected_counts']
    
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

    
    #check how many cells have zero "counts" for all genes
    #calculated on .X, which is log-normalized (not counts)
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

    
    # add atlas key (reference or query dataset)
    # unset roche cell_type_level1 obs
    
    print("Adding atlas key")
    adata_final.obs["atlas_key"] = ["query" if cell_type == "unknown" else "ref" for cell_type in adata_final.obs["cell_type_level1"]]
    
    print("Preprocessing finished!")  

    return adata_final
print("preprocess")
adata_pp = preprocess_adata(adata_final, mode = "auto", batch_key="sample")
adata_pp.write_h5ad("./output_1226/big-concat-ensembl-nodub-aggr-scranlog1p-annot-nodoublets-pp.h5ad")
print(adata_pp)
print("Integration with scPoli")
adata = sc.read_h5ad("./output_1226/big-concat-ensembl-nodub-aggr-scranlog1p-annot-nodoublets-pp.h5ad")
# Subset where atlas_key is 'ref'
adata_ref = adata[adata.obs['atlas_key'] == 'ref'].copy()

# Subset where atlas_key is 'query'
adata_query = adata[adata.obs['atlas_key'] == 'query'].copy()
print(adata_ref)
adata_query ##unknown
print(adata_ref.obs["cell_type_level1"].value_counts())
print(adata_query.obs["cell_type_level1"].value_counts())


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
adata=adata_ref,
condition_keys="sample",
cell_type_keys="cell_type_level1",
embedding_dims=10,
recon_loss='mse',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
scpoli_model.save("./output_1226/model_percorrect/refmse", overwrite=True)

print("Reference mapping of unlabeled query dataset")
scpoli_query = scPoli.load_query_data(
    adata=adata_query,
    reference_model=scpoli_model,
    labeled_indices=[],
)
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)
scpoli_query.save("./output_1226/model_percorrect/querymse", overwrite=True)