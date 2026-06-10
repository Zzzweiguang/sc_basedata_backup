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


adata = sc.read_h5ad("./output_1226/Atlas_level1-corrected-hvg-preint.h5ad") # for new pipeline with no level 2 with new pdc
print("prepare.....")
# Subset where atlas_key is 'ref'
adata_ref = adata[adata.obs['atlas_key'] == 'ref'].copy()
# Subset where atlas_key is 'query'
adata_query = adata[adata.obs['atlas_key'] == 'query'].copy()
#adata_ref.obs["cell_type_level2"].value_counts()
print(adata_ref.obs["cell_type_level1_corrected"].value_counts())
#adata_query.obs["cell_type_level2"].value_counts()
print(adata_query.obs["cell_type_level1_corrected"].value_counts())
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
#cell_type_keys="cell_type_level2",
cell_type_keys="cell_type_level1_corrected", # for new pipeline
embedding_dims=10,
recon_loss='mse',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
print("save ref model")
scpoli_model.save("./output_1226/model_corrected/ref-level1_corrected", overwrite=True) 
print("Train with query")
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
print("save query model")
scpoli_query.save("./output_1226/model_corrected/query-level1_corrected", overwrite=True) # with new pdcs
print("Reload")
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
cell_type_keys="cell_type_level1_corrected",
embedding_dims=10,
recon_loss='mse',
)
scpoli_query = scPoli.load_query_data(
    adata=adata_query,
    reference_model=scpoli_model,
    labeled_indices=[],
)
adata_query.X = adata_query.X.astype("float32")
adata_ref.X = adata_ref.X.astype("float32")
scpoli_loaded = scpoli_query.load(dir_path="./output_1226/model_corrected/query-level1_corrected", adata=adata_query)
results_dict = scpoli_loaded.classify(adata_query, scale_uncertainties=True)
#get latent representation of reference data
scpoli_loaded.model.eval()
data_latent_source = scpoli_loaded.get_latent(
    adata_ref,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = adata_ref.obs.copy()

#get latent representation of query data
data_latent= scpoli_loaded.get_latent(
    adata_query,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = adata_query.obs.copy()

#get label annotations
adata_latent.obs['cell_type_pred'] = results_dict['cell_type_level1_corrected']['preds'].tolist()
adata_latent.obs['cell_type_uncert'] = results_dict['cell_type_level1_corrected']['uncert'].tolist()

#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent],
    batch_key='query'
)
adata_latent_full.obs['cell_type_pred'][adata_latent_full.obs['query'].isin(['0'])] = "Reference" #np.nan
print("save data before umap")
adata_latent_full.write("./output_1226/Atlas_level1-corrected-hvg-preumap.h5ad")
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
adata_latent_full.obs['cell_type_level1_all'] = [
    row['cell_type_level1_corrected'] if row['cell_type_pred'] == 'Reference' else row['cell_type_pred'] 
    for _, row in adata_latent_full.obs.iterrows()
]
print("save data after umap")
adata_latent_full.write_h5ad("./output_1226/Atlas_level1-corrected-hvg-integrated.h5ad")

print("finish without error!!!!!!")