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
import torch

print("CUDA available:", torch.cuda.is_available())
print("CUDA device count:", torch.cuda.device_count())

if torch.cuda.is_available():
    print("Using GPU:", torch.cuda.get_device_name(0))
else:
    print("Using CPU")
print(os.getcwd())

adata_query= sc.read_h5ad("./output_mouse/Mouse-aggred-normalized-hvg.h5ad")
scpoli_model = scPoli.load("/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/train_ref_model_noBasophil_0521/output/ref_model_noBasophil", adata=adata_query)##没有correct的
print("start train query model")
scpoli_query = scPoli.load_query_data(
    adata=adata_query,
    reference_model=scpoli_model,
    labeled_indices=[],
)
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10,
    use_gpu=True,
)

print("save query model")
scpoli_query.save("./output_mouse/mouse_model", overwrite=True)

print("Create embedding adata")
adata_query = sc.read_h5ad("./output_mouse/Mouse-aggred-normalized-hvg.h5ad")
scpoli_query = scPoli.load("./output_mouse/mouse_model", adata=adata_query)
adata_query.X = adata_query.X.astype("float32")
results_dict = scpoli_query.classify(adata_query, scale_uncertainties=True)
adata = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/4_big_integration/train_ref_model_noBasophil_0521/output/big-concat-ensembl-nodub-pp.h5ad")
#get whole atlas including reference

#get latent representation of reference data
scpoli_query.model.eval()
data_latent_source = scpoli_query.get_latent(
    adata,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = adata.obs.copy()

#get latent representation of query data
data_latent= scpoli_query.get_latent(
    adata_query,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = adata_query.obs.copy()

#get label annotations
adata_latent.obs['cell_type_pred'] = results_dict['cell_type_level1']['preds'].tolist()
adata_latent.obs['cell_type_uncert'] = results_dict['cell_type_level1']['uncert'].tolist()
#adata_latent.obs['classifier_outcome'] = (adata_latent.obs['cell_type_pred'] == adata_latent.obs['cell_type_level1'])

#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent],
    batch_key='query'
)
adata_latent_full.obs['cell_type_pred_ref'] = np.where(
    adata_latent_full.obs['query'].isin(['0']),  # Condition
    "Reference",                                # Value if condition is true
    adata_latent_full.obs['cell_type_pred']     # Value if condition is false
)
print(adata_latent_full)
print(adata_latent_full.obs)
print(adata_latent_full.obs["cell_type_level1"])
# Create 'cell_type_level1_combined' based on the 'query' condition
adata_latent_full.obs['cell_type_level1_human(measure)&mouse'] = np.where(
    adata_latent_full.obs['query'] == "0",
    adata_latent_full.obs['cell_type_level1'],
    adata_latent_full.obs['cell_type_pred_ref']
)
print(adata_latent_full.obs)
print(adata_latent_full.obs["cell_type_level1_human(measure)&mouse"])
print(adata_latent_full)
print("save pre umap")
adata_latent_full.write_h5ad("./output_mouse/human_mouse-aggred-normalized-hvg-nogene-preumap.h5ad")
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
print("save umap ed")
adata_latent_full.write_h5ad("./output_mouse/human_mouse-aggred-normalized-hvg-nogene-umap.h5ad")
print("finished!!!")
