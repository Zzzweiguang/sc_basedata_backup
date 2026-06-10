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
import torch

print("CUDA available:", torch.cuda.is_available())
print("CUDA device count:", torch.cuda.device_count())

if torch.cuda.is_available():
    print("Using GPU:", torch.cuda.get_device_name(0))
else:
    print("Using CPU")
print(os.getcwd())

adata_pp = sc.read_h5ad("./output_260420/big-concat-ensembl-nodub-aggr-scranlog1p-annot-nodoublets-pp_3.h5ad")
print(adata_pp)
print("Integration with scPoli")
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
cell_type_keys="cell_type_level1",
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
scpoli_model.save("./output_260420/model_percorrect/all_measure_no_IAISR", overwrite=True)


scpoli_model.model.eval()

data_latent = scpoli_model.get_latent(
    adata_pp,
    mean=True
)

adata_latent_full = sc.AnnData(data_latent)
adata_latent_full.obs = adata_pp.obs.copy()
adata_latent_full.write_h5ad("./output_260420/Altas_level1_measure_preumap_2.h5ad")
# adata_latent_full.write_h5ad("./output_260420/Altas_level1_measure_preumap.h5ad")

adata_latent_full = sc.read_h5ad("./output_260420/Altas_level1_measure_preumap_2.h5ad")
adata_latent_full

sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)

adata_latent_full.write_h5ad("./output_260420/Altas_level1_measure_nogene_noIAISR_2.h5ad")
print("finished")