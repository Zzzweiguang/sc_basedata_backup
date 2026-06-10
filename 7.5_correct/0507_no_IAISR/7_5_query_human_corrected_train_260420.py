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


# adata_pp = sc.read_h5ad("./output_260420/all_human_corrected_counts_train_pp.h5ad")
# adata_pp = sc.read_h5ad("./output_260420/all_human_corrected_counts_train_pp_2.h5ad")
# adata_pp = sc.read_h5ad("./output_query_human/all_human_corrected_counts_no_IAISR_pp.h5ad")
adata_pp = sc.read_h5ad("./output_query_human/all_human_corrected_counts_no_IAISR_pp_rename.h5ad")


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
cell_type_keys="cell_type_level1",
embedding_dims=20,
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
scpoli_model.save("./output_query_human/model_corrected/all_human", overwrite=True)

print("get latent")
scpoli_model.model.eval()
data_latent = scpoli_model.get_latent(
    adata_pp,
    mean=True
)

adata_latent_full = sc.AnnData(data_latent)
adata_latent_full.obs = adata_pp.obs.copy()

# adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_preumap_corrected.h5ad")###之前校正1.0版本的好像覆盖了。。。
# adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_preumap_corrected_3.h5ad")###之前校正1.0版本的好像覆盖了。。。
adata_latent_full.write_h5ad("./output_query_human/Altas_level1_all_human_preumap_corrected_no_IAISR_renamed.h5ad")###这个是修正后+transfer后+矫正后
adata_latent_full = sc.read_h5ad("./output_query_human/Altas_level1_all_human_preumap_corrected_no_IAISR_renamed.h5ad")###这个是修正后+transfer后+矫正后
print("umap!!!")
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
# adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_umap_nogene_noIAISR_corrected.h5ad")
# adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_umap_nogene_noIAISR_corrected_2.h5ad")
adata_latent_full.write_h5ad("./output_query_human/Altas_level1_all_human_umap_nogene_noIAISR_corrected_no_IAISR_renamed.h5ad")###这个是想看校正后的数据的10维

print("finished!!!")