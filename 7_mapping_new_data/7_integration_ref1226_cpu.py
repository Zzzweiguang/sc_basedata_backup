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
print(os.getcwd())


adata = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7_mapping_new_data/output_all_human1226/Atlas_reference_preint-hvg.h5ad")
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
adata=adata,
condition_keys="sample",
cell_type_keys="cell_type_level1", 
embedding_dims=10,
recon_loss='mse',
)

scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=4
)
print("save ref model")
scpoli_model.save("./output_all_human1226/model_1226/reference_level1", overwrite=True) # new pipeline with pdcs
print("finished!!!!")