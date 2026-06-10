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

adata_query= sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7_mapping_new_data/output_human_as_atlas/mergerd_aggred-log-normalized-hvg.h5ad")
scpoli_model = scPoli.load("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7_mapping_new_data/model/reference_level1", adata=adata_query)
print(scpoli_model)

scpoli_query = scPoli.load_query_data(
    adata=adata_query,
    reference_model=scpoli_model,
    labeled_indices=[],
)
print("start training")
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)

print("saving")
scpoli_query.save("/home/lixiangyu/zr/Annotate/ANNOTATE_new/7_mapping_new_data/model/human_as_atlas_query-level1", overwrite=True)
print("finished!")