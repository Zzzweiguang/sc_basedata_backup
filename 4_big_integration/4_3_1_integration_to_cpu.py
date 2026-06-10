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


adata = sc.read_h5ad("./output_1226/big-concat-ensembl-nodub-aggr-scranlog1p-annot-nodoublets-pp.h5ad")
# Subset where atlas_key is 'ref'
adata_ref = adata[adata.obs['atlas_key'] == 'ref'].copy()

# Subset where atlas_key is 'query'
adata_query = adata[adata.obs['atlas_key'] == 'query'].copy()


print("load model....")
scpoli_model = scPoli.load(
    dir_path="./output_1226/model_percorrect/refmse",
    adata=adata_ref
)

scpoli_query = scPoli.load(
    dir_path="./output_1226/model_percorrect/querymse",
    adata=adata_query
)
print("transfer to cpu")
scpoli_model.model.to("cpu")
scpoli_query.model.to("cpu")


print("save model_cpu....")
scpoli_model.save(
    "./output_1226/model_percorrect/refmse_cpu",
    overwrite=True
)
scpoli_query.save(
    "./output_1226/model_percorrect/querymse_cpu",
    overwrite=True
)
print("finished!")