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

print("load data")
adata_latent_all = sc.read("./output_all_data1226/Atlas-level1-all_data-preumap.h5ad")
print("umap start!!")
sc.pp.neighbors(adata_latent_all, n_neighbors=15)
sc.tl.umap(adata_latent_all)
print("umap finish!!")

print("save umap ed")
adata_latent_all.write_h5ad("./output_all_data1226/Atlas-level1-all_data-nogene.h5ad")
print("finished!!!")