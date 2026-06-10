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


adata_latent_full = sc.read_h5ad("./output_260420/Altas_level1_all_human_preumap_corrected_3.h5ad")
print("umap!!!")
# adata_latent_full = adata_latent_full[adata_latent_full.obs["dataset"] != "IAISR"].copy()
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
# adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_umap_nogene_noIAISR_corrected.h5ad")
adata_latent_full.write_h5ad("./output_260420/Altas_level1_all_human_umap_nogene_corrected_3.h5ad")##有IAISR

print("finished!!!")