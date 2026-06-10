import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import jax
import os

print("Jax GPU:", jax.devices())

adata = sc.read_h5ad("./output_260420/small-concat-ensembl-aggr-sparse-method-noCxG-mse-dim10-eta4-noIAISR.h5ad")

adata.obs_names_make_unique()

print("Starting evaluations...")

embedding_obsm_keys_list = ['PCA', 'scVI', 'scANVI', 'scPoli', 'Harmony', 'LIGER']

bm = Benchmarker(
    adata,
    batch_key="sample",
    label_key="cell_type_level1",
    embedding_obsm_keys=embedding_obsm_keys_list,
    n_jobs=7,
)

bm.benchmark()

print("Finished benchmarker..")
print("Writing data..")

save_dir = "./output_260420/evaluation_dim10_eta4_noIAISR"
os.makedirs(save_dir, exist_ok=True)

df = bm.get_results(min_max_scale=False)
df.to_csv(f"{save_dir}/results_scIB.csv")
df.transpose().to_csv(f"{save_dir}/results_transposed_sciB.csv")

try:
    bm.plot_results_table(save_dir=save_dir)
except Exception as e:
    print("Plot failed, but CSV files were saved:", e)

print("Finished without errors")