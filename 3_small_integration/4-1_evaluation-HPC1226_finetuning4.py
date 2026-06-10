import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import jax
import os

print("Jax GPU:", jax.devices() )

# ## Evaluation

adata = sc.read_h5ad("./output_1226/small-concat-ensembl-aggr-sparse-method-noCxG-mse-dim10-eta6.h5ad")

print("Starting evaluations...")

embedding_obsm_keys_list = ['PCA', 'scVI', 'scANVI','scPoli','Harmony', 'LIGER']   ##noscgen


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
bm.plot_results_table(save_dir="./output_1226/evaluation_dim10_eta6")
df = bm.get_results(min_max_scale=False)
df.to_csv("./output_1226/evaluation_dim10_eta6/results_scIB.csv")
df1 = df.transpose()
df1.to_csv("./output_1226/evaluation_dim10_eta6/results_transposed_sciB.csv")
print("Finished without errors")



