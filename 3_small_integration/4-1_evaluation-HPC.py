import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import jax
import os

print("Jax GPU:", jax.devices() )

# ## Evaluation

adata = sc.read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output/small-concat-ensembl-aggr-sparse-method-mse.h5ad")

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
bm.plot_results_table(save_dir="/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output/evaluation")
df = bm.get_results(min_max_scale=False)
df.to_csv("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output/evaluation/results_scIB.csv")
df1 = df.transpose()
df1.to_csv("/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/output/evaluation/results_transposed_sciB.csv")
print("Finished without errors")



