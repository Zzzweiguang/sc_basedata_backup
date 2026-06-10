import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import jax
import os
import gc

print("Jax GPU:", jax.devices())

embedding_obsm_keys_list = [
    "PCA",
    "scVI",
    "scANVI",
    "scPoli",
    "Harmony",
    "LIGER",
]

input_dir = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/3_small_integration_mouse_0518/output_0518"
save_root = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/3_small_integration/3_small_integration_mouse_0518/output_0518"
os.makedirs(save_root, exist_ok=True)

scpoli_grid = [
    {"embedding_dims": 5,  "eta": 2},
    {"embedding_dims": 5,  "eta": 5},
    {"embedding_dims": 5,  "eta": 10},

    {"embedding_dims": 10, "eta": 2},
    {"embedding_dims": 10, "eta": 5},
    {"embedding_dims": 10, "eta": 10},
    {"embedding_dims": 10, "eta": 15},

    {"embedding_dims": 15, "eta": 2},
    {"embedding_dims": 15, "eta": 5},
    {"embedding_dims": 15, "eta": 10},
    {"embedding_dims": 15, "eta": 15},

    {"embedding_dims": 20, "eta": 5},
    {"embedding_dims": 20, "eta": 10},
]

for params in scpoli_grid:
    dim = params["embedding_dims"]
    eta = params["eta"]

    input_file = os.path.join(
        input_dir,
        f"mouse-small_concat-ensembl-aggr-sparse-scranlog1p-methed-dim{dim}-eta{eta}.h5ad",
    )

    save_dir = os.path.join(
        save_root,
        f"evaluation_dim{dim}_eta{eta}",
    )

    os.makedirs(save_dir, exist_ok=True)

    result_file = os.path.join(save_dir, "results_scIB.csv")

    if os.path.exists(result_file):
        print("=" * 80)
        print(f"Result already exists, skipping dim={dim}, eta={eta}")
        print("Existing result:", result_file)
        print("=" * 80)
        continue

    if not os.path.exists(input_file):
        print("=" * 80)
        print(f"Input file not found, skipping dim={dim}, eta={eta}")
        print("Missing file:", input_file)
        print("=" * 80)
        continue

    print("=" * 80)
    print(f"Evaluating dim={dim}, eta={eta}")
    print("Input:", input_file)
    print("Output:", save_dir)
    print("=" * 80)

    adata = sc.read_h5ad(input_file)
    adata.obs_names_make_unique()

    missing_keys = [
        key for key in embedding_obsm_keys_list
        if key not in adata.obsm
    ]

    if len(missing_keys) > 0:
        print(f"Missing embeddings for dim={dim}, eta={eta}: {missing_keys}")
        print("Skipping this file.")
        del adata
        gc.collect()
        continue

    print("Starting evaluations...")

    bm = Benchmarker(
        adata,
        batch_key="sample",
        label_key="cell_type_level1",
        embedding_obsm_keys=embedding_obsm_keys_list,
        n_jobs=7,
    )

    bm.benchmark()

    print("Finished benchmarker.")
    print("Writing data.")

    df = bm.get_results(min_max_scale=False)

    df.to_csv(f"{save_dir}/results_scIB.csv")
    df.transpose().to_csv(f"{save_dir}/results_transposed_sciB.csv")

    try:
        bm.plot_results_table(save_dir=save_dir)
    except Exception as e:
        print("Plot failed, but CSV files were saved:", e)

    print(f"Finished without errors: dim={dim}, eta={eta}")

    del adata
    del bm
    del df
    gc.collect()

print("Finished all evaluations.")