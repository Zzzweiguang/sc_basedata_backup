import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import calinski_harabasz_score, silhouette_score
import matplotlib.pyplot as plt


def read_gene_list(path, gene_col=None):
    import pandas as pd
    import os

    ext = os.path.splitext(path)[1].lower()

    if ext == ".csv":
        df = pd.read_csv(path, sep=";")

        if gene_col is not None:
            genes = df[gene_col].astype(str).str.strip().tolist()
        else:
            # 默认优先找这些常见列名
            candidate_cols = ["gene", "Gene", "gene_symbol", "GeneSymbol", "symbol", "SYMBOL"]
            found = [c for c in candidate_cols if c in df.columns]
            if found:
                genes = df[found[0]].astype(str).str.strip().tolist()
            else:
                # 如果没找到，就默认取第一列
                genes = df.iloc[:, 2].astype(str).str.strip().tolist()

    else:
        with open(path, "r", encoding="utf-8") as f:
            genes = [x.strip() for x in f if x.strip()]

    genes = [g for g in genes if g and g.lower() != "nan"]
    return list(dict.fromkeys(genes))


def ensure_log1p_layer(adata, layer_name="log1p"):
    if layer_name in adata.layers:
        return adata

    adata = adata.copy()
    X = adata.X
    max_val = X.max() if not sparse.issparse(X) else X.max()

    if max_val > 20:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers[layer_name] = adata.X.copy()
    else:
        adata.layers[layer_name] = adata.X.copy()

    return adata


def subset_genes(adata, gene_list):
    keep = [g for g in gene_list if g in adata.var_names]
    if len(keep) < 50:
        raise ValueError(f"匹配到的基因太少：{len(keep)}。请检查 gene symbol 是否一致。")
    return adata[:, keep].copy(), keep


def filter_expressed_genes(adata, min_cells=3):
    adata = adata.copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def sample_cells(adata, n_cells=5000, random_state=0, stratify_key=None):
    rng = np.random.default_rng(random_state)

    if adata.n_obs <= n_cells:
        return adata.copy()

    if stratify_key is None or stratify_key not in adata.obs.columns:
        idx = rng.choice(adata.obs_names, size=n_cells, replace=False)
        return adata[idx].copy()

    obs = adata.obs.copy()
    groups = obs[stratify_key].astype(str)
    sampled_idx = []

    for g in groups.unique():
        cells_g = obs.index[groups == g].tolist()
        prop = len(cells_g) / adata.n_obs
        k = max(1, int(round(prop * n_cells)))
        k = min(k, len(cells_g))
        sampled_idx.extend(rng.choice(cells_g, size=k, replace=False).tolist())

    sampled_idx = list(dict.fromkeys(sampled_idx))

    if len(sampled_idx) > n_cells:
        sampled_idx = rng.choice(sampled_idx, size=n_cells, replace=False).tolist()
    elif len(sampled_idx) < n_cells:
        remain = list(set(adata.obs_names) - set(sampled_idx))
        extra = rng.choice(remain, size=n_cells - len(sampled_idx), replace=False).tolist()
        sampled_idx.extend(extra)

    return adata[sampled_idx].copy()


def run_pipeline(
    adata,
    use_layer="log1p",
    n_pcs=30,
    n_neighbors=15,
    resolution=1.0,
    random_state=0,
):
    adata = adata.copy()

    if use_layer is not None:
        X_backup = adata.X.copy()
        adata.X = adata.layers[use_layer].copy()

    sc.pp.scale(adata, max_value=10)
    n_comps = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps, random_state=random_state)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", random_state=random_state)

    if use_layer is not None:
        adata.X = X_backup

    return adata


def compute_metrics(adata, emb_key="X_pca", label_key="leiden"):
    X = adata.obsm[emb_key]
    labels = adata.obs[label_key].astype(str).values
    n_clusters = len(np.unique(labels))

    if n_clusters < 2:
        return {"CH": np.nan, "Silhouette": np.nan, "n_clusters": n_clusters}

    ch = calinski_harabasz_score(X, labels)
    sil = silhouette_score(X, labels)

    return {"CH": ch, "Silhouette": sil, "n_clusters": n_clusters}


def save_umap(adata, outpath, color_keys):
    valid_keys = [k for k in color_keys if k in adata.obs.columns]
    if "leiden" not in valid_keys:
        valid_keys = ["leiden"] + valid_keys

    sc.pl.umap(adata, color=valid_keys, show=False)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def evaluate_one_dataset(
    adata_path,
    hkg_path,
    seg_path,
    species_name,
    outdir,
    artery_col=None,
    celltype_col=None,
    n_repeat=20,
    n_cells=5000,
    n_pcs=30,
    n_neighbors=15,
    resolution=1.0,
):
    os.makedirs(outdir, exist_ok=True)

    print(f"\n===== Processing {species_name} =====")
    adata = sc.read_h5ad(adata_path)
    adata = ensure_log1p_layer(adata, layer_name="log1p")

    hkg = read_gene_list(hkg_path)
    seg = read_gene_list(seg_path)

    records = []
    saved_sets = set()

    for rep in range(n_repeat):
        print(f"[{species_name}] repeat {rep+1}/{n_repeat}")

        ad_sub = sample_cells(
            adata,
            n_cells=n_cells,
            random_state=rep,
            stratify_key=artery_col if artery_col in adata.obs.columns else None
        )

        # all genes
        ad_all = filter_expressed_genes(ad_sub, min_cells=3)
        ad_all = run_pipeline(
            ad_all,
            use_layer="log1p",
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            resolution=resolution,
            random_state=rep,
        )
        m_all = compute_metrics(ad_all)
        records.append({
            "species": species_name,
            "repeat": rep,
            "gene_set": "all_genes",
            "n_genes_used": ad_all.n_vars,
            **m_all
        })
        if "all_genes" not in saved_sets:
            save_umap(
                ad_all,
                os.path.join(outdir, f"{species_name}_UMAP_all_genes.png"),
                color_keys=["leiden", artery_col, celltype_col]
            )
            saved_sets.add("all_genes")

        # HKG
        ad_hkg, keep_hkg = subset_genes(ad_sub, hkg)
        ad_hkg = run_pipeline(
            ad_hkg,
            use_layer="log1p",
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            resolution=resolution,
            random_state=rep,
        )
        m_hkg = compute_metrics(ad_hkg)
        records.append({
            "species": species_name,
            "repeat": rep,
            "gene_set": "HKG",
            "n_genes_used": len(keep_hkg),
            **m_hkg
        })
        if "HKG" not in saved_sets:
            save_umap(
                ad_hkg,
                os.path.join(outdir, f"{species_name}_UMAP_HKG.png"),
                color_keys=["leiden", artery_col, celltype_col]
            )
            saved_sets.add("HKG")

        # SEG
        ad_seg, keep_seg = subset_genes(ad_sub, seg)
        ad_seg = run_pipeline(
            ad_seg,
            use_layer="log1p",
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            resolution=resolution,
            random_state=rep,
        )
        m_seg = compute_metrics(ad_seg)
        records.append({
            "species": species_name,
            "repeat": rep,
            "gene_set": "SEG",
            "n_genes_used": len(keep_seg),
            **m_seg
        })
        if "SEG" not in saved_sets:
            save_umap(
                ad_seg,
                os.path.join(outdir, f"{species_name}_UMAP_SEG.png"),
                color_keys=["leiden", artery_col, celltype_col]
            )
            saved_sets.add("SEG")

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(outdir, f"{species_name}_metrics_all_repeats.csv"), index=False)

    summary = (
        df.groupby(["species", "gene_set"], as_index=False)
          .agg(
              CH_mean=("CH", "mean"),
              CH_std=("CH", "std"),
              Silhouette_mean=("Silhouette", "mean"),
              Silhouette_std=("Silhouette", "std"),
              n_clusters_mean=("n_clusters", "mean"),
              n_genes_used_mean=("n_genes_used", "mean"),
          )
    )
    summary.to_csv(os.path.join(outdir, f"{species_name}_metrics_summary.csv"), index=False)

    return df, summary


def plot_summary(summary, species_name, outdir):
    order = ["all_genes", "HKG", "SEG"]
    df = summary.copy()
    df["gene_set"] = pd.Categorical(df["gene_set"], categories=order, ordered=True)
    df = df.sort_values("gene_set")

    plt.figure(figsize=(5, 4))
    plt.bar(df["gene_set"], df["CH_mean"], yerr=df["CH_std"], capsize=4)
    plt.ylabel("Calinski-Harabasz index")
    plt.title(species_name)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{species_name}_CH_barplot.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(5, 4))
    plt.bar(df["gene_set"], df["Silhouette_mean"], yerr=df["Silhouette_std"], capsize=4)
    plt.ylabel("Silhouette score")
    plt.title(species_name)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{species_name}_Silhouette_barplot.png"), dpi=300)
    plt.close()


if __name__ == "__main__":
    # ===== 这里改成你的文件名 =====
    HUMAN_H5AD = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/output_0121/Big-atlas-withlevel3-human_withUMAP2.h5ad"
    MOUSE_H5AD = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/output_0121/Big-atlas-withlevel3-mouse_withUMAP2.h5ad"

    HKG_HUMAN = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/Housekeeping_GenesHuman.csv"
    HKG_MOUSE = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/Housekeeping_GenesMouse.csv"
    SEG_HUMAN = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/output/SEGs/SEGs_results/human_SEG_result.csv"
    SEG_MOUSE = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/output/SEGs/SEGs_results/mouse_SEG_result.csv"

    # ===== 这里改成你 obs 里的列名 =====
    ARTERY_COL_HUMAN = "tissue"   # 没有就改成 None
    CELLTYPE_COL_HUMAN = "cell_type_level2"   # 没有就改成 None

    ARTERY_COL_MOUSE = "tissue"
    CELLTYPE_COL_MOUSE = "cell_type_level2"

    OUTDIR = "results_seg_eval"

    # human
    df_h, sum_h = evaluate_one_dataset(
        adata_path=HUMAN_H5AD,
        hkg_path=HKG_HUMAN,
        seg_path=SEG_HUMAN,
        species_name="human",
        outdir=OUTDIR,
        artery_col=ARTERY_COL_HUMAN,
        celltype_col=CELLTYPE_COL_HUMAN,
        n_repeat=20,
        n_cells=5000,
        n_pcs=30,
        n_neighbors=15,
        resolution=1.0,
    )
    plot_summary(sum_h, "human", OUTDIR)

    # mouse
    df_m, sum_m = evaluate_one_dataset(
        adata_path=MOUSE_H5AD,
        hkg_path=HKG_MOUSE,
        seg_path=SEG_MOUSE,
        species_name="mouse",
        outdir=OUTDIR,
        artery_col=ARTERY_COL_MOUSE,
        celltype_col=CELLTYPE_COL_MOUSE,
        n_repeat=20,
        n_cells=5000,
        n_pcs=30,
        n_neighbors=15,
        resolution=1.0,
    )
    plot_summary(sum_m, "mouse", OUTDIR)

    print("\nAll done.")