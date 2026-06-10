# =========================================
# SEG 计算完整流程（human + mouse）
# 作者：已根据你的数据路径定制
# =========================================

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os

# -------------------------------
# 配置输入路径
# -------------------------------
CONFIG = {
    "human": {
        "h5ad": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/output_0121/Big-atlas-withlevel3-human_withUMAP2.h5ad",
        "hkg": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/Housekeeping_GenesHuman.csv",
        "gene_length": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/FPKM/human_gene_length.tsv"
    },
    "mouse": {
        "h5ad": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/output_0121/Big-atlas-withlevel3-mouse_withUMAP2.h5ad",
        "hkg": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/Housekeeping_GenesMouse.csv",
        "gene_length": "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/SEGs/FPKM/mouse_gene_length.tsv"
    }
}

# -------------------------------
# 主函数
# -------------------------------
def run_seg_pipeline(species, cfg):

    print(f"\n========== 处理 {species} ==========")

    # -------------------------------
    # Step 1: 读取数据
    # -------------------------------
    adata = sc.read_h5ad(cfg["h5ad"])

    print("细胞数:", adata.n_obs)
    print("基因数:", adata.n_vars)

    # -------------------------------
    # Step 2: zero rate
    # -------------------------------
    print("计算 zero rate...")

    if sparse.issparse(adata.X):
        zero_count = adata.n_obs - np.array((adata.X != 0).sum(axis=0)).flatten()
    else:
        zero_count = (adata.X == 0).sum(axis=0)
    zero_rate = np.array(zero_count).flatten() / adata.n_obs

    zero_df = pd.DataFrame({
        "gene": adata.var_names,
        "zero_rate": zero_rate
    })

    zero_df["zero_rank"] = zero_df["zero_rate"].rank(method="average")

    # -------------------------------
    # Step 3: pseudo-bulk（按 tissue）
    # -------------------------------
    print("构建 pseudo-bulk (按 tissue)...")

    group_ids = adata.obs["tissue"].astype("category")
    groups = group_ids.cat.categories

    pseudo_bulk = []

    for g in groups:
        sub = adata[adata.obs["tissue"] == g]
        summed = np.array(sub.X.sum(axis=0)).flatten()
        pseudo_bulk.append(summed)

    pseudo_bulk = np.vstack(pseudo_bulk)

    pseudobulk_df = pd.DataFrame(
        pseudo_bulk,
        index=groups,
        columns=adata.var_names
    )

    # -------------------------------
    # Step 4: FPKM
    # -------------------------------
    print("计算 FPKM...")

    gene_len_df = pd.read_csv(cfg["gene_length"], sep="\t")
    gene_len_dict = dict(zip(gene_len_df["gene"], gene_len_df["length"]))

    lengths = np.array([gene_len_dict.get(g, np.nan) for g in pseudobulk_df.columns])

    valid = ~np.isnan(lengths)
    pseudobulk_df = pseudobulk_df.loc[:, valid]
    lengths = lengths[valid]

    counts = pseudobulk_df.values
    libsize = counts.sum(axis=1, keepdims=True)

    fpkm = counts / (lengths / 1000) / (libsize / 1e6)

    fpkm_df = pd.DataFrame(
        fpkm,
        index=pseudobulk_df.index,
        columns=pseudobulk_df.columns
    )

    # -------------------------------
    # Step 5: CV
    # -------------------------------
    print("计算 CV...")

    mean = fpkm_df.mean(axis=0)
    std = fpkm_df.std(axis=0)

    cv = std / (mean + 1e-8)

    cv_df = pd.DataFrame({
        "gene": fpkm_df.columns,
        "cv": cv
    })

    cv_df["cv_rank"] = cv_df["cv"].rank(method="average")

    # -------------------------------
    # Step 6: merge ranking
    # -------------------------------
    print("合并 ranking...")

    merged = pd.merge(zero_df, cv_df, on="gene")

    merged["final_rank"] = merged["zero_rank"] + merged["cv_rank"]
    merged = merged.sort_values("final_rank")

    # -------------------------------
    # Step 7: SEG数量 = HKG数量
    # -------------------------------
    print("读取 HKG...")
    hkg_df = pd.read_csv(cfg["hkg"], sep=";")
    #去重
    hkg_genes = set(hkg_df.iloc[:, 1])  # 或用列名
    N = len(hkg_genes)

    print("HKG唯一基因数:", N)

    SEG = merged.head(N)

    # -------------------------------
    # 输出
    # -------------------------------
    out_file = f"/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/output/SEGs/SEGs_results/{species}_SEG_result.csv"
    SEG.to_csv(out_file, index=False)

    print(f"输出完成: {out_file}")


# -------------------------------
# 执行
# -------------------------------
for sp, cfg in CONFIG.items():
    run_seg_pipeline(sp, cfg)