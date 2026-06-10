#!/usr/bin/env python
"""Create full-atlas dotplots for candidate marker genes grouped by all cell types."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata.utils import make_index_unique
from scipy import sparse


BASE_DIR = Path("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/cell_type_rename/test_0530")
NOTEBOOK = BASE_DIR / "test_0530_2.ipynb"
HUMAN_H5AD = Path(
    "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/0521_no_Basophil/"
    "output_allhuman/scPoli_concat_level3_marker_all_metadata.h5ad"
)
FIG_DIR = BASE_DIR / "figures_human"
TABLE_DIR = BASE_DIR / "tables_human"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def load_marker_config():
    """Reuse HUMAN_MARKER_SETS and HUMAN_VALIDATION_PLAN from test_0530_2.ipynb."""
    nb = json.loads(NOTEBOOK.read_text())
    ns = {"pd": pd}
    exec("".join(nb["cells"][9]["source"]), ns)
    exec("".join(nb["cells"][13]["source"]), ns)
    return ns["HUMAN_MARKER_SETS"], ns["HUMAN_VALIDATION_PLAN"]


def decode_array(values):
    """Decode h5py byte/object arrays into plain strings."""
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8"))
        else:
            out.append(str(value))
    return out


def safe_state_name(state):
    """Convert a state name to a stable filename fragment."""
    return (
        state.replace("/", "_")
        .replace(" ", "_")
        .replace("+", "pos")
        .replace("-", "_")
        .replace(",", "")
        .replace(";", "")
        .replace("(", "")
        .replace(")", "")
    )


def read_categorical_obs(handle, column):
    """Read an AnnData categorical obs column from h5ad."""
    group = handle["obs"][column]
    categories = np.asarray(decode_array(group["categories"][:]), dtype=object)
    codes = group["codes"][:]
    values = np.asarray(["NA" if code < 0 else categories[code] for code in codes], dtype=object)
    return values, categories


def marker_coverage(marker_dict, var_names, min_genes=2):
    """Filter marker sets by availability in the h5ad var_names."""
    var_index = pd.Index(var_names)
    filtered = {}
    rows = []
    for state, genes in marker_dict.items():
        genes = list(dict.fromkeys(genes))
        present = [gene for gene in genes if gene in var_index]
        missing = [gene for gene in genes if gene not in var_index]
        rows.append(
            {
                "state": state,
                "n_total": len(genes),
                "n_present": len(present),
                "pct_present": len(present) / len(genes) if genes else 0,
                "present": ",".join(present),
                "missing": ",".join(missing),
            }
        )
        if len(present) >= min_genes:
            filtered[state] = present
    return filtered, pd.DataFrame(rows)


def read_h5ad_for_dotplot(marker_sets, groupby):
    """Load the full sparse X, gene names, and one obs grouping column."""
    all_markers = sorted(
        {
            gene
            for parent_states in marker_sets.values()
            for genes in parent_states.values()
            for gene in genes
        }
    )
    print(f"读取 h5ad 主矩阵: {HUMAN_H5AD}", flush=True)
    with h5py.File(HUMAN_H5AD, "r") as handle:
        shape = tuple(handle["X"].attrs["shape"])
        print(f"X shape: {shape}", flush=True)
        var_names = pd.Index(decode_array(handle["var/original_gene_names"][:]))
        var_names = pd.Index(make_index_unique(var_names))
        group_values, group_categories = read_categorical_obs(handle, groupby)

        var_to_idx = {gene: i for i, gene in enumerate(var_names)}
        present_markers = [gene for gene in all_markers if gene in var_to_idx]
        marker_indices = np.array([var_to_idx[gene] for gene in present_markers], dtype=np.int32)
        print(f"用于 dotplot 的 marker: {len(present_markers)}/{len(all_markers)}", flush=True)

        data = handle["X/data"][:]
        indices = handle["X/indices"][:]
        indptr = handle["X/indptr"][:]
        x = sparse.csr_matrix((data, indices, indptr), shape=shape)
        x_marker = x[:, marker_indices].tocsc()

    return x_marker, pd.Index(present_markers), group_values, pd.Index(group_categories)


def compute_group_stats(x_marker, marker_names, group_values, group_order):
    """Compute mean expression and fraction positive for each marker by group."""
    group_values = pd.Series(group_values, dtype="object")
    rows = []
    for group in group_order:
        mask = group_values.values == group
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        sub = x_marker[mask, :]
        mean_expr = np.asarray(sub.mean(axis=0)).ravel()
        pct_expr = np.asarray((sub > 0).mean(axis=0)).ravel()
        for gene, mean_value, pct_value in zip(marker_names, mean_expr, pct_expr):
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "mean_expr": float(mean_value),
                    "pct_expr": float(pct_value),
                    "n_cells": n_cells,
                }
            )
    return pd.DataFrame(rows)


def build_parent_dotplot_data(stats, marker_dict, parent):
    """Prepare one parent-class dotplot table."""
    gene_to_state = {}
    ordered_genes = []
    for state, genes in marker_dict.items():
        for gene in genes:
            if gene not in gene_to_state:
                gene_to_state[gene] = state
                ordered_genes.append(gene)
    df = stats[stats["gene"].isin(ordered_genes)].copy()
    df["state"] = df["gene"].map(gene_to_state)
    df["gene"] = pd.Categorical(df["gene"], categories=ordered_genes, ordered=True)
    df["parent"] = parent
    max_mean = df["mean_expr"].max()
    df["scaled_mean_expr"] = df["mean_expr"] / max_mean if max_mean > 0 else 0.0
    return df.sort_values(["group", "gene"])


def plot_parent_dotplot(dot_df, parent, group_order, save_prefix):
    """Draw a Scanpy-style dotplot from precomputed group statistics."""
    genes = list(dot_df["gene"].cat.categories)
    x_lookup = {gene: idx for idx, gene in enumerate(genes)}
    y_order = [group for group in group_order if group in set(dot_df["group"])]
    y_lookup = {group: idx for idx, group in enumerate(y_order)}

    fig_w = max(8, min(26, len(genes) * 0.32))
    fig_h = max(5, len(y_order) * 0.34)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    x = dot_df["gene"].astype(str).map(x_lookup)
    y = dot_df["group"].map(y_lookup)
    sizes = 8 + dot_df["pct_expr"].clip(0, 1) * 95
    points = ax.scatter(
        x,
        y,
        s=sizes,
        c=dot_df["scaled_mean_expr"],
        cmap="viridis",
        vmin=0,
        vmax=1,
        linewidths=0.15,
        edgecolors="0.45",
    )

    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order, fontsize=8)
    ax.set_title(f"human full atlas: {parent} candidate markers by cell_type_level2", fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.invert_yaxis()
    ax.grid(axis="x", color="0.92", linewidth=0.5)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    cbar = fig.colorbar(points, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("scaled mean expression", fontsize=8)

    for pct in [0.25, 0.50, 0.75]:
        ax.scatter([], [], s=8 + pct * 95, c="0.7", edgecolors="0.45", label=f"{pct:.0%}")
    ax.legend(title="pct cells", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    fig.savefig(save_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(save_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(groupby="cell_type_level2", min_genes=2):
    """Run all seven parent-class full-atlas marker dotplots."""
    marker_sets, validation_plan = load_marker_config()
    x_marker, marker_names, group_values, group_order = read_h5ad_for_dotplot(marker_sets, groupby)
    stats = compute_group_stats(x_marker, marker_names, group_values, group_order)
    stats_path = TABLE_DIR / f"human_full_atlas_marker_stats_by_{groupby}.tsv"
    stats.to_csv(stats_path, sep="\t", index=False)
    print(f"保存全 marker 分组统计: {stats_path}", flush=True)

    completed = []
    for parent in validation_plan:
        print(f"\n## full atlas marker dotplot: {parent}", flush=True)
        marker_dict, report = marker_coverage(marker_sets[parent], marker_names, min_genes=min_genes)
        report_path = TABLE_DIR / f"human_full_atlas_{safe_state_name(parent)}_all_candidate_dotplot_marker_coverage_by_{groupby}.tsv"
        report.to_csv(report_path, sep="\t", index=False)

        dot_df = build_parent_dotplot_data(stats, marker_dict, parent)
        dot_df_path = TABLE_DIR / f"human_full_atlas_{safe_state_name(parent)}_candidate_marker_dotplot_data_by_{groupby}.tsv"
        dot_df.to_csv(dot_df_path, sep="\t", index=False)

        save_prefix = FIG_DIR / f"dotplot_human_full_atlas_{safe_state_name(parent)}_all_candidate_markers_by_{groupby}"
        plot_parent_dotplot(dot_df, parent, group_order, save_prefix)
        print(f"保存: {save_prefix.with_suffix('.png')}", flush=True)
        print(f"保存: {save_prefix.with_suffix('.pdf')}", flush=True)
        completed.append(parent)

    print("\n已完成这些父类的完整图谱 DotPlot:", completed, flush=True)


if __name__ == "__main__":
    run()
