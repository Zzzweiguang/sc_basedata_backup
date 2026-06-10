#!/usr/bin/env python
"""Run full-atlas marker-high signature score UMAPs without loading all h5ad layers."""

from __future__ import annotations

import gc
import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
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
    """Reuse HUMAN_MARKER_SETS and HUMAN_VALIDATION_PLAN from the notebook."""
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


def filter_marker_dict(marker_dict, var_names, min_genes=2):
    """Keep markers present in the AnnData object and save a coverage report."""
    var_names = pd.Index(var_names.astype(str) if hasattr(var_names, "astype") else var_names)
    filtered = {}
    rows = []
    for state, genes in marker_dict.items():
        genes = list(dict.fromkeys(genes))
        present = [gene for gene in genes if gene in var_names]
        missing = [gene for gene in genes if gene not in var_names]
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


def score_column_name(parent, state):
    """Build the same score-column style used in the notebook."""
    return f"human_full_atlas_{parent.replace(' ', '_')}_{safe_state_name(state)}"


def marker_mean_score(adata, genes):
    """Calculate a fast marker mean-expression score for full-atlas plotting."""
    values = adata[:, genes].X.mean(axis=1)
    return np.asarray(values).ravel()


def load_minimal_anndata():
    """Load only X, var names, and UMAP; skip extra h5ad layers to reduce memory."""
    print(f"读取 h5ad 主矩阵: {HUMAN_H5AD}", flush=True)
    with h5py.File(HUMAN_H5AD, "r") as handle:
        shape = tuple(handle["X"].attrs["shape"])
        print(f"X shape: {shape}", flush=True)

        data = handle["X/data"][:]
        indices = handle["X/indices"][:]
        indptr = handle["X/indptr"][:]
        x = sparse.csr_matrix((data, indices, indptr), shape=shape)

        var_names = pd.Index(decode_array(handle["var/original_gene_names"][:]))
        var_names = make_index_unique(var_names)
        umap = handle["obsm/X_umap"][:]

    adata = AnnData(X=x, var=pd.DataFrame(index=var_names))
    adata.obsm["X_umap"] = umap
    return adata


def plot_parent_scores(adata, parent, marker_dict, ncols=3):
    """Plot all candidate-state scores for one parent class on the full UMAP."""
    cols = list(marker_dict)
    n_panels = len(cols)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.0 * ncols, 3.6 * nrows),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()
    umap = adata.obsm["X_umap"]

    for ax, state in zip(axes, cols):
        score = np.asarray(adata.obs[score_column_name(parent, state)], dtype=float)
        order = np.argsort(score)
        points = ax.scatter(
            umap[order, 0],
            umap[order, 1],
            c=score[order],
            s=0.12,
            cmap="viridis",
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(state, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.colorbar(points, ax=ax, fraction=0.035, pad=0.01)

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.suptitle(f"human full atlas: {parent} candidate signature scores", fontsize=11)
    prefix = FIG_DIR / f"umap_human_full_atlas_{safe_state_name(parent)}_all_signature_scores"
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"保存: {prefix.with_suffix('.png')}", flush=True)
    print(f"保存: {prefix.with_suffix('.pdf')}", flush=True)


def run():
    marker_sets, validation_plan = load_marker_config()
    adata = load_minimal_anndata()
    completed = []

    for parent in validation_plan:
        print(f"\n## full atlas: {parent}", flush=True)
        marker_dict, report = filter_marker_dict(marker_sets[parent], adata.var_names, min_genes=2)
        report_path = TABLE_DIR / f"human_full_atlas_{parent.replace(' ', '_')}_marker_coverage.tsv"
        report.to_csv(report_path, sep="\t", index=False)
        print(f"marker coverage: {report_path}", flush=True)

        for state, genes in marker_dict.items():
            col = score_column_name(parent, state)
            print(f"computing marker mean score '{col}'", flush=True)
            adata.obs[col] = marker_mean_score(adata, genes)

        plot_parent_scores(adata, parent, marker_dict, ncols=3)

        for state in marker_dict:
            col = score_column_name(parent, state)
            if col in adata.obs:
                del adata.obs[col]
        gc.collect()
        completed.append(parent)

    print("\n已完成这些父类的完整图谱 UMAP:", completed, flush=True)


if __name__ == "__main__":
    run()
