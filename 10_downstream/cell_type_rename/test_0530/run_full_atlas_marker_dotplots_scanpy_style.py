#!/usr/bin/env python
"""Create Scanpy-style full-atlas dotplots for candidate marker genes."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
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


def read_categorical_obs(handle, column):
    """Read an AnnData categorical obs column from h5ad."""
    group = handle["obs"][column]
    categories = pd.Index(decode_array(group["categories"][:]))
    codes = group["codes"][:]
    values = ["NA" if code < 0 else categories[code] for code in codes]
    return pd.Categorical(values, categories=list(categories), ordered=False)


def filter_marker_dict(marker_dict, var_names, min_genes=2, verbose=True):
    """Keep genes present in the marker-only AnnData object."""
    var_names = pd.Index(var_names)
    out = {}
    report = []
    for state, genes in marker_dict.items():
        genes = list(dict.fromkeys(genes))
        present = [gene for gene in genes if gene in var_names]
        missing = [gene for gene in genes if gene not in var_names]
        report.append(
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
            out[state] = present
        if verbose:
            print(f"{state}: {len(present)}/{len(genes)} 个 marker 存在", flush=True)
            if missing:
                print("  缺失:", ", ".join(missing), flush=True)
    return out, pd.DataFrame(report)


def load_marker_only_anndata(marker_sets, groupby):
    """Load only candidate marker genes plus the grouping column."""
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
        group_values = read_categorical_obs(handle, groupby)

        var_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
        present_markers = [gene for gene in all_markers if gene in var_to_idx]
        marker_indices = [var_to_idx[gene] for gene in present_markers]
        print(f"用于 dotplot 的 marker: {len(present_markers)}/{len(all_markers)}", flush=True)

        data = handle["X/data"][:]
        indices = handle["X/indices"][:]
        indptr = handle["X/indptr"][:]
        x = sparse.csr_matrix((data, indices, indptr), shape=shape)
        x_marker = x[:, marker_indices].copy()

    adata = AnnData(
        X=x_marker,
        obs=pd.DataFrame({groupby: group_values}),
        var=pd.DataFrame(index=pd.Index(present_markers)),
    )
    return adata


def plot_scanpy_dotplot(adata, marker_dict, groupby, title, save_prefix):
    """Call Scanpy dotplot with the same options used in test_0530_2.ipynb."""
    n_genes = sum(len(genes) for genes in marker_dict.values())
    n_groups = adata.obs[groupby].nunique()
    figsize = (max(6, min(24, n_genes * 0.30)), max(4, 0.30 * n_groups))
    print(title, flush=True)
    dp = sc.pl.dotplot(
        adata,
        groupby=groupby,
        var_names=marker_dict,
        standard_scale="var",
        swap_axes=False,
        figsize=figsize,
        return_fig=True,
    )
    dp.legend(width=2.2)
    dp.savefig(str(save_prefix.with_suffix(".pdf")))
    dp.savefig(str(save_prefix.with_suffix(".png")))
    dp.show()
    plt.close("all")


def run(groupby="cell_type_level2", min_genes=2):
    """Run all seven parent-class full-atlas Scanpy-style dotplots."""
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white", fontsize=8)

    marker_sets, validation_plan = load_marker_config()
    adata = load_marker_only_anndata(marker_sets, groupby)

    completed = []
    for parent in validation_plan:
        print(f"\n## full atlas Scanpy dotplot: {parent}", flush=True)
        marker_dict, report = filter_marker_dict(marker_sets[parent], adata.var_names, min_genes=min_genes)
        report_path = TABLE_DIR / (
            f"human_full_atlas_{safe_state_name(parent)}"
            f"_all_candidate_scanpy_dotplot_marker_coverage_by_{groupby}.tsv"
        )
        report.to_csv(report_path, sep="\t", index=False)

        save_prefix = FIG_DIR / (
            f"dotplot_human_full_atlas_{safe_state_name(parent)}"
            f"_all_candidate_markers_by_{groupby}_scanpy_style"
        )
        plot_scanpy_dotplot(
            adata,
            marker_dict,
            groupby=groupby,
            title=f"human full atlas: {parent} candidate markers by {groupby}",
            save_prefix=save_prefix,
        )
        print(f"保存: {save_prefix.with_suffix('.pdf')}", flush=True)
        print(f"保存: {save_prefix.with_suffix('.png')}", flush=True)
        completed.append(parent)

    print("\n已完成这些父类的 Scanpy 格式完整图谱 DotPlot:", completed, flush=True)


if __name__ == "__main__":
    run()
