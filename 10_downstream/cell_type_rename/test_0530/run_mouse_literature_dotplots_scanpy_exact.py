#!/usr/bin/env python
"""Run mouse literature-marker dotplots with the same Scanpy path as the notebook."""

from __future__ import annotations

import gc
import json
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


BASE_DIR = Path("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/cell_type_rename/test_0530")
FIG_DIR = BASE_DIR / "figures_mouse_literature_dotplot_scanpy_exact"
TABLE_DIR = BASE_DIR / "tables_mouse_literature_dotplot_scanpy_exact"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

MOUSE_H5AD = Path(
    "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/0521_no_Basophil/output_mouse/"
    "scPoli_concat_level3_marker_allmouse.h5ad"
)
NOTEBOOK = BASE_DIR / "test_0530_mouse_literature_markers.ipynb"
LITERATURE_MARKER_TSV = BASE_DIR / "tables_mouse/mouse_candidate_celltypes_marker_table_simplified.tsv"

RANDOM_STATE = 7
MAX_CELLS_FOR_PLOT = 120_000
MIN_GENES = 2
PARENT_COL = "cell_type_level1_corrected"

LITERATURE_PARENT_TO_ROWS = {
    "T cell": [23, 25, 26, 27],
    "Macrophage": [8, 17, 24, 25, 28, 30, 36, 37, 38],
    "Monocyte": [17, 24, 25, 28, 30, 36],
    "Smooth muscle cell": [1, 3, 4, 5, 6, 7, 10, 12, 13, 14, 15, 16, 19, 20, 22, 29, 31, 32, 33, 34, 35],
    "Endothelial cell": [2, 11, 14, 17, 18],
    "Mast cell": [25],
    "Dendritic cell": [21, 25, 30],
    "B cell": [25],
    "Natural killer cell": [25],
    "Fibroblast": [9, 11, 14, 23],
}

VALIDATION_GROUPBY = {
    "T cell": "cell_type_level2",
    "Macrophage": "cell_type_level3",
    "Monocyte": "cell_type_level2",
    "Smooth muscle cell": "cell_type_level2",
    "Endothelial cell": "cell_type_level3",
    "Mast cell": "cell_type_level2",
    "Dendritic cell": "cell_type_level3",
    "B cell": "cell_type_level2",
    "Natural killer cell": "cell_type_level2",
    "Fibroblast": "cell_type_level2",
}


def load_marker_sets_from_notebook() -> dict[str, dict[str, list[str]]]:
    notebook = json.loads(NOTEBOOK.read_text())
    for cell in notebook["cells"]:
        source = "".join(cell.get("source", []))
        if (
            cell.get("cell_type") == "code"
            and "MOUSE_MARKER_SETS = {" in source
            and "Cd3e-Il7r-high CD4 T cell" in source
        ):
            namespace: dict[str, object] = {}
            exec(source, namespace)
            return namespace["MOUSE_MARKER_SETS"]
    raise RuntimeError(f"MOUSE_MARKER_SETS cell not found in {NOTEBOOK}")


def parse_marker_list(value) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def state_name_from_literature_row(row: pd.Series) -> str:
    return f"{row['主要来源GSE']} | {row['候选细胞类型/状态']}"


def build_marker_sets_from_literature(table: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    row_index = table.set_index("序号", drop=False)
    marker_sets: dict[str, dict[str, list[str]]] = {}
    for parent, row_ids in LITERATURE_PARENT_TO_ROWS.items():
        states: dict[str, list[str]] = {}
        for row_id in row_ids:
            row = row_index.loc[row_id]
            genes = parse_marker_list(row["marker基因"])
            if genes:
                states[state_name_from_literature_row(row)] = list(dict.fromkeys(genes))
        marker_sets[parent] = states
    return marker_sets


def normalize_values(values):
    if values is None:
        return None
    if isinstance(values, str):
        return [values]
    return list(values)


def subset_parent(adata, parent_values, parent_col=PARENT_COL, copy=True):
    parent_values = normalize_values(parent_values)
    if parent_values is None:
        return adata.copy() if copy else adata
    if parent_col not in adata.obs.columns:
        raise KeyError(f"{parent_col!r} is not in adata.obs")
    mask = adata.obs[parent_col].astype(str).isin(parent_values)
    out = adata[mask]
    return out.copy() if copy else out


def filter_marker_dict(marker_dict, var_names, min_genes=MIN_GENES, verbose=True):
    var_names = pd.Index(var_names.astype(str) if hasattr(var_names, "astype") else var_names)
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
            print(f"{state}: {len(present)}/{len(genes)} markers present", flush=True)
            if missing:
                print("  missing:", ", ".join(missing), flush=True)
    return out, pd.DataFrame(report)


def downsample_for_plot(adata, max_cells=MAX_CELLS_FOR_PLOT, random_state=RANDOM_STATE):
    if adata.n_obs <= max_cells:
        return adata
    return sc.pp.subsample(adata, n_obs=max_cells, random_state=random_state, copy=True)


def split_marker_dict(marker_dict, max_genes=60):
    chunks = []
    current = {}
    n_genes = 0
    for state, genes in marker_dict.items():
        if current and n_genes + len(genes) > max_genes:
            chunks.append(current)
            current = {}
            n_genes = 0
        current[state] = genes
        n_genes += len(genes)
    if current:
        chunks.append(current)
    return chunks


def safe_file_name(value: str) -> str:
    return value.replace(" ", "_")


def plot_dotplot(adata, marker_dict, groupby, title=None, save_path=None):
    if groupby not in adata.obs.columns:
        raise KeyError(f"{groupby!r} is not in adata.obs")
    chunks = split_marker_dict(marker_dict, max_genes=60)
    saved = []
    for idx, chunk in enumerate(chunks, start=1):
        n_genes = sum(len(values) for values in chunk.values())
        figsize = (max(6, min(24, n_genes * 0.30)), max(4, 0.30 * adata.obs[groupby].nunique()))
        if title:
            print(f"{title} - part {idx}/{len(chunks)}", flush=True)
        dotplot = sc.pl.dotplot(
            adata,
            groupby=groupby,
            var_names=chunk,
            standard_scale="var",
            swap_axes=False,
            figsize=figsize,
            return_fig=True,
        )
        dotplot.legend(width=2.2)
        if save_path:
            save_path = Path(save_path)
            suffix = f"_part{idx}" if len(chunks) > 1 else ""
            out = save_path.with_name(save_path.stem + suffix + save_path.suffix)
            dotplot.savefig(str(out))
            saved.append(out)
            print(f"saved {out}", flush=True)
        plt.close("all")
    return saved


def clear_unused_layers(adata) -> None:
    for layer_name in list(adata.layers.keys()):
        del adata.layers[layer_name]
    gc.collect()


def main() -> None:
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white", fontsize=8)

    marker_sets = load_marker_sets_from_notebook()

    print(f"reading {MOUSE_H5AD}", flush=True)
    adata_mouse = sc.read_h5ad(MOUSE_H5AD)
    adata_mouse.var["ensemble_id"] = adata_mouse.var_names
    adata_mouse.var_names = adata_mouse.var["original_gene_names"].astype(str).values
    adata_mouse.var_names_make_unique()
    clear_unused_layers(adata_mouse)
    print(f"adata_mouse shape: {adata_mouse.shape}", flush=True)

    completed = []
    skipped = []
    for parent, groupby in VALIDATION_GROUPBY.items():
        print(f"\n## {parent} / {groupby}", flush=True)
        adata_parent = subset_parent(adata_mouse, parent_values=[parent], parent_col=PARENT_COL, copy=True)
        print(f"parent shape: {adata_parent.shape}", flush=True)
        if groupby in adata_parent.obs.columns:
            print(adata_parent.obs[groupby].astype(str).value_counts().head(50).to_string(), flush=True)

        marker_dict, report = filter_marker_dict(marker_sets[parent], adata_parent.var_names, min_genes=MIN_GENES)
        report.to_csv(TABLE_DIR / f"mouse_{safe_file_name(parent)}_marker_coverage.tsv", sep="\t", index=False)
        if not marker_dict:
            skipped.append({"parent": parent, "reason": f"no marker set passed min_genes={MIN_GENES}"})
            del adata_parent
            gc.collect()
            continue

        adata_plot = downsample_for_plot(adata_parent, max_cells=MAX_CELLS_FOR_PLOT, random_state=RANDOM_STATE).copy()
        print(f"plot shape: {adata_plot.shape}", flush=True)
        saved = plot_dotplot(
            adata_plot,
            marker_dict,
            groupby=groupby,
            title=f"mouse {parent}: literature marker candidates",
            save_path=FIG_DIR / f"dotplot_mouse_{safe_file_name(parent)}_{groupby}.pdf",
        )
        completed.append({"parent": parent, "groupby": groupby, "n_files": len(saved)})
        del adata_parent, adata_plot
        gc.collect()

    pd.DataFrame(completed).to_csv(TABLE_DIR / "mouse_literature_dotplot_completed.tsv", sep="\t", index=False)
    pd.DataFrame(skipped, columns=["parent", "reason"]).to_csv(
        TABLE_DIR / "mouse_literature_dotplot_skipped.tsv",
        sep="\t",
        index=False,
    )
    print("completed:", completed, flush=True)
    print("skipped:", skipped, flush=True)
    print(f"figure dir: {FIG_DIR}", flush=True)
    print(f"table dir: {TABLE_DIR}", flush=True)


if __name__ == "__main__":
    main()
