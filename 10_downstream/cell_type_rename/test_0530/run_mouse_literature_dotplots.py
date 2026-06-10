#!/usr/bin/env python
"""Run mouse literature-marker dotplots without modifying existing notebooks."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata.utils import make_index_unique
from scipy import sparse


BASE_DIR = Path("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/cell_type_rename/test_0530")
H5AD = Path(
    "/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/0521_no_Basophil/"
    "output_mouse/scPoli_concat_level3_marker_allmouse.h5ad"
)
LITERATURE_MARKER_TSV = BASE_DIR / "tables_mouse/mouse_candidate_celltypes_marker_table_simplified.tsv"
FIG_DIR = BASE_DIR / "figures_mouse_literature_dotplot"
TABLE_DIR = BASE_DIR / "tables_mouse_literature_dotplot"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

PARENT_COL = "cell_type_level1_corrected"
MAX_GENES_PER_DOTPLOT = 60
ROW_CHUNK_SIZE = 5000

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

LOW_SPECIFICITY_PARENT_NOTES = {
    "Mast cell": "from GSE245373 CD45+ leukocyte meta-clusters; not a mast-specific marker panel.",
    "B cell": "from GSE245373 CD45+ leukocyte meta-clusters; B-cell-specific evidence mainly includes Cd79a.",
    "Natural killer cell": "from GSE245373 CD45+ leukocyte meta-clusters; NK/NKT evidence mainly includes Nkg7/Cd3e/Cd8a.",
}


def decode_array(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_categorical_obs(handle: h5py.File, column: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    group = handle["obs"][column]
    categories = decode_array(group["categories"][:])
    codes = group["codes"][:].astype(np.int64)
    values = np.array(["NA" if code < 0 else categories[code] for code in codes], dtype=object)
    return values, categories, codes


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


def split_marker_dict(marker_dict: dict[str, list[str]], max_genes: int = MAX_GENES_PER_DOTPLOT) -> list[dict[str, list[str]]]:
    chunks: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] = {}
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


def safe_name(value: str) -> str:
    out = value.replace("/", "_").replace(" ", "_").replace("+", "pos").replace("-", "_")
    for char in ",;()[]|:":
        out = out.replace(char, "")
    return out


def filter_marker_dict(
    marker_dict: dict[str, list[str]],
    var_names: pd.Index,
    min_genes: int = 1,
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    filtered: dict[str, list[str]] = {}
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


def flatten_marker_dict(marker_dict: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for values in marker_dict.values():
        genes.extend(values)
    return list(dict.fromkeys(genes))


def standard_scale_var(df: pd.DataFrame) -> pd.DataFrame:
    scaled = df.copy()
    scaled = scaled - scaled.min(axis=0)
    denom = scaled.max(axis=0).replace(0, np.nan)
    return scaled.divide(denom, axis=1).fillna(0)


def plot_dotplot_from_summary(
    mean_df: pd.DataFrame,
    fraction_df: pd.DataFrame,
    marker_dict: dict[str, list[str]],
    groupby: str,
    parent: str,
) -> list[Path]:
    saved: list[Path] = []
    chunks = split_marker_dict(marker_dict)
    for part, chunk in enumerate(chunks, start=1):
        genes = flatten_marker_dict(chunk)
        color_df = standard_scale_var(mean_df.loc[:, genes])
        size_df = fraction_df.loc[:, genes]

        obs = pd.DataFrame({groupby: pd.Categorical(color_df.index, categories=list(color_df.index), ordered=True)})
        dummy = AnnData(
            X=np.zeros((len(color_df.index), len(genes)), dtype=np.float32),
            obs=obs,
            var=pd.DataFrame(index=pd.Index(genes)),
        )

        n_genes = len(genes)
        n_groups = len(color_df.index)
        figsize = (max(6, min(24, n_genes * 0.30)), max(4, 0.30 * n_groups))
        title = f"mouse {parent}: literature marker candidates"
        if len(chunks) > 1:
            title = f"{title} - part {part}/{len(chunks)}"

        dotplot = sc.pl.DotPlot(
            dummy,
            var_names=chunk,
            groupby=groupby,
            dot_color_df=color_df,
            dot_size_df=size_df,
            figsize=figsize,
            title=title,
        )
        dotplot.legend(
            width=2.2,
            colorbar_title="Mean expression\nstandardized by gene",
            size_title="Fraction of cells\nin group (%)",
        )
        suffix = f"_part{part}" if len(chunks) > 1 else ""
        prefix = FIG_DIR / f"dotplot_mouse_{safe_name(parent)}_{groupby}_literature_markers{suffix}"
        pdf_path = prefix.with_suffix(".pdf")
        png_path = prefix.with_suffix(".png")
        dotplot.savefig(str(pdf_path))
        dotplot.savefig(str(png_path))
        plt.close("all")
        saved.extend([pdf_path, png_path])
        print(f"saved {pdf_path}", flush=True)
        print(f"saved {png_path}", flush=True)
    return saved


def init_summary(parent: str, groupby: str, group_categories: list[str], group_values: np.ndarray, parent_values: np.ndarray, n_markers: int):
    in_parent = parent_values == parent
    observed = [group for group in group_categories if np.any(in_parent & (group_values == group))]
    if not observed:
        observed = sorted(pd.unique(group_values[in_parent]).tolist())
    group_to_idx = {group: idx for idx, group in enumerate(observed)}
    group_codes = np.full(group_values.shape[0], -1, dtype=np.int64)
    for group, idx in group_to_idx.items():
        group_codes[group_values == group] = idx
    counts = np.bincount(group_codes[in_parent & (group_codes >= 0)], minlength=len(observed)).astype(np.int64)
    return {
        "parent": parent,
        "groupby": groupby,
        "groups": observed,
        "group_to_idx": group_to_idx,
        "group_codes": group_codes,
        "parent_mask": in_parent,
        "counts": counts,
        "sum": np.zeros((len(observed), n_markers), dtype=np.float64),
        "nz": np.zeros((len(observed), n_markers), dtype=np.float64),
    }


def update_summary_from_chunk(summary: dict, chunk_marker: sparse.csr_matrix, start: int, end: int) -> None:
    parent_mask = summary["parent_mask"][start:end]
    group_codes = summary["group_codes"][start:end]
    selected = parent_mask & (group_codes >= 0)
    if not np.any(selected):
        return
    sub = chunk_marker[selected]
    codes = group_codes[selected]
    n_selected = sub.shape[0]
    indicator = sparse.csr_matrix(
        (np.ones(n_selected, dtype=np.float64), (codes, np.arange(n_selected))),
        shape=(len(summary["groups"]), n_selected),
    )
    summary["sum"] += (indicator @ sub).toarray()
    summary["nz"] += (indicator @ (sub > 0)).toarray()


def scan_marker_summaries(
    handle: h5py.File,
    marker_indices: list[int],
    summaries: dict[str, dict],
    n_obs: int,
    n_vars: int,
) -> None:
    x_group = handle["X"]
    indptr = x_group["indptr"][:]
    data_ds = x_group["data"]
    indices_ds = x_group["indices"]
    marker_indices_array = np.array(marker_indices, dtype=np.int64)

    for start in range(0, n_obs, ROW_CHUNK_SIZE):
        end = min(start + ROW_CHUNK_SIZE, n_obs)
        data_start = int(indptr[start])
        data_end = int(indptr[end])
        row_data = data_ds[data_start:data_end]
        row_indices = indices_ds[data_start:data_end]
        row_indptr = indptr[start : end + 1] - data_start
        chunk_full = sparse.csr_matrix((row_data, row_indices, row_indptr), shape=(end - start, n_vars))
        chunk_marker = chunk_full[:, marker_indices_array].astype(np.float64)
        for summary in summaries.values():
            update_summary_from_chunk(summary, chunk_marker, start, end)
        if start == 0 or end == n_obs or (start // ROW_CHUNK_SIZE) % 10 == 0:
            print(f"processed rows {end}/{n_obs}", flush=True)


def save_expression_tables(
    parent: str,
    summary: dict,
    marker_dict: dict[str, list[str]],
    all_marker_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = summary["counts"].astype(float)
    denom = np.where(counts[:, None] > 0, counts[:, None], np.nan)
    mean = pd.DataFrame(summary["sum"] / denom, index=summary["groups"], columns=all_marker_genes).fillna(0)
    fraction = pd.DataFrame(summary["nz"] / denom, index=summary["groups"], columns=all_marker_genes).fillna(0)

    genes = flatten_marker_dict(marker_dict)
    mean.loc[:, genes].to_csv(
        TABLE_DIR / f"mouse_{safe_name(parent)}_mean_expression_by_{summary['groupby']}.tsv",
        sep="\t",
    )
    fraction.loc[:, genes].to_csv(
        TABLE_DIR / f"mouse_{safe_name(parent)}_fraction_expressing_by_{summary['groupby']}.tsv",
        sep="\t",
    )
    pd.DataFrame({"group": summary["groups"], "n_cells": summary["counts"]}).to_csv(
        TABLE_DIR / f"mouse_{safe_name(parent)}_group_counts_by_{summary['groupby']}.tsv",
        sep="\t",
        index=False,
    )
    return mean, fraction


def main() -> None:
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, facecolor="white", fontsize=8)

    literature_table = pd.read_csv(LITERATURE_MARKER_TSV, sep="\t")
    marker_sets = build_marker_sets_from_literature(literature_table)

    for parent, note in LOW_SPECIFICITY_PARENT_NOTES.items():
        print(f"low specificity note - {parent}: {note}", flush=True)

    with h5py.File(H5AD, "r") as handle:
        n_obs, n_vars = tuple(handle["X"].attrs["shape"])
        var_names = pd.Index(make_index_unique(pd.Index(decode_array(handle["var/original_gene_names"][:]))))
        parent_values, _, _ = read_categorical_obs(handle, PARENT_COL)

        group_values_by_col: dict[str, np.ndarray] = {}
        group_categories_by_col: dict[str, list[str]] = {}
        for groupby in sorted(set(VALIDATION_GROUPBY.values())):
            values, categories, _ = read_categorical_obs(handle, groupby)
            group_values_by_col[groupby] = values
            group_categories_by_col[groupby] = categories

        filtered_marker_sets: dict[str, dict[str, list[str]]] = {}
        coverage_rows = []
        for parent, marker_dict in marker_sets.items():
            filtered, coverage = filter_marker_dict(marker_dict, var_names, min_genes=1)
            coverage.insert(0, "parent", parent)
            coverage.to_csv(TABLE_DIR / f"mouse_{safe_name(parent)}_marker_coverage.tsv", sep="\t", index=False)
            coverage_rows.append(coverage)
            filtered_marker_sets[parent] = filtered
            print(f"{parent}: {len(filtered)}/{len(marker_dict)} states pass marker coverage", flush=True)

        pd.concat(coverage_rows, ignore_index=True).to_csv(
            TABLE_DIR / "mouse_literature_marker_coverage_all_parents.tsv",
            sep="\t",
            index=False,
        )

        all_marker_genes = sorted({gene for states in filtered_marker_sets.values() for genes in states.values() for gene in genes})
        marker_to_var_idx = {gene: int(var_names.get_loc(gene)) for gene in all_marker_genes}
        marker_indices = [marker_to_var_idx[gene] for gene in all_marker_genes]
        print(f"using {len(all_marker_genes)} unique literature markers present in adata", flush=True)

        summaries = {
            parent: init_summary(
                parent=parent,
                groupby=VALIDATION_GROUPBY[parent],
                group_categories=group_categories_by_col[VALIDATION_GROUPBY[parent]],
                group_values=group_values_by_col[VALIDATION_GROUPBY[parent]],
                parent_values=parent_values,
                n_markers=len(all_marker_genes),
            )
            for parent in VALIDATION_GROUPBY
        }

        scan_marker_summaries(handle, marker_indices, summaries, n_obs=n_obs, n_vars=n_vars)

    completed = []
    skipped = []
    for parent, summary in summaries.items():
        marker_dict = filtered_marker_sets[parent]
        if not marker_dict:
            skipped.append({"parent": parent, "reason": "no marker set passed coverage filtering"})
            continue
        if np.sum(summary["counts"]) == 0:
            skipped.append({"parent": parent, "reason": f"no cells in {PARENT_COL}"})
            continue
        mean, fraction = save_expression_tables(parent, summary, marker_dict, all_marker_genes)
        plot_dotplot_from_summary(mean, fraction, marker_dict, summary["groupby"], parent)
        completed.append(parent)

    pd.DataFrame({"parent": completed}).to_csv(TABLE_DIR / "mouse_literature_dotplot_completed.tsv", sep="\t", index=False)
    pd.DataFrame(skipped, columns=["parent", "reason"]).to_csv(
        TABLE_DIR / "mouse_literature_dotplot_skipped.tsv",
        sep="\t",
        index=False,
    )
    print("completed parents:", completed, flush=True)
    print("skipped:", skipped, flush=True)
    print(f"figure dir: {FIG_DIR}", flush=True)
    print(f"table dir: {TABLE_DIR}", flush=True)


if __name__ == "__main__":
    main()
