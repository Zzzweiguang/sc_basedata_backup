#!/usr/bin/env python
"""Plot marker-celltype Spearman correlations in a diagonal paired layout."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


MARKER_GENES = ["ZNF555", "AMFR", "FOXK1", "WDR91", "IL12RB2", "ICAM3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corr-csv",
        default="output_ml_marker_biological_validation/marker_immune_celltype_correlations.csv",
    )
    parser.add_argument("--out-dir", default="output_ml_marker_biological_validation")
    parser.add_argument("--prefix", default="marker_immune_celltype_correlation_diagonal")
    return parser.parse_args()


def setup_plot_style() -> None:
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"
    sns.set_theme(style="ticks")


def save_figure(fig: plt.Figure, out_path: Path) -> None:
    fig.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def load_matrix(corr_csv: Path) -> pd.DataFrame:
    corr = pd.read_csv(corr_csv)
    matrix = corr.pivot(index="gene", columns="celltype_feature", values="spearman_r")
    return matrix.reindex(MARKER_GENES).fillna(0.0)


def assign_unique_celltypes(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign one unique cell type to each marker by maximizing abs(Spearman r)."""
    cost = -matrix.abs().to_numpy()
    row_idx, col_idx = linear_sum_assignment(cost)
    assignment = (
        pd.DataFrame(
            {
                "gene": matrix.index[row_idx],
                "assigned_celltype": matrix.columns[col_idx],
                "spearman_r": matrix.to_numpy()[row_idx, col_idx],
                "abs_spearman_r": np.abs(matrix.to_numpy()[row_idx, col_idx]),
            }
        )
        .set_index("gene")
        .reindex(MARKER_GENES)
        .reset_index()
    )
    selected_columns = assignment["assigned_celltype"].tolist()
    diagonal_matrix = matrix.loc[MARKER_GENES, selected_columns]
    return diagonal_matrix, assignment


def order_all_celltypes(matrix: pd.DataFrame, assignment: pd.DataFrame) -> pd.DataFrame:
    """Keep all cell types while placing assigned marker-celltype pairs first."""
    assigned_columns = assignment["assigned_celltype"].tolist()
    remaining_columns = [col for col in matrix.columns if col not in assigned_columns]
    return matrix.loc[MARKER_GENES, assigned_columns + remaining_columns]


def plot_diagonal_heatmap(matrix: pd.DataFrame, assignment: pd.DataFrame, out_dir: Path, prefix: str) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    sns.heatmap(
        matrix,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.4,
        linecolor="white",
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 6},
        cbar_kws={"label": "Spearman r"},
        ax=ax,
    )
    for i in range(len(MARKER_GENES)):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=False, edgecolor="#111111", lw=1.2))
    ax.set_xlabel("Assigned cell type", fontsize=8, fontweight="bold")
    ax.set_ylabel("ML marker", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save_figure(fig, out_dir / f"{prefix}_heatmap")

    diag = pd.DataFrame(
        np.nan,
        index=matrix.index,
        columns=matrix.columns,
    )
    for i, gene in enumerate(matrix.index):
        diag.iloc[i, i] = matrix.iloc[i, i]

    fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    sns.heatmap(
        diag,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.4,
        linecolor="white",
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 6},
        cbar_kws={"label": "Spearman r"},
        ax=ax,
    )
    ax.set_facecolor("#f3f3f3")
    ax.set_xlabel("Assigned cell type", fontsize=8, fontweight="bold")
    ax.set_ylabel("ML marker", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save_figure(fig, out_dir / f"{prefix}_only")


def plot_all_celltype_heatmap(
    matrix: pd.DataFrame,
    assignment: pd.DataFrame,
    out_dir: Path,
    prefix: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.2), constrained_layout=True)
    sns.heatmap(
        matrix,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "Spearman r"},
        ax=ax,
    )
    assigned_columns = assignment["assigned_celltype"].tolist()
    for row_i, celltype in enumerate(assigned_columns):
        col_i = matrix.columns.get_loc(celltype)
        ax.add_patch(plt.Rectangle((col_i, row_i), 1, 1, fill=False, edgecolor="#111111", lw=1.2))
    ax.axvline(len(assigned_columns), color="#111111", lw=0.8)
    ax.set_xlabel("Cell type", fontsize=8, fontweight="bold")
    ax.set_ylabel("ML marker", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save_figure(fig, out_dir / f"{prefix}_all_celltypes")


def main() -> None:
    args = parse_args()
    corr_csv = Path(args.corr_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_plot_style()
    matrix = load_matrix(corr_csv)
    diagonal_matrix, assignment = assign_unique_celltypes(matrix)
    all_celltype_matrix = order_all_celltypes(matrix, assignment)
    assignment.to_csv(out_dir / f"{args.prefix}_assigned_pairs.csv", index=False)
    diagonal_matrix.to_csv(out_dir / f"{args.prefix}_selected_matrix.csv")
    all_celltype_matrix.to_csv(out_dir / f"{args.prefix}_all_celltypes_matrix.csv")
    plot_diagonal_heatmap(diagonal_matrix, assignment, out_dir, args.prefix)
    plot_all_celltype_heatmap(all_celltype_matrix, assignment, out_dir, args.prefix)

    print("Done.")
    print(f"Assigned pairs: {(out_dir / f'{args.prefix}_assigned_pairs.csv').resolve()}")
    print(f"Figure: {(out_dir / f'{args.prefix}_heatmap.png').resolve()}")
    print(f"Diagonal-only figure: {(out_dir / f'{args.prefix}_only.png').resolve()}")
    print(f"All-celltype figure: {(out_dir / f'{args.prefix}_all_celltypes.png').resolve()}")


if __name__ == "__main__":
    main()
