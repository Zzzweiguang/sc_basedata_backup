#!/usr/bin/env python
"""Sample-level biomarker discovery adapted from Wang et al. 2022.

Reference:
Wang J. et al. Identification of immune cell infiltration and diagnostic
biomarkers in unstable atherosclerotic plaques by integrated bioinformatics
analysis and machine learning. Front Immunol. 2022;13:956078.

The original paper used bulk plaque expression, CIBERSORT immune infiltration,
DEGs, LASSO, random forest, ROC analysis, and scRNA/clinical validation.
This script adapts the workflow to a labelled single-cell AnnData object by
aggregating cells to sample-level pseudobulk expression before ML.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ARTICLE_DIAGNOSTIC_GENES = ["CD68", "PAM", "IGFBP6"]
ARTICLE_SIGNAL_GENES = {
    "macrophage_m1_instability": ["CD68", "CD80", "CD86", "IL1B", "TNF", "NLRP3", "CCL2", "CXCL8"],
    "antigen_inflammation_lipid": ["CD74", "B2M", "HLA-DRA", "MMP9", "CTSL", "IFI30", "CD36", "APOE"],
    "smc_contractile_stability": ["TAGLN", "ACAT2", "MYH10", "MYH11"],
    "cd8_t_nk_protective": ["CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "IFNG", "KLRD1"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        default="/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/0521_no_Basophil/output_allhuman/scPoli_concat_level3_marker_all_metadata.h5ad",
        help="Input AnnData .h5ad file.",
    )
    parser.add_argument(
        "--out-dir",
        default="output_article_biomarker_pipeline",
        help="Output directory.",
    )
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="Plaque_type")
    parser.add_argument("--celltype-col", default="cell_type_level1_corrected")
    parser.add_argument("--layer", default="raw_decontXcounts")
    parser.add_argument("--target-sum", type=float, default=1e6, help="Pseudobulk CPM target sum.")
    parser.add_argument("--min-cells-per-sample", type=int, default=50)
    parser.add_argument("--max-ml-genes", type=int, default=2000)
    parser.add_argument("--min-ml-genes", type=int, default=50)
    parser.add_argument("--deg-padj", type=float, default=0.05)
    parser.add_argument("--deg-log2fc", type=float, default=0.25)
    parser.add_argument("--rf-top-n", type=int, default=50)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=8)
    return parser.parse_args()


def ensure_plot_style() -> None:
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"
    sns.set_theme(style="ticks")


def clean_label_series(values: pd.Series) -> pd.Series:
    labels = values.astype(object).where(values.notna(), "unknown").astype(str)
    labels = labels.str.strip().str.lower()
    return labels.replace(
        {
            "": "unknown",
            "na": "unknown",
            "nan": "unknown",
            "none": "unknown",
            "<na>": "unknown",
        }
    )


def read_obs_var(h5ad_path: str) -> tuple[pd.DataFrame, pd.DataFrame, int, int]:
    print(f"Reading metadata from {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    var = adata.var.copy()
    n_obs, n_vars = adata.shape
    adata.file.close()
    return obs, var, n_obs, n_vars


def resolve_gene_names(var: pd.DataFrame) -> pd.Index:
    if "original_gene_names" in var.columns:
        gene_names = var["original_gene_names"].astype(str).str.strip()
    else:
        gene_names = pd.Index(var.index.astype(str)).to_series(index=var.index)
    gene_names = gene_names.replace({"": np.nan, "nan": np.nan, "None": np.nan})
    fallback = pd.Index(var.index.astype(str)).to_series(index=var.index)
    gene_names = gene_names.fillna(fallback)
    return pd.Index(gene_names.astype(str))


def prepare_labelled_obs(
    obs: pd.DataFrame,
    sample_col: str,
    label_col: str,
    celltype_col: str,
    min_cells_per_sample: int,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    required = [sample_col, label_col, celltype_col]
    missing = [col for col in required if col not in obs.columns]
    if missing:
        raise ValueError(f"Missing required obs columns: {missing}")

    obs = obs.copy()
    obs["_label_clean"] = clean_label_series(obs[label_col])
    labelled_mask = obs["_label_clean"].isin(["stable", "unstable"])
    labelled_obs = obs.loc[labelled_mask].copy()
    if labelled_obs.empty:
        raise ValueError("No labelled cells found after keeping stable/unstable Plaque_type values.")

    label_nunique = labelled_obs.groupby(sample_col, observed=True)["_label_clean"].nunique()
    conflict_samples = label_nunique[label_nunique > 1].index.astype(str).tolist()
    if conflict_samples:
        preview = ", ".join(conflict_samples[:20])
        raise ValueError(f"Samples contain both stable and unstable labels: {preview}")

    sample_cell_counts = labelled_obs.groupby(sample_col, observed=True).size().rename("n_labelled_cells")
    keep_samples = sample_cell_counts[sample_cell_counts >= min_cells_per_sample].index
    dropped = sample_cell_counts[sample_cell_counts < min_cells_per_sample]
    if not dropped.empty:
        print(f"Dropping {dropped.size} labelled samples with < {min_cells_per_sample} cells.")
    labelled_obs = labelled_obs[labelled_obs[sample_col].isin(keep_samples)].copy()

    sample_label = (
        labelled_obs[[sample_col, "_label_clean"]]
        .drop_duplicates()
        .set_index(sample_col)["_label_clean"]
        .sort_index()
    )
    y = sample_label.map({"stable": 0, "unstable": 1}).astype(int)
    if y.nunique() != 2:
        raise ValueError(f"Need both stable and unstable samples. Current sample labels: {sample_label.value_counts().to_dict()}")
    return labelled_obs, sample_label, y


def save_sample_metadata(
    out_dir: Path,
    labelled_obs: pd.DataFrame,
    sample_col: str,
    sample_label: pd.Series,
) -> None:
    meta = labelled_obs.groupby(sample_col, observed=True).agg(n_labelled_cells=(sample_col, "size"))
    for col in ["dataset", "tissue", "species", "symptoms", "gender", "age"]:
        if col in labelled_obs.columns:
            meta[col] = labelled_obs.groupby(sample_col, observed=True)[col].agg(lambda x: ";".join(sorted(map(str, pd.unique(x)))))
    meta["label"] = sample_label
    meta["label_numeric"] = sample_label.map({"stable": 0, "unstable": 1}).astype(int)
    meta.sort_index().to_csv(out_dir / "sample_metadata.csv")


def build_celltype_proportions(
    labelled_obs: pd.DataFrame,
    sample_col: str,
    celltype_col: str,
) -> pd.DataFrame:
    counts = pd.crosstab(labelled_obs[sample_col], labelled_obs[celltype_col])
    props = counts.div(counts.sum(axis=1), axis=0)
    props.columns = [f"{celltype_col}__{col}" for col in props.columns.astype(str)]
    return props.sort_index()


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def differential_table(features: pd.DataFrame, y: pd.Series, feature_name: str) -> pd.DataFrame:
    y = y.loc[features.index]
    stable = features.loc[y == 0]
    unstable = features.loc[y == 1]
    stable_values = stable.to_numpy(dtype=float)
    unstable_values = unstable.to_numpy(dtype=float)

    stat, pvalue = stats.ttest_ind(
        unstable_values,
        stable_values,
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    padj = benjamini_hochberg(pvalue)
    mean_stable = np.nanmean(stable_values, axis=0)
    mean_unstable = np.nanmean(unstable_values, axis=0)
    diff = mean_unstable - mean_stable

    return pd.DataFrame(
        {
            feature_name: features.columns.astype(str),
            "mean_stable": mean_stable,
            "mean_unstable": mean_unstable,
            "unstable_minus_stable": diff,
            "t_stat": stat,
            "pvalue": pvalue,
            "padj": padj,
        }
    ).sort_values(["padj", "pvalue", "unstable_minus_stable"], ascending=[True, True, False])


def aggregate_csr_layer_by_sample(
    h5ad_path: str,
    layer: str,
    n_obs: int,
    n_vars: int,
    sample_codes_by_row: np.ndarray,
    n_samples: int,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    if layer == "X":
        group_path = "X"
    else:
        group_path = f"layers/{layer}"

    sample_gene_counts = np.zeros((n_samples, n_vars), dtype=np.float64)
    sample_cell_counts = np.zeros(n_samples, dtype=np.int64)

    with h5py.File(h5ad_path, "r") as handle:
        if group_path not in handle:
            raise ValueError(f"H5AD does not contain sparse matrix group: {group_path}")
        group = handle[group_path]
        encoding = group.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise NotImplementedError(f"Only csr_matrix layers are supported. {group_path} encoding={encoding!r}")

        indptr_ds = group["indptr"]
        indices_ds = group["indices"]
        data_ds = group["data"]
        n_blocks = math.ceil(n_obs / chunk_size)

        for block_id, start in enumerate(range(0, n_obs, chunk_size), start=1):
            end = min(start + chunk_size, n_obs)
            row_sample_codes = sample_codes_by_row[start:end]
            labelled = row_sample_codes >= 0
            if not labelled.any():
                continue

            data_start = int(indptr_ds[start])
            data_end = int(indptr_ds[end])
            local_indptr = np.asarray(indptr_ds[start : end + 1], dtype=np.int64) - data_start
            block = sparse.csr_matrix(
                (
                    np.asarray(data_ds[data_start:data_end], dtype=np.float64),
                    np.asarray(indices_ds[data_start:data_end], dtype=np.int32),
                    local_indptr,
                ),
                shape=(end - start, n_vars),
            )
            block = block[labelled]
            labelled_codes = row_sample_codes[labelled]

            for sample_code in np.unique(labelled_codes):
                rows = labelled_codes == sample_code
                sample_gene_counts[sample_code] += np.asarray(block[rows].sum(axis=0)).ravel()
                sample_cell_counts[sample_code] += int(rows.sum())

            if block_id == 1 or block_id == n_blocks or block_id % 25 == 0:
                done = min(end, n_obs)
                print(f"  processed rows {done:,}/{n_obs:,}")

    return sample_gene_counts, sample_cell_counts


def build_pseudobulk_expression(
    h5ad_path: str,
    layer: str,
    obs: pd.DataFrame,
    labelled_obs: pd.DataFrame,
    sample_col: str,
    sample_order: pd.Index,
    gene_names: pd.Index,
    n_obs: int,
    n_vars: int,
    target_sum: float,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_to_code = {sample: idx for idx, sample in enumerate(sample_order)}
    sample_codes_by_row = np.full(n_obs, -1, dtype=np.int32)
    row_positions = obs.index.get_indexer(labelled_obs.index)
    if (row_positions < 0).any():
        raise ValueError("Unable to align labelled obs rows back to full obs index.")
    sample_codes_by_row[row_positions] = labelled_obs[sample_col].map(sample_to_code).to_numpy(dtype=np.int32)

    print(f"Aggregating layer {layer!r} to sample-level pseudobulk counts.")
    sample_gene_counts, sample_cell_counts = aggregate_csr_layer_by_sample(
        h5ad_path=h5ad_path,
        layer=layer,
        n_obs=n_obs,
        n_vars=n_vars,
        sample_codes_by_row=sample_codes_by_row,
        n_samples=len(sample_order),
        chunk_size=chunk_size,
    )

    counts = pd.DataFrame(sample_gene_counts, index=sample_order, columns=gene_names)
    if counts.columns.duplicated().any():
        n_dup = int(counts.columns.duplicated().sum())
        print(f"Collapsing {n_dup} duplicated gene-symbol columns by summing counts.")
        counts = counts.T.groupby(level=0).sum().T

    library_size = counts.sum(axis=1)
    log2cpm = np.log2(counts.div(library_size.replace(0, np.nan), axis=0) * target_sum + 1.0)
    log2cpm = log2cpm.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    cell_counts = pd.DataFrame({"n_cells_aggregated": sample_cell_counts}, index=sample_order)
    return log2cpm, cell_counts


def make_deg_table(log2cpm: pd.DataFrame, y: pd.Series, target_sum: float) -> pd.DataFrame:
    y = y.loc[log2cpm.index]
    stable = log2cpm.loc[y == 0]
    unstable = log2cpm.loc[y == 1]
    stat, pvalue = stats.ttest_ind(
        unstable.to_numpy(dtype=float),
        stable.to_numpy(dtype=float),
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    padj = benjamini_hochberg(pvalue)
    mean_log2_stable = stable.mean(axis=0).to_numpy(dtype=float)
    mean_log2_unstable = unstable.mean(axis=0).to_numpy(dtype=float)
    diff_log2 = mean_log2_unstable - mean_log2_stable
    mean_cpm_stable = np.maximum(np.power(2.0, mean_log2_stable) - 1.0, 0.0)
    mean_cpm_unstable = np.maximum(np.power(2.0, mean_log2_unstable) - 1.0, 0.0)
    log2fc = np.log2((mean_cpm_unstable + 1.0) / (mean_cpm_stable + 1.0))

    table = pd.DataFrame(
        {
            "gene": log2cpm.columns.astype(str),
            "mean_log2cpm_stable": mean_log2_stable,
            "mean_log2cpm_unstable": mean_log2_unstable,
            "mean_log2cpm_diff_unstable_minus_stable": diff_log2,
            "mean_cpm_stable": mean_cpm_stable,
            "mean_cpm_unstable": mean_cpm_unstable,
            "log2FC_mean_cpm": log2fc,
            "t_stat": stat,
            "pvalue": pvalue,
            "padj": padj,
            "target_sum": target_sum,
        }
    )
    table["direction"] = np.where(table["log2FC_mean_cpm"] >= 0, "unstable_up", "stable_up")
    return table.sort_values(["padj", "pvalue", "log2FC_mean_cpm"], ascending=[True, True, False])


def select_ml_genes(
    deg: pd.DataFrame,
    log2cpm: pd.DataFrame,
    deg_padj: float,
    deg_log2fc: float,
    min_ml_genes: int,
    max_ml_genes: int,
) -> tuple[list[str], pd.DataFrame]:
    nonzero_var = log2cpm.var(axis=0) > 0
    usable_genes = set(nonzero_var[nonzero_var].index.astype(str))
    significant = deg[
        (deg["gene"].isin(usable_genes))
        & (deg["padj"] <= deg_padj)
        & (deg["log2FC_mean_cpm"].abs() >= deg_log2fc)
    ].copy()
    significant["candidate_source"] = "deg_threshold"

    if significant.shape[0] < min_ml_genes:
        fallback = deg[deg["gene"].isin(usable_genes)].head(max_ml_genes).copy()
        fallback["candidate_source"] = "top_pvalue_fallback"
        candidate = fallback
    else:
        candidate = significant.head(max_ml_genes)

    article_genes = sorted(
        set(ARTICLE_DIAGNOSTIC_GENES + [gene for genes in ARTICLE_SIGNAL_GENES.values() for gene in genes])
        & set(log2cpm.columns.astype(str))
    )
    article_rows = deg[deg["gene"].isin(article_genes)].copy()
    article_rows["candidate_source"] = "article_signal_added"
    candidate = pd.concat([candidate, article_rows], axis=0)
    candidate = candidate.drop_duplicates("gene")

    if candidate.shape[0] > max_ml_genes:
        article_mask = candidate["gene"].isin(article_genes)
        article_part = candidate[article_mask]
        non_article_part = candidate[~article_mask].head(max_ml_genes - article_part.shape[0])
        candidate = pd.concat([non_article_part, article_part], axis=0)

    return candidate["gene"].astype(str).tolist(), candidate


def cv_splits_for_y(y: pd.Series, random_state: int) -> StratifiedKFold:
    min_class_count = int(y.value_counts().min())
    if min_class_count < 2:
        raise ValueError("Each class needs at least two samples for cross-validation.")
    return StratifiedKFold(n_splits=min(5, min_class_count), shuffle=True, random_state=random_state)


def fit_lasso_and_rf(
    x: pd.DataFrame,
    y: pd.Series,
    cv: StratifiedKFold,
    random_state: int,
    n_jobs: int,
) -> tuple[pd.DataFrame, dict]:
    y = y.loc[x.index]
    info: dict[str, object] = {"n_samples": int(x.shape[0]), "n_features": int(x.shape[1])}

    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    lasso_model = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegressionCV(
                    Cs=np.logspace(-3, 1, 20),
                    cv=cv,
                    penalty="l1",
                    solver="saga",
                    scoring="roc_auc",
                    class_weight="balanced",
                    max_iter=20000,
                    random_state=random_state,
                    n_jobs=n_jobs,
                    refit=True,
                ),
            ),
        ]
    )
    lasso_model.fit(x, y)
    lasso_coef = lasso_model.named_steps["clf"].coef_[0]
    info["lasso_selected_c"] = float(lasso_model.named_steps["clf"].C_[0])
    info["lasso_refit_note"] = "LogisticRegressionCV l1"

    if np.count_nonzero(np.abs(lasso_coef) > 1e-8) == 0:
        fallback_model = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        C=10.0,
                        penalty="l1",
                        solver="saga",
                        class_weight="balanced",
                        max_iter=20000,
                        random_state=random_state,
                    ),
                ),
            ]
        )
        fallback_model.fit(x, y)
        lasso_coef = fallback_model.named_steps["clf"].coef_[0]
        info["lasso_refit_note"] = "CV selected zero genes; refit l1 with C=10.0 for ranking"

    rf_model = RandomForestClassifier(
        n_estimators=1000,
        random_state=random_state,
        class_weight="balanced_subsample",
        max_features="sqrt",
        n_jobs=n_jobs,
    )
    rf_model.fit(x, y)

    ranking = pd.DataFrame(
        {
            "gene": x.columns.astype(str),
            "lasso_coef": lasso_coef,
            "abs_lasso_coef": np.abs(lasso_coef),
            "rf_importance": rf_model.feature_importances_,
        }
    )
    ranking["selected_by_lasso"] = ranking["abs_lasso_coef"] > 1e-8
    return ranking, info


def add_biomarker_annotations(
    ranking: pd.DataFrame,
    deg: pd.DataFrame,
    rf_top_n: int,
) -> pd.DataFrame:
    deg_cols = [
        "gene",
        "mean_log2cpm_stable",
        "mean_log2cpm_unstable",
        "mean_log2cpm_diff_unstable_minus_stable",
        "mean_cpm_stable",
        "mean_cpm_unstable",
        "log2FC_mean_cpm",
        "pvalue",
        "padj",
        "direction",
    ]
    result = ranking.merge(deg[deg_cols], on="gene", how="left")
    rf_rank = result["rf_importance"].rank(ascending=False, method="min")
    result["rf_rank"] = rf_rank.astype(int)
    result["selected_by_rf_top"] = result["rf_rank"] <= rf_top_n
    result["selected_by_both_ml"] = result["selected_by_lasso"] & result["selected_by_rf_top"]
    result["selected_by_article_diagnostic"] = result["gene"].isin(ARTICLE_DIAGNOSTIC_GENES)

    signal_map = {}
    for signal, genes in ARTICLE_SIGNAL_GENES.items():
        for gene in genes:
            signal_map.setdefault(gene, []).append(signal)
    result["article_signal_group"] = result["gene"].map(lambda gene: ";".join(signal_map.get(gene, [])))

    p_score = -np.log10(result["padj"].fillna(1.0).clip(lower=1e-300))
    for col, new_col in [("abs_lasso_coef", "_lasso_scaled"), ("rf_importance", "_rf_scaled")]:
        values = result[col].fillna(0.0)
        max_value = values.max()
        result[new_col] = values / max_value if max_value > 0 else 0.0
    result["_p_scaled"] = p_score / p_score.max() if p_score.max() > 0 else 0.0
    result["combined_rank_score"] = result["_lasso_scaled"] + result["_rf_scaled"] + result["_p_scaled"]
    result = result.drop(columns=["_lasso_scaled", "_rf_scaled", "_p_scaled"])
    return result.sort_values(
        ["selected_by_both_ml", "selected_by_article_diagnostic", "combined_rank_score", "rf_importance"],
        ascending=[False, False, False, False],
    )


def model_metrics(y_true: pd.Series, y_pred: np.ndarray, proba: np.ndarray) -> dict[str, float]:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    return {
        "n_samples": int(len(y_true)),
        "n_stable": int((y_true == 0).sum()),
        "n_unstable": int((y_true == 1).sum()),
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision_unstable": precision_score(y_true, y_pred, zero_division=0),
        "recall_unstable_sensitivity": recall_score(y_true, y_pred, zero_division=0),
        "specificity_stable": specificity,
        "f1_unstable": f1_score(y_true, y_pred, zero_division=0),
        "roc_auc": roc_auc_score(y_true, proba) if pd.Series(y_true).nunique() == 2 else np.nan,
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def evaluate_cv_models(
    x: pd.DataFrame,
    y: pd.Series,
    selected_genes: list[str],
    out_dir: Path,
    cv: StratifiedKFold,
    random_state: int,
    n_jobs: int,
) -> pd.DataFrame:
    y = y.loc[x.index]
    model_specs = {
        "logistic_l2_candidate_genes": (
            x,
            Pipeline(
                [
                    ("scaler", StandardScaler()),
                    (
                        "clf",
                        LogisticRegression(
                            penalty="l2",
                            C=1.0,
                            class_weight="balanced",
                            max_iter=10000,
                            random_state=random_state,
                        ),
                    ),
                ]
            ),
        ),
        "random_forest_candidate_genes": (
            x,
            RandomForestClassifier(
                n_estimators=1000,
                random_state=random_state,
                class_weight="balanced_subsample",
                max_features="sqrt",
                n_jobs=n_jobs,
            ),
        ),
    }
    if len(selected_genes) >= 2:
        model_specs["logistic_l2_selected_biomarkers"] = (
            x.loc[:, selected_genes],
            Pipeline(
                [
                    ("scaler", StandardScaler()),
                    (
                        "clf",
                        LogisticRegression(
                            penalty="l2",
                            C=1.0,
                            class_weight="balanced",
                            max_iter=10000,
                            random_state=random_state,
                        ),
                    ),
                ]
            ),
        )

    metric_rows = []
    prediction_rows = []
    for model_name, (model_x, model) in model_specs.items():
        pred = cross_val_predict(model, model_x, y, cv=cv, method="predict")
        proba = cross_val_predict(model, model_x, y, cv=cv, method="predict_proba")[:, 1]
        metrics = model_metrics(y, pred, proba)
        metrics["model"] = model_name
        metrics["n_features"] = int(model_x.shape[1])
        metrics["cv_note"] = "CV uses feature set selected from full data; use external validation for final diagnostic claims."
        metric_rows.append(metrics)
        pred_df = pd.DataFrame(
            {
                "sample_id": model_x.index,
                "model": model_name,
                "true_label_numeric": y.values,
                "true_label": y.map({0: "stable", 1: "unstable"}).values,
                "prob_unstable": proba,
                "pred_label_numeric": pred,
                "pred_label": pd.Series(pred).map({0: "stable", 1: "unstable"}).values,
            }
        )
        prediction_rows.append(pred_df)

    predictions = pd.concat(prediction_rows, axis=0, ignore_index=True)
    predictions.to_csv(out_dir / "sample_diagnostic_cv_predictions.csv", index=False)
    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(out_dir / "sample_diagnostic_cv_metrics.csv", index=False)
    return metrics_df


def fit_immune_ml_ranking(
    proportions: pd.DataFrame,
    y: pd.Series,
    cv: StratifiedKFold,
    random_state: int,
    n_jobs: int,
) -> pd.DataFrame:
    ranking, _ = fit_lasso_and_rf(proportions, y.loc[proportions.index], cv, random_state, n_jobs)
    ranking = ranking.rename(columns={"gene": "celltype_feature"})
    ranking["rf_rank"] = ranking["rf_importance"].rank(ascending=False, method="min").astype(int)
    return ranking.sort_values(["selected_by_lasso", "rf_importance"], ascending=[False, False])


def save_json(data: dict, path: Path) -> None:
    with path.open("w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def save_plots(out_dir: Path, deg: pd.DataFrame, biomarker_table: pd.DataFrame, log2cpm: pd.DataFrame, y: pd.Series) -> None:
    ensure_plot_style()
    plot_volcano(out_dir, deg, biomarker_table)
    plot_rf_importance(out_dir, biomarker_table)
    plot_lasso(out_dir, biomarker_table)
    plot_heatmap(out_dir, biomarker_table, log2cpm, y)


def save_figure(fig: plt.Figure, out_path: Path) -> None:
    fig.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_volcano(out_dir: Path, deg: pd.DataFrame, biomarker_table: pd.DataFrame) -> None:
    df = deg.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["is_sig"] = (df["padj"] <= 0.05) & (df["log2FC_mean_cpm"].abs() >= 0.25)
    fig, ax = plt.subplots(figsize=(3.6, 3.2), constrained_layout=True)
    colors = np.where(df["is_sig"], "#C74440", "#8A8F98")
    ax.scatter(df["log2FC_mean_cpm"], df["neg_log10_padj"], c=colors, s=5, alpha=0.65, linewidths=0)
    ax.axvline(0, color="#333333", lw=0.6)
    ax.axhline(-np.log10(0.05), color="#333333", lw=0.6, ls="--")
    top_labels = (
        biomarker_table[biomarker_table["selected_by_both_ml"] | biomarker_table["selected_by_article_diagnostic"]]
        .head(15)["gene"]
        .tolist()
    )
    label_df = df[df["gene"].isin(top_labels)]
    texts = []
    for _, row in label_df.iterrows():
        texts.append(ax.text(row["log2FC_mean_cpm"], row["neg_log10_padj"], row["gene"], fontsize=5))
    try:
        from adjustText import adjust_text

        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.4), lim=300)
    except Exception:
        pass
    ax.set_xlabel("log2FC mean CPM (unstable / stable)", fontsize=8, fontweight="bold")
    ax.set_ylabel("-log10 adjusted P", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "deg_volcano")


def plot_rf_importance(out_dir: Path, biomarker_table: pd.DataFrame) -> None:
    top = biomarker_table.sort_values("rf_importance", ascending=False).head(30).iloc[::-1]
    fig, ax = plt.subplots(figsize=(4.0, 4.8), constrained_layout=True)
    ax.barh(top["gene"], top["rf_importance"], color="#4B7A9F")
    ax.set_xlabel("Random forest importance", fontsize=8, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=7)
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "random_forest_top_genes")


def plot_lasso(out_dir: Path, biomarker_table: pd.DataFrame) -> None:
    top = biomarker_table[biomarker_table["selected_by_lasso"]].copy()
    if top.empty:
        top = biomarker_table.sort_values("abs_lasso_coef", ascending=False).head(30)
    top = top.sort_values("lasso_coef").tail(40)
    fig, ax = plt.subplots(figsize=(4.2, max(2.6, 0.14 * len(top))), constrained_layout=True)
    colors = np.where(top["lasso_coef"] >= 0, "#C74440", "#2F6B4F")
    ax.barh(top["gene"], top["lasso_coef"], color=colors)
    ax.axvline(0, color="#333333", lw=0.6)
    ax.set_xlabel("LASSO logistic coefficient", fontsize=8, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=7)
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "lasso_selected_gene_coefficients")


def plot_heatmap(out_dir: Path, biomarker_table: pd.DataFrame, log2cpm: pd.DataFrame, y: pd.Series) -> None:
    genes = biomarker_table[biomarker_table["selected_by_both_ml"]]["gene"].tolist()
    if len(genes) < 2:
        genes = biomarker_table.head(25)["gene"].tolist()
    genes = [gene for gene in genes if gene in log2cpm.columns][:40]
    if len(genes) < 2:
        return
    sample_order = y.sort_values().index
    data = log2cpm.loc[sample_order, genes]
    z = (data - data.mean(axis=0)) / data.std(axis=0, ddof=0).replace(0, np.nan)
    z = z.fillna(0.0).T
    fig, ax = plt.subplots(figsize=(max(4.0, 0.12 * z.shape[1]), max(2.8, 0.16 * z.shape[0])), constrained_layout=True)
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0, cbar_kws={"label": "Z-score"}, xticklabels=False, yticklabels=True)
    ax.set_xlabel("Samples ordered by Plaque_type", fontsize=8, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=6)
    save_figure(fig, out_dir / "selected_biomarker_heatmap")


def write_methods_summary(
    out_dir: Path,
    args: argparse.Namespace,
    sample_label: pd.Series,
    candidate_genes: list[str],
    biomarker_table: pd.DataFrame,
    metrics: pd.DataFrame,
) -> None:
    core = biomarker_table[biomarker_table["selected_by_both_ml"]]["gene"].tolist()
    if not core:
        core = biomarker_table.head(20)["gene"].tolist()
    lines = [
        "# Article-adapted biomarker pipeline summary",
        "",
        "Reference: Wang J. et al., Front Immunol. 2022, doi:10.3389/fimmu.2022.956078.",
        "",
        "Adaptation used here:",
        "- Plaque_type values stable/unstable were treated as labelled samples; unknown/missing labels were excluded.",
        "- Single cells were aggregated by sample before differential expression and machine learning.",
        "- Annotated cell-type proportions were used as the single-cell analogue of immune infiltration.",
        "- Pseudobulk log2(CPM+1) expression was used for DEG, LASSO logistic regression, and random forest.",
        "",
        "Inputs:",
        f"- h5ad: `{args.h5ad}`",
        f"- layer: `{args.layer}`",
        f"- sample column: `{args.sample_col}`",
        f"- label column: `{args.label_col}`",
        f"- cell-type column: `{args.celltype_col}`",
        "",
        "Labelled sample counts:",
        sample_label.value_counts().to_string(),
        "",
        f"ML candidate genes: {len(candidate_genes)}",
        "",
        "Top data-derived biomarkers:",
        ", ".join(core[:30]),
        "",
        "Cross-validation metrics:",
        metrics.to_string(index=False),
        "",
        "Caution: this is discovery on the provided dataset. Diagnostic performance should be confirmed on an external cohort.",
        "",
    ]
    (out_dir / "README_article_biomarker_pipeline.md").write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    save_json(vars(args), out_dir / "run_parameters.json")
    obs, var, n_obs, n_vars = read_obs_var(args.h5ad)
    gene_names = resolve_gene_names(var)
    labelled_obs, sample_label, y = prepare_labelled_obs(
        obs=obs,
        sample_col=args.sample_col,
        label_col=args.label_col,
        celltype_col=args.celltype_col,
        min_cells_per_sample=args.min_cells_per_sample,
    )
    sample_order = y.index
    print("Labelled sample counts:")
    print(sample_label.value_counts())
    print("Labelled cell counts:")
    print(labelled_obs["_label_clean"].value_counts())

    save_sample_metadata(out_dir, labelled_obs, args.sample_col, sample_label)

    candidate_celltype_cols = [
        args.celltype_col,
        "cell_type_level1_corrected",
        "cell_type_level2",
        "cell_type_level3",
    ]
    seen_cols = set()
    for col in candidate_celltype_cols:
        if col in seen_cols or col not in labelled_obs.columns:
            continue
        seen_cols.add(col)
        props = build_celltype_proportions(labelled_obs, args.sample_col, col).reindex(sample_order).fillna(0.0)
        props.to_csv(out_dir / f"immune_celltype_proportions__{col}.csv")
        diff = differential_table(props, y, feature_name="celltype_feature")
        diff.to_csv(out_dir / f"immune_celltype_differential__{col}.csv", index=False)

    main_props = build_celltype_proportions(labelled_obs, args.sample_col, args.celltype_col).reindex(sample_order).fillna(0.0)
    cv = cv_splits_for_y(y, args.random_state)
    immune_ranking = fit_immune_ml_ranking(main_props, y, cv, args.random_state, args.n_jobs)
    immune_ranking.to_csv(out_dir / "immune_infiltration_ml_ranking.csv", index=False)

    log2cpm, agg_cell_counts = build_pseudobulk_expression(
        h5ad_path=args.h5ad,
        layer=args.layer,
        obs=obs,
        labelled_obs=labelled_obs,
        sample_col=args.sample_col,
        sample_order=sample_order,
        gene_names=gene_names,
        n_obs=n_obs,
        n_vars=n_vars,
        target_sum=args.target_sum,
        chunk_size=args.chunk_size,
    )
    log2cpm.to_csv(out_dir / "pseudobulk_log2cpm.csv.gz")
    agg_cell_counts.to_csv(out_dir / "pseudobulk_aggregated_cell_counts.csv")

    deg = make_deg_table(log2cpm, y, args.target_sum)
    deg.to_csv(out_dir / "pseudobulk_deg_stable_vs_unstable.csv", index=False)

    candidate_genes, candidate_table = select_ml_genes(
        deg=deg,
        log2cpm=log2cpm,
        deg_padj=args.deg_padj,
        deg_log2fc=args.deg_log2fc,
        min_ml_genes=args.min_ml_genes,
        max_ml_genes=args.max_ml_genes,
    )
    candidate_table.to_csv(out_dir / "ml_candidate_genes.csv", index=False)
    x_ml = log2cpm.loc[:, candidate_genes]

    ranking, lasso_info = fit_lasso_and_rf(x_ml, y, cv, args.random_state, args.n_jobs)
    save_json(lasso_info, out_dir / "ml_fit_info.json")
    biomarker_table = add_biomarker_annotations(ranking, deg, args.rf_top_n)
    biomarker_table.to_csv(out_dir / "biomarker_lasso_rf_ranked.csv", index=False)

    core = biomarker_table[biomarker_table["selected_by_both_ml"]]["gene"].tolist()
    if len(core) < 2:
        core = biomarker_table[biomarker_table["selected_by_lasso"] | biomarker_table["selected_by_rf_top"]].head(30)["gene"].tolist()
    metrics = evaluate_cv_models(x_ml, y, core, out_dir, cv, args.random_state, args.n_jobs)

    save_plots(out_dir, deg, biomarker_table, log2cpm, y)
    write_methods_summary(out_dir, args, sample_label, candidate_genes, biomarker_table, metrics)

    print("\nDone.")
    print(f"Output directory: {out_dir.resolve()}")
    print("Top biomarkers:")
    print(biomarker_table.head(20)[["gene", "selected_by_both_ml", "selected_by_lasso", "selected_by_rf_top", "log2FC_mean_cpm", "padj", "lasso_coef", "rf_importance"]].to_string(index=False))
    print("\nCV metrics:")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
