#!/usr/bin/env python
"""Strict nested sample-level biomarker validation.

Each outer fold performs all supervised steps on the training samples only:
DEG -> LASSO -> random forest -> biomarker selection -> final classifier.
The held-out samples are used only for evaluation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from run_article_biomarker_pipeline import (
    add_biomarker_annotations,
    cv_splits_for_y,
    evaluate_cv_models,
    fit_lasso_and_rf,
    make_deg_table,
    model_metrics,
    save_figure,
    save_json,
    select_ml_genes,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        default="output_article_biomarker_pipeline",
        help="Directory containing pseudobulk_log2cpm.csv.gz and sample_metadata.csv.",
    )
    parser.add_argument(
        "--out-dir",
        default="output_article_biomarker_pipeline/nested_cv",
        help="Output directory for strict nested CV results.",
    )
    parser.add_argument("--outer-splits", type=int, default=5)
    parser.add_argument("--deg-padj", type=float, default=0.05)
    parser.add_argument("--deg-log2fc", type=float, default=0.25)
    parser.add_argument("--min-ml-genes", type=int, default=50)
    parser.add_argument("--max-ml-genes", type=int, default=2000)
    parser.add_argument("--rf-top-n", type=int, default=50)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=8)
    return parser.parse_args()


def load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.Series]:
    x_path = input_dir / "pseudobulk_log2cpm.csv.gz"
    meta_path = input_dir / "sample_metadata.csv"
    if not x_path.exists():
        raise FileNotFoundError(f"Missing pseudobulk expression: {x_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing sample metadata: {meta_path}")

    x = pd.read_csv(x_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    y = meta["label_numeric"].astype(int)
    common = x.index.intersection(y.index)
    x = x.loc[common].sort_index()
    y = y.loc[common].sort_index()
    return x, y


def fit_final_classifier(x_train: pd.DataFrame, y_train: pd.Series) -> Pipeline:
    model = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="l2",
                    C=1.0,
                    class_weight="balanced",
                    max_iter=10000,
                    random_state=42,
                ),
            ),
        ]
    )
    model.fit(x_train, y_train)
    return model


def choose_fold_biomarkers(biomarker_table: pd.DataFrame) -> tuple[list[str], str]:
    selected = biomarker_table.loc[biomarker_table["selected_by_both_ml"], "gene"].astype(str).tolist()
    if selected:
        return selected, "lasso_and_rf_top"

    fallback = biomarker_table.loc[
        biomarker_table["selected_by_lasso"] | biomarker_table["selected_by_rf_top"],
        "gene",
    ].astype(str).head(30).tolist()
    if fallback:
        return fallback, "fallback_lasso_or_rf_top30"

    return biomarker_table["gene"].astype(str).head(30).tolist(), "fallback_combined_rank_top30"


def run_nested_cv(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_json(vars(args), out_dir / "nested_cv_parameters.json")

    x_all, y_all = load_inputs(input_dir)
    min_class_count = int(y_all.value_counts().min())
    if min_class_count < 2:
        raise ValueError("Each class needs at least two samples for nested CV.")
    n_splits = min(args.outer_splits, min_class_count)
    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.random_state)

    fold_metric_rows = []
    prediction_rows = []
    selected_rows = []

    for fold_id, (train_idx, test_idx) in enumerate(outer_cv.split(x_all, y_all), start=1):
        print(f"\n===== outer fold {fold_id}/{n_splits} =====")
        x_train_all = x_all.iloc[train_idx]
        x_test_all = x_all.iloc[test_idx]
        y_train = y_all.iloc[train_idx]
        y_test = y_all.iloc[test_idx]

        print("train label counts:", y_train.value_counts().to_dict())
        print("test label counts:", y_test.value_counts().to_dict())

        deg_train = make_deg_table(x_train_all, y_train, target_sum=1e6)
        candidate_genes, candidate_table = select_ml_genes(
            deg=deg_train,
            log2cpm=x_train_all,
            deg_padj=args.deg_padj,
            deg_log2fc=args.deg_log2fc,
            min_ml_genes=args.min_ml_genes,
            max_ml_genes=args.max_ml_genes,
        )

        inner_cv = cv_splits_for_y(y_train, args.random_state + fold_id)
        ranking, lasso_info = fit_lasso_and_rf(
            x=x_train_all.loc[:, candidate_genes],
            y=y_train,
            cv=inner_cv,
            random_state=args.random_state + fold_id,
            n_jobs=args.n_jobs,
        )
        biomarker_table = add_biomarker_annotations(ranking, deg_train, args.rf_top_n)
        selected_genes, selection_rule = choose_fold_biomarkers(biomarker_table)

        x_train = x_train_all.loc[:, selected_genes]
        x_test = x_test_all.loc[:, selected_genes]
        model = fit_final_classifier(x_train, y_train)
        proba = model.predict_proba(x_test)[:, 1]
        pred = model.predict(x_test)
        metrics = model_metrics(y_test, pred, proba)
        metrics.update(
            {
                "fold": fold_id,
                "n_train": int(len(y_train)),
                "n_test": int(len(y_test)),
                "n_candidate_genes": int(len(candidate_genes)),
                "n_selected_biomarkers": int(len(selected_genes)),
                "selection_rule": selection_rule,
                "lasso_selected_c": float(lasso_info.get("lasso_selected_c", np.nan)),
                "lasso_refit_note": str(lasso_info.get("lasso_refit_note", "")),
            }
        )
        fold_metric_rows.append(metrics)

        fold_pred = pd.DataFrame(
            {
                "fold": fold_id,
                "sample_id": x_test.index,
                "true_label_numeric": y_test.values,
                "true_label": y_test.map({0: "stable", 1: "unstable"}).values,
                "prob_unstable": proba,
                "pred_label_numeric": pred,
                "pred_label": pd.Series(pred).map({0: "stable", 1: "unstable"}).values,
                "n_selected_biomarkers": len(selected_genes),
                "selection_rule": selection_rule,
            }
        )
        prediction_rows.append(fold_pred)

        fold_selected = biomarker_table[biomarker_table["gene"].isin(selected_genes)].copy()
        fold_selected.insert(0, "fold", fold_id)
        fold_selected.insert(1, "selection_rule", selection_rule)
        fold_selected["n_train"] = len(y_train)
        fold_selected["n_test"] = len(y_test)
        selected_rows.append(fold_selected)

        candidate_table.to_csv(out_dir / f"fold_{fold_id}_candidate_genes.csv", index=False)
        biomarker_table.to_csv(out_dir / f"fold_{fold_id}_train_biomarker_ranking.csv", index=False)
        print(f"selected {len(selected_genes)} biomarkers by {selection_rule}")

    fold_metrics = pd.DataFrame(fold_metric_rows)
    predictions = pd.concat(prediction_rows, axis=0, ignore_index=True)
    selected = pd.concat(selected_rows, axis=0, ignore_index=True)

    overall = model_metrics(
        predictions["true_label_numeric"].astype(int),
        predictions["pred_label_numeric"].astype(int),
        predictions["prob_unstable"].astype(float),
    )
    overall_df = pd.DataFrame([overall])
    overall_df.insert(0, "model", "nested_fold_internal_selection_logistic")
    overall_df["n_outer_splits"] = n_splits
    overall_df["note"] = "Each fold performs DEG, LASSO, RF, biomarker selection, and final model fitting on training samples only."

    frequency = (
        selected.groupby("gene", as_index=False)
        .agg(
            n_folds_selected=("fold", "nunique"),
            mean_lasso_coef=("lasso_coef", "mean"),
            mean_abs_lasso_coef=("abs_lasso_coef", "mean"),
            mean_rf_importance=("rf_importance", "mean"),
            mean_rf_rank=("rf_rank", "mean"),
            mean_log2fc=("log2FC_mean_cpm", "mean"),
            median_padj=("padj", "median"),
        )
        .sort_values(["n_folds_selected", "mean_rf_importance", "mean_abs_lasso_coef"], ascending=[False, False, False])
    )

    fold_metrics.to_csv(out_dir / "nested_cv_fold_metrics.csv", index=False)
    predictions.to_csv(out_dir / "nested_cv_predictions.csv", index=False)
    selected.to_csv(out_dir / "nested_cv_selected_biomarkers_by_fold.csv", index=False)
    frequency.to_csv(out_dir / "nested_cv_selected_biomarker_frequency.csv", index=False)
    overall_df.to_csv(out_dir / "nested_cv_overall_metrics.csv", index=False)
    write_summary(out_dir, y_all, overall_df, fold_metrics, frequency)
    plot_nested_results(out_dir, predictions, frequency)

    return overall_df, fold_metrics, predictions, frequency


def write_summary(
    out_dir: Path,
    y_all: pd.Series,
    overall_df: pd.DataFrame,
    fold_metrics: pd.DataFrame,
    frequency: pd.DataFrame,
) -> None:
    lines = [
        "# Strict nested biomarker CV",
        "",
        "Validation design:",
        "- Outer folds split sample IDs, not cells.",
        "- In each outer fold, DEG, LASSO, RF and biomarker selection use training samples only.",
        "- Held-out samples are used only for final prediction and metrics.",
        "",
        "Label counts:",
        y_all.map({0: "stable", 1: "unstable"}).value_counts().to_string(),
        "",
        "Overall held-out metrics:",
        overall_df.to_string(index=False),
        "",
        "Fold metrics:",
        fold_metrics.to_string(index=False),
        "",
        "Most repeatedly selected biomarkers:",
        frequency.head(30).to_string(index=False),
        "",
    ]
    (out_dir / "README_nested_cv.md").write_text("\n".join(lines))


def plot_nested_results(out_dir: Path, predictions: pd.DataFrame, frequency: pd.DataFrame) -> None:
    sns.set_theme(style="ticks")

    cm = confusion_matrix(
        predictions["true_label_numeric"].astype(int),
        predictions["pred_label_numeric"].astype(int),
        labels=[0, 1],
    )
    fig, ax = plt.subplots(figsize=(2.6, 2.4), constrained_layout=True)
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar=False,
        xticklabels=["stable", "unstable"],
        yticklabels=["stable", "unstable"],
        ax=ax,
    )
    ax.set_xlabel("Predicted", fontsize=8, fontweight="bold")
    ax.set_ylabel("True", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    save_figure(fig, out_dir / "nested_cv_confusion_matrix")

    plot_df = predictions.copy()
    plot_df["true_label"] = pd.Categorical(plot_df["true_label"], categories=["stable", "unstable"], ordered=True)
    fig, ax = plt.subplots(figsize=(3.4, 2.7), constrained_layout=True)
    sns.stripplot(data=plot_df, x="true_label", y="prob_unstable", hue="pred_label", dodge=True, size=5, ax=ax)
    ax.axhline(0.5, color="#333333", lw=0.6, ls="--")
    ax.set_xlabel("True label", fontsize=8, fontweight="bold")
    ax.set_ylabel("Held-out probability unstable", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    ax.legend(title="Predicted", fontsize=6, title_fontsize=7, frameon=False)
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "nested_cv_prob_unstable")

    top = frequency.head(30).iloc[::-1]
    if not top.empty:
        fig, ax = plt.subplots(figsize=(4.0, max(2.6, 0.16 * len(top))), constrained_layout=True)
        ax.barh(top["gene"], top["n_folds_selected"], color="#4B7A9F")
        ax.set_xlabel("Folds selected", fontsize=8, fontweight="bold")
        ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=7)
        sns.despine(ax=ax)
        save_figure(fig, out_dir / "nested_cv_biomarker_selection_frequency")


def main() -> None:
    args = parse_args()
    overall_df, fold_metrics, _, frequency = run_nested_cv(args)
    print("\nDone.")
    print(f"Output directory: {Path(args.out_dir).resolve()}")
    print("\nOverall held-out metrics:")
    print(overall_df.to_string(index=False))
    print("\nFold metrics:")
    print(fold_metrics.to_string(index=False))
    print("\nMost repeatedly selected biomarkers:")
    print(frequency.head(30).to_string(index=False))


if __name__ == "__main__":
    main()
