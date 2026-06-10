#!/usr/bin/env python
"""Validate biological states represented by selected ML marker genes.

The validation mirrors the article's logic at sample level:
1) stable/unstable expression direction,
2) single-marker ROC,
3) marker-panel state score,
4) correlation with immune cell infiltration,
5) correlation with article-inspired plaque biology genes,
6) cross-sample heatmap and a Markdown report.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


MARKER_GENES = ["ZNF555", "AMFR", "FOXK1", "WDR91", "IL12RB2", "ICAM3"]
MARKER_DIRECTIONS = {
    "ZNF555": "stable_up",
    "AMFR": "stable_up",
    "FOXK1": "stable_up",
    "WDR91": "stable_up",
    "IL12RB2": "unstable_up",
    "ICAM3": "unstable_up",
}

STATE_AXIS = {
    "immune_activation_adhesion": ["IL12RB2", "ICAM3"],
    "stable_like_regulatory_homeostasis": ["ZNF555", "AMFR", "FOXK1", "WDR91"],
}

ARTICLE_BIOLOGY_GENES = {
    "macrophage_m1_instability": ["CD68", "CD80", "CD86", "IL1B", "TNF", "NLRP3", "CCL2", "CXCL8"],
    "antigen_inflammation_lipid": ["CD74", "B2M", "HLA-DRA", "MMP9", "CTSL", "IFI30", "CD36", "APOE"],
    "smc_contractile_stability": ["TAGLN", "ACAT2", "MYH10", "MYH11"],
    "cd8_t_nk_protective": ["CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "IFNG", "KLRD1"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="output_article_biomarker_pipeline")
    parser.add_argument("--out-dir", default="output_ml_marker_biological_validation")
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


def load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(input_dir / "pseudobulk_log2cpm.csv.gz", index_col=0)
    meta = pd.read_csv(input_dir / "sample_metadata.csv", index_col=0)
    prop = pd.read_csv(input_dir / "immune_celltype_proportions__cell_type_level1_corrected.csv", index_col=0)
    common = expr.index.intersection(meta.index).intersection(prop.index)
    return expr.loc[common].sort_index(), meta.loc[common].sort_index(), prop.loc[common].sort_index()


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    return (df - df.mean(axis=0)) / df.std(axis=0, ddof=0).replace(0, np.nan)


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def marker_expression_stats(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    y = meta["label_numeric"].astype(int)
    rows = []
    for gene in MARKER_GENES:
        stable = expr.loc[y == 0, gene]
        unstable = expr.loc[y == 1, gene]
        stat, pvalue = stats.ttest_ind(unstable, stable, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "gene": gene,
                "expected_direction": MARKER_DIRECTIONS[gene],
                "mean_stable": stable.mean(),
                "mean_unstable": unstable.mean(),
                "unstable_minus_stable": unstable.mean() - stable.mean(),
                "t_stat": stat,
                "pvalue": pvalue,
                "direction_in_data": "unstable_up" if unstable.mean() >= stable.mean() else "stable_up",
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = bh_adjust(out["pvalue"].to_numpy())
    out["direction_matches_expected"] = out["expected_direction"] == out["direction_in_data"]
    return out.sort_values(["direction_matches_expected", "padj"], ascending=[False, True])


def marker_roc(expr: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    y = meta["label_numeric"].astype(int)
    rows = []
    fig, ax = plt.subplots(figsize=(3.4, 3.2), constrained_layout=True)
    for gene in MARKER_GENES:
        score = expr[gene].astype(float)
        if MARKER_DIRECTIONS[gene] == "stable_up":
            score = -score
        auc = roc_auc_score(y, score)
        fpr, tpr, _ = roc_curve(y, score)
        rows.append({"gene": gene, "expected_direction": MARKER_DIRECTIONS[gene], "auc_for_unstable": auc})
        ax.plot(fpr, tpr, lw=1.2, label=f"{gene} ({auc:.2f})")
    ax.plot([0, 1], [0, 1], color="#333333", lw=0.6, ls="--")
    ax.set_xlabel("False positive rate", fontsize=8, fontweight="bold")
    ax.set_ylabel("True positive rate", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    ax.legend(frameon=False, fontsize=5.5, title="Marker AUC", title_fontsize=6)
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "single_marker_roc")
    return pd.DataFrame(rows).sort_values("auc_for_unstable", ascending=False)


def panel_scores(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    marker_z = zscore(expr[MARKER_GENES]).fillna(0.0)
    scores = pd.DataFrame(index=expr.index)
    scores["immune_activation_adhesion_score"] = marker_z[STATE_AXIS["immune_activation_adhesion"]].mean(axis=1)
    scores["stable_like_homeostasis_score"] = marker_z[STATE_AXIS["stable_like_regulatory_homeostasis"]].mean(axis=1)
    scores["directional_unstable_panel_score"] = (
        marker_z["IL12RB2"]
        + marker_z["ICAM3"]
        - marker_z["ZNF555"]
        - marker_z["AMFR"]
        - marker_z["FOXK1"]
        - marker_z["WDR91"]
    ) / len(MARKER_GENES)
    scores["label"] = meta["label"]
    scores["label_numeric"] = meta["label_numeric"].astype(int)
    return scores


def score_stats(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in ["immune_activation_adhesion_score", "stable_like_homeostasis_score", "directional_unstable_panel_score"]:
        stable = scores.loc[scores["label_numeric"] == 0, col]
        unstable = scores.loc[scores["label_numeric"] == 1, col]
        stat, pvalue = stats.ttest_ind(unstable, stable, equal_var=False, nan_policy="omit")
        auc_score = scores[col] if col != "stable_like_homeostasis_score" else -scores[col]
        rows.append(
            {
                "score": col,
                "mean_stable": stable.mean(),
                "mean_unstable": unstable.mean(),
                "unstable_minus_stable": unstable.mean() - stable.mean(),
                "pvalue": pvalue,
                "auc_for_unstable": roc_auc_score(scores["label_numeric"].astype(int), auc_score),
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = bh_adjust(out["pvalue"].to_numpy())
    return out


def immune_correlations(expr: pd.DataFrame, prop: pd.DataFrame) -> pd.DataFrame:
    marker_z = zscore(expr[MARKER_GENES]).fillna(0.0)
    rows = []
    for gene in MARKER_GENES:
        for celltype_feature in prop.columns:
            rho, pvalue = stats.spearmanr(marker_z[gene], prop[celltype_feature], nan_policy="omit")
            rows.append(
                {
                    "gene": gene,
                    "celltype_feature": celltype_feature.replace("cell_type_level1_corrected__", ""),
                    "spearman_r": rho,
                    "pvalue": pvalue,
                }
            )
    out = pd.DataFrame(rows)
    out["padj"] = bh_adjust(out["pvalue"].to_numpy())
    return out.sort_values(["gene", "padj", "spearman_r"], ascending=[True, True, False])


def article_axis_scores(expr: pd.DataFrame) -> pd.DataFrame:
    z = zscore(expr).fillna(0.0)
    out = pd.DataFrame(index=expr.index)
    for axis, genes in ARTICLE_BIOLOGY_GENES.items():
        present = [gene for gene in genes if gene in z.columns]
        if present:
            out[axis] = z[present].mean(axis=1)
    return out


def article_axis_correlations(scores: pd.DataFrame, article_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    score_cols = ["immune_activation_adhesion_score", "stable_like_homeostasis_score", "directional_unstable_panel_score"]
    for score_col in score_cols:
        for axis in article_scores.columns:
            rho, pvalue = stats.spearmanr(scores[score_col], article_scores[axis], nan_policy="omit")
            rows.append({"marker_score": score_col, "article_axis": axis, "spearman_r": rho, "pvalue": pvalue})
    out = pd.DataFrame(rows)
    out["padj"] = bh_adjust(out["pvalue"].to_numpy())
    return out.sort_values(["marker_score", "padj"], ascending=[True, True])


def plot_marker_expression(expr: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> None:
    plot_df = expr[MARKER_GENES].copy()
    plot_df["label"] = meta["label"]
    long_df = plot_df.melt(id_vars="label", var_name="gene", value_name="log2cpm")
    fig, ax = plt.subplots(figsize=(5.2, 3.0), constrained_layout=True)
    sns.boxplot(data=long_df, x="gene", y="log2cpm", hue="label", width=0.65, fliersize=0, ax=ax)
    sns.stripplot(data=long_df, x="gene", y="log2cpm", hue="label", dodge=True, size=3, alpha=0.75, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, fontsize=7, title="Plaque type", title_fontsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("Pseudobulk log2(CPM+1)", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "marker_expression_by_plaque_type")


def plot_panel_scores(scores: pd.DataFrame, out_dir: Path) -> None:
    score_df = scores.reset_index().rename(columns={scores.index.name or "index": "sample_id"})
    long_df = score_df.melt(
        id_vars=["sample_id", "label", "label_numeric"],
        value_vars=["immune_activation_adhesion_score", "stable_like_homeostasis_score", "directional_unstable_panel_score"],
        var_name="score",
        value_name="value",
    )
    fig, ax = plt.subplots(figsize=(5.0, 3.0), constrained_layout=True)
    sns.boxplot(data=long_df, x="score", y="value", hue="label", fliersize=0, width=0.65, ax=ax)
    sns.stripplot(data=long_df, x="score", y="value", hue="label", dodge=True, size=3, alpha=0.75, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, fontsize=7, title="Plaque type", title_fontsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("Z-score derived state score", fontsize=8, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    sns.despine(ax=ax)
    save_figure(fig, out_dir / "marker_panel_scores_by_plaque_type")


def plot_heatmap(expr: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> None:
    sample_order = meta.sort_values(["label_numeric", "tissue", "sample" if "sample" in meta.columns else "label"]).index
    z = zscore(expr.loc[sample_order, MARKER_GENES]).fillna(0.0).T
    fig, ax = plt.subplots(figsize=(4.8, 2.6), constrained_layout=True)
    sns.heatmap(z, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "Z-score"}, xticklabels=False, yticklabels=True)
    ax.set_xlabel("Samples ordered by Plaque_type", fontsize=8, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=7)
    save_figure(fig, out_dir / "marker_zscore_heatmap")


def plot_immune_correlation(corr: pd.DataFrame, out_dir: Path) -> None:
    matrix = corr.pivot(index="gene", columns="celltype_feature", values="spearman_r").reindex(MARKER_GENES)
    fig, ax = plt.subplots(figsize=(6.0, 2.8), constrained_layout=True)
    sns.heatmap(matrix, cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={"label": "Spearman r"})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save_figure(fig, out_dir / "marker_immune_celltype_correlation_heatmap")


def write_report(
    out_dir: Path,
    expression_stats: pd.DataFrame,
    roc: pd.DataFrame,
    score_summary: pd.DataFrame,
    immune_corr: pd.DataFrame,
    article_corr: pd.DataFrame,
) -> None:
    top_immune = immune_corr.reindex(immune_corr["spearman_r"].abs().sort_values(ascending=False).index).head(20)

    def as_markdown_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        text_df = df.copy()
        for col in text_df.columns:
            if pd.api.types.is_float_dtype(text_df[col]):
                text_df[col] = text_df[col].map(lambda value: f"{value:.4g}" if pd.notna(value) else "")
            else:
                text_df[col] = text_df[col].astype(str)
        header = "| " + " | ".join(text_df.columns) + " |"
        sep = "| " + " | ".join(["---"] * len(text_df.columns)) + " |"
        rows = ["| " + " | ".join(map(str, row)) + " |" for row in text_df.to_numpy()]
        return "\n".join([header, sep] + rows)

    lines = [
        "# ML Marker Biological State Validation",
        "",
        "## 验证 1：Marker 来源和状态假设",
        "",
        "本文件验证 6 个预测型 marker 对应的生物学状态：`ZNF555`, `AMFR`, `FOXK1`, `WDR91`, `IL12RB2`, `ICAM3`。",
        "",
        "- Stable-like/homeostasis markers: `ZNF555`, `AMFR`, `FOXK1`, `WDR91`",
        "- Immune activation/adhesion markers: `IL12RB2`, `ICAM3`",
        "",
        "这些基因来自当前数据的 sample-level pseudobulk + DEG + LASSO + Random Forest + nested-CV 稳定性筛选。这里验证的是它们是否构成可解释的 unstable-plaque-associated signature，而不是证明因果调控关系。",
        "",
        "## 验证 2：Stable vs Unstable 表达方向",
        "",
        "输出表：[`marker_expression_stats.csv`](marker_expression_stats.csv)",
        "",
        as_markdown_table(expression_stats),
        "",
        "图：[`marker_expression_by_plaque_type.png`](marker_expression_by_plaque_type.png)",
        "",
        "## 验证 3：单基因诊断能力 ROC",
        "",
        "输出表：[`single_marker_roc_auc.csv`](single_marker_roc_auc.csv)",
        "",
        as_markdown_table(roc),
        "",
        "图：[`single_marker_roc.png`](single_marker_roc.png)",
        "",
        "## 验证 4：组合 marker panel 状态分数",
        "",
        "这里构建三个 score：",
        "",
        "- `immune_activation_adhesion_score = mean(z(IL12RB2), z(ICAM3))`",
        "- `stable_like_homeostasis_score = mean(z(ZNF555), z(AMFR), z(FOXK1), z(WDR91))`",
        "- `directional_unstable_panel_score = mean(z(IL12RB2), z(ICAM3), -z(ZNF555), -z(AMFR), -z(FOXK1), -z(WDR91))`",
        "",
        "输出表：[`marker_panel_score_stats.csv`](marker_panel_score_stats.csv)",
        "",
        as_markdown_table(score_summary),
        "",
        "图：[`marker_panel_scores_by_plaque_type.png`](marker_panel_scores_by_plaque_type.png)，[`marker_zscore_heatmap.png`](marker_zscore_heatmap.png)",
        "",
        "## 验证 5：与免疫细胞浸润的相关性",
        "",
        "借鉴原文 CIBERSORT 免疫浸润分析，这里用单细胞注释得到的 sample-level cell type proportions，与 marker 表达做 Spearman 相关。",
        "",
        "输出表：[`marker_immune_celltype_correlations.csv`](marker_immune_celltype_correlations.csv)",
        "",
        "Top absolute correlations:",
        "",
        as_markdown_table(top_immune),
        "",
        "图：[`marker_immune_celltype_correlation_heatmap.png`](marker_immune_celltype_correlation_heatmap.png)",
        "",
        "## 验证 6：与论文启发的 plaque biology axes 的相关性",
        "",
        "这里计算 marker panel score 与论文中提到的 macrophage/M1、antigen/inflammation/lipid、SMC contractile、CD8/NK protective gene-axis score 的相关性。",
        "",
        "输出表：[`marker_score_article_axis_correlations.csv`](marker_score_article_axis_correlations.csv)",
        "",
        as_markdown_table(article_corr),
        "",
        "## 解释结论",
        "",
        "如果结果方向稳定，可以将这 6 个 marker 解释为一个预测型状态 signature：`IL12RB2/ICAM3` 代表免疫激活和细胞黏附/互作增强，`ZNF555/AMFR/FOXK1/WDR91` 代表 stable-like 转录调控、蛋白/脂质稳态和内吞运输程序。需要注意，这仍是当前数据内的候选解释，最终需要外部队列或实验验证。",
        "",
    ]
    (out_dir / "ML_marker_biological_state_validation.md").write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_plot_style()

    expr, meta, prop = load_inputs(input_dir)
    missing = [gene for gene in MARKER_GENES if gene not in expr.columns]
    if missing:
        raise ValueError(f"Missing marker genes in pseudobulk expression: {missing}")

    expression_stats = marker_expression_stats(expr, meta)
    roc = marker_roc(expr, meta, out_dir)
    scores = panel_scores(expr, meta)
    score_summary = score_stats(scores)
    immune_corr = immune_correlations(expr, prop)
    article_scores = article_axis_scores(expr)
    article_corr = article_axis_correlations(scores, article_scores)

    expression_stats.to_csv(out_dir / "marker_expression_stats.csv", index=False)
    roc.to_csv(out_dir / "single_marker_roc_auc.csv", index=False)
    scores.to_csv(out_dir / "marker_panel_scores_by_sample.csv")
    score_summary.to_csv(out_dir / "marker_panel_score_stats.csv", index=False)
    immune_corr.to_csv(out_dir / "marker_immune_celltype_correlations.csv", index=False)
    article_scores.to_csv(out_dir / "article_biology_axis_scores_by_sample.csv")
    article_corr.to_csv(out_dir / "marker_score_article_axis_correlations.csv", index=False)

    plot_marker_expression(expr, meta, out_dir)
    plot_panel_scores(scores, out_dir)
    plot_heatmap(expr, meta, out_dir)
    plot_immune_correlation(immune_corr, out_dir)
    write_report(out_dir, expression_stats, roc, score_summary, immune_corr, article_corr)

    print("Done.")
    print(f"Output directory: {out_dir.resolve()}")
    print(f"Report: {(out_dir / 'ML_marker_biological_state_validation.md').resolve()}")


if __name__ == "__main__":
    main()
