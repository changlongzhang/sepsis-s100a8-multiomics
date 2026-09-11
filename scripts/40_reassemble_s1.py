from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(r"D:\桌面\sepsis")
DATA = ROOT / "review_revision" / "09_figures"
OUT = ROOT / "文章" / "返修投稿新主图"
PREVIEW = DATA / "S1_preview.png"

NAVY, TEAL, CORAL = "#315B78", "#5B958C", "#D56A54"
INK, MID, LIGHT, PALE = "#20252A", "#7A858D", "#D9DEE2", "#F2F4F5"
SHAP_CMAP = LinearSegmentedColormap.from_list("shap", ["#3569A8", "#F3F0E8", "#CB4C64"])


def setup_style():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": .7,
                         "xtick.major.width": .6, "ytick.major.width": .6,
                         "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.facecolor": "white"})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(colors=INK, pad=2)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def panel_label(ax, text, x=-.13, y=1.08):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="left", color=INK, clip_on=False)


def load_data():
    return {"cv": pd.read_csv(DATA / "S1_lasso_cv.csv"),
            "path": pd.read_csv(DATA / "S1_lasso_path.csv"),
            "rank": pd.read_csv(DATA / "S1_model_rankings.csv"),
            "shap": pd.read_csv(DATA / "S1_shap_values.csv"),
            "x": pd.read_csv(DATA / "S1_feature_values.csv"),
            "base": float(pd.read_csv(DATA / "S1_shap_baseline.csv").iloc[0, 0])}


def lambda_positions(d):
    imin = d.cvm.idxmin()
    lam_min = d.loc[imin, "ll"]
    threshold = d.loc[imin, "cvm"] + d.loc[imin, "cvsd"]
    return lam_min, d[d.cvm <= threshold].ll.max()


def plot_cv(ax, d):
    d = d.sort_values("ll")
    ax.vlines(d.ll, d.cvlo, d.cvup, color=LIGHT, lw=.7, zorder=1)
    ax.scatter(d.ll, d.cvm, s=9, facecolor=NAVY, edgecolor="white", linewidth=.25, zorder=2)
    lam_min, lam_1se = lambda_positions(d)
    ax.axvline(lam_min, color=CORAL, lw=1, ls=(0, (3, 2)))
    ax.axvline(lam_1se, color=TEAL, lw=1, ls=(0, (3, 2)))
    ymin, ymax = ax.get_ylim()
    ax.text(lam_min, ymin+.02*(ymax-ymin), r"$\lambda_{min}$", color=CORAL,
            ha="right", va="bottom", fontsize=8)
    ax.text(lam_1se, ymin+.02*(ymax-ymin), r"$\lambda_{1se}$", color=TEAL,
            ha="left", va="bottom", fontsize=8)
    ax.set(xlabel=r"log($\lambda$)", ylabel="CV deviance")
    clean(ax); panel_label(ax, "A", x=-.14, y=1.10)


def plot_paths(ax, d, cv):
    for gene, g in d.groupby("coef"):
        color, lw, alpha, z = LIGHT, .65, .85, 1
        if gene == "S100A8": color, lw, alpha, z = NAVY, 1.5, 1, 3
        elif gene == "MAFG": color, lw, alpha, z = CORAL, 1.2, 1, 2
        ax.plot(np.log(g["lambda"]), g["value"], color=color, lw=lw, alpha=alpha, zorder=z)
    lam_min, lam_1se = lambda_positions(cv)
    ax.axvline(lam_min, color=CORAL, lw=1, ls=(0, (3, 2)))
    ax.axvline(lam_1se, color=TEAL, lw=1, ls=(0, (3, 2)))
    ax.plot([], [], color=NAVY, lw=1.5, label=r"$S100A8$")
    ax.plot([], [], color=CORAL, lw=1.2, label=r"$MAFG$")
    ax.legend(frameon=False, fontsize=8, loc="upper right", handlelength=1.6)
    ax.set(xlabel=r"log($\lambda$)", ylabel="Regularized coefficient")
    clean(ax); panel_label(ax, "B")


def plot_rank_matrix(ax, d):
    models = ["LASSO", "glmBoost", "XGBoost", "Ranger"]
    labels = ["LASSO", "glmBoost", "XGBoost", "Random forest"]
    d = d.copy()
    d["Rank"] = d.groupby("Model").Value.rank(method="first", ascending=False).astype(int)
    genes = d.groupby("Gene").Rank.agg(["min", "mean"]).sort_values(["min", "mean"]).index.tolist()
    for i, gene in enumerate(genes):
        for j, model in enumerate(models):
            hit = d[(d.Gene == gene) & (d.Model == model)]
            if hit.empty:
                color, rank = PALE, None
            else:
                rank = int(hit.Rank.iloc[0]); t = (5-rank)/4
                rgb = np.array(mpl.colors.to_rgb(TEAL))
                color = tuple(1-(1-rgb)*(.32+.68*t))
            ax.add_patch(plt.Rectangle((j-.42, i-.42), .84, .84, fc=color, ec="white", lw=.7))
            if rank:
                ax.text(j, i, str(rank), ha="center", va="center", fontsize=8,
                        color="white" if rank <= 2 else INK, fontweight="bold")
    ax.set(xlim=(-.5, 3.5), ylim=(len(genes)-.5, -.5))
    ax.set_xticks(range(4), labels); ax.set_yticks(range(len(genes)), genes, fontstyle="italic")
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel("Within-model rank (1 = highest)")
    panel_label(ax, "C")


def plot_beeswarm(ax, shap, x):
    order = shap.abs().mean().sort_values(ascending=False).index.tolist()
    rng = np.random.default_rng(123)
    for yi, gene in enumerate(order):
        vals, feat = shap[gene].to_numpy(), x[gene].to_numpy()
        z = (feat-np.nanmean(feat))/(np.nanstd(feat) or 1)
        ax.scatter(vals, yi+rng.normal(0, .085, len(vals)), c=np.clip(z, -2, 2),
                   cmap=SHAP_CMAP, vmin=-2, vmax=2, s=10, alpha=.86, edgecolors="none")
    ax.axvline(0, color=MID, lw=.7)
    ax.set_yticks(range(len(order)), order, fontstyle="italic")
    ax.set_ylim(len(order)-.6, -.6); ax.set(xlabel="SHAP value", ylabel="Feature")
    clean(ax); panel_label(ax, "D")
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-2, 2), cmap=SHAP_CMAP)
    cb = plt.colorbar(sm, ax=ax, fraction=.045, pad=.03)
    cb.set_ticks([-2, 2], labels=["Low", "High"]); cb.ax.tick_params(length=0, labelsize=8)
    cb.set_label("Feature value", fontsize=8); cb.outline.set_visible(False)


def plot_waterfall(ax, shap, x, baseline, row=1):
    vals, feat = shap.iloc[row], x.iloc[row]
    order = vals.abs().sort_values(ascending=True).index.tolist()
    current = baseline
    for i, gene in enumerate(order):
        v = vals[gene]
        ax.barh(i, v, left=current, height=.62, color=CORAL if v >= 0 else NAVY,
                edgecolor="white", linewidth=.5)
        ax.text(current+v/2, i, f"{v:+.3f}", ha="center", va="center", fontsize=8,
                color="white" if abs(v) > .009 else INK)
        current += v
    ax.axvline(baseline, color=MID, lw=.8, ls=(0, (2, 2)))
    ax.axvline(current, color=INK, lw=.8, ls=(0, (2, 2)))
    ax.set_yticks(range(len(order)), [f"{g}  ({feat[g]:.2f})" for g in order], fontstyle="italic")
    ax.set_xlabel("Predicted probability contribution")
    ax.text(baseline, len(order)-.18, f"Baseline = {baseline:.3f}", ha="left", va="bottom", fontsize=8, color=MID)
    ax.text(current, len(order)-.18, f"Prediction = {current:.3f}", ha="right", va="bottom", fontsize=8, color=INK)
    clean(ax); panel_label(ax, "E", x=-.07, y=1.17)


def plot_dependence(axes, shap, x):
    genes = ["C5orf32", "S100A8", "MAFG", "CD177", "MYL6", "EXOSC4"]
    partners = {"C5orf32":"S100A8", "S100A8":"C5orf32", "MAFG":"S100A8",
                "CD177":"MYL6", "MYL6":"S100A8", "EXOSC4":"S100A8"}
    for ax, gene in zip(axes, genes):
        partner = partners[gene]
        ax.scatter(x[gene], shap[gene], c=np.clip(x[partner], -2, 2), cmap=SHAP_CMAP,
                   vmin=-2, vmax=2, s=10, alpha=.88, edgecolors="none")
        ax.axhline(0, color=LIGHT, lw=.65)
        ax.set_xlabel(gene, fontstyle="italic"); ax.set_ylabel("SHAP value")
        ax.text(.98, .08, f"Color: {partner}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color=MID, bbox=dict(facecolor="white", edgecolor="none", alpha=.72, pad=.5))
        clean(ax)
    panel_label(axes[0], "F", x=-.28, y=1.20)


def main():
    setup_style(); d = load_data(); OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.5, 8.75))
    outer = fig.add_gridspec(4, 1, height_ratios=[1.9, 1.95, 1.4, 2.4],
                             left=.11, right=.975, bottom=.065, top=.975, hspace=.58)
    top = outer[0].subgridspec(1, 2, wspace=.34)
    mid = outer[1].subgridspec(1, 2, wspace=.45)
    ax_a, ax_b = fig.add_subplot(top[0,0]), fig.add_subplot(top[0,1])
    ax_c, ax_d = fig.add_subplot(mid[0,0]), fig.add_subplot(mid[0,1])
    ax_e = fig.add_subplot(outer[2])
    dep = outer[3].subgridspec(2, 3, wspace=.42, hspace=.62)
    dep_axes = [fig.add_subplot(dep[i,j]) for i in range(2) for j in range(3)]
    plot_cv(ax_a, d["cv"]); plot_paths(ax_b, d["path"], d["cv"])
    plot_rank_matrix(ax_c, d["rank"]); plot_beeswarm(ax_d, d["shap"], d["x"])
    plot_waterfall(ax_e, d["shap"], d["x"], d["base"]); plot_dependence(dep_axes, d["shap"], d["x"])
    fig.savefig(OUT / "S1_fig_python_master.pdf", format="pdf")
    fig.savefig(PREVIEW, dpi=300)
    plt.close(fig)


if __name__ == "__main__": main()
