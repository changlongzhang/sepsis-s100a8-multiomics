"""Rebuild Fig. 2 source panels as editable, publication-grade Python graphics.

Each panel is exported as an independent figure before final composition.  The
script deliberately avoids rasterizing PDF/SVG source panels.  Run with
``--panel A`` through ``--panel F`` or ``--panel all``.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


PROJECT = Path(r"D:\桌面\sepsis")
REVISION = PROJECT / "review_revision"
FIG_DIR = REVISION / "09_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
FINAL_DIR = PROJECT / "文章" / "返修投稿新主图"
FINAL_DIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 9,
        "text.color": "#26323C",
        "axes.labelcolor": "#26323C",
        "xtick.color": "#26323C",
        "ytick.color": "#26323C",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
    }
)

INK = "#26323C"
MUTED = "#65727E"
PALE = "#EEF2F5"
CONTROL = "#4C78A8"
DISEASE = "#D66A5E"
TPL = "#4F938D"
PRECURSOR = "#D49A45"
TARGET = "#C97832"
COLORS = [CONTROL, "#7897B5", TPL, PRECURSOR, "#7A7180"]
FILLS = ["#E8EDF1", "#EDF1F4", "#E2EEF2", "#F5EBDD", "#EEE9F2"]


def export(fig: plt.Figure, stem: str, dpi: int = 600) -> None:
    """Overwrite one independent panel in editable and submission formats."""
    target_source_size = {
        "Fig_ClassBalance_StudyDesign": 9.0,
        "Fig_ModelPerformance_BalanceStrategies": 11.0,
        "Fig_S100A8_SelectionFrequency": 12.0,
        "Fig_S100A8_RankDistribution": 12.0,
        "Fig_AUC_Distribution_Downsampling": 12.0,
        "Fig_Bootstrap_EffectSize": 12.0,
    }[stem]
    for item in fig.findobj(match=mpl.text.Text):
        if item.get_text().strip() and item.get_fontsize() < target_source_size:
            item.set_fontsize(target_source_size)
    try:
        fig.tight_layout(pad=0.35)
    except Exception:
        pass
    base = FIG_DIR / stem
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=dpi,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(base.with_suffix(".png"), dpi=300)
    plt.close(fig)


def rounded_box(ax, x, y, w, h, face, edge="none", radius=0.02, lw=0.8, z=1):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
        transform=ax.transAxes,
        clip_on=False,
        zorder=z,
    )
    ax.add_patch(patch)
    return patch


def panel_a() -> None:
    """Minimal audit trail used in top ML papers: one rule, five explicit stages."""
    fig, ax = plt.subplots(figsize=(7.20, 1.88))
    fig.patch.set_facecolor("white")
    ax.set_axis_off()

    titles = [
        "Discovery\ncohort",
        "Fold-internal\npre-processing",
        "Repeated nested\nvalidation",
        "Class-balance\nstress test",
        "Uncertainty\naudit",
    ]
    details = [
        "760 sepsis\n42 normal",
        "screening + scaling\ninside every\ntraining fold",
        "weighted 5-fold CV\n20 repeated cycles",
        "1:1 downsampling\n1,000 iterations",
        "stratified bootstrap\n1,000 iterations",
    ]
    xs = np.linspace(0.075, 0.925, 5)
    line_y = 0.61
    ax.annotate("", xy=(0.965, line_y), xytext=(0.035, line_y),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#9AA6AF",
                                mutation_scale=9))

    for i, (x, title, detail) in enumerate(zip(xs, titles, details), start=1):
        face = CONTROL if i in (1, 5) else "white"
        edge = CONTROL if i in (1, 5) else "#7F8B94"
        text_color = "white" if i in (1, 5) else INK
        ax.scatter(x, line_y, s=430, facecolor=face, edgecolor=edge,
                   linewidth=1.0, transform=ax.transAxes, zorder=3)
        ax.text(x, line_y, str(i), transform=ax.transAxes, ha="center",
                va="center", fontsize=7.4, fontweight="bold", color=text_color,
                zorder=4)
        ax.text(x, 0.88, title, transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, fontweight="bold", color=INK, linespacing=1.02)
        ax.text(x, 0.37, detail, transform=ax.transAxes, ha="center", va="top",
                fontsize=6.1, color=MUTED, linespacing=1.10)

    # The only pictorial element is the imbalance itself: an honest proportion bar.
    sepsis_frac = 760 / 802
    bar_x, bar_y, bar_w, bar_h = 0.020, 0.075, 0.165, 0.032
    ax.add_patch(plt.Rectangle((bar_x, bar_y), bar_w, bar_h,
                               transform=ax.transAxes, color="#DCE2E6", lw=0))
    ax.add_patch(plt.Rectangle((bar_x, bar_y), bar_w * sepsis_frac, bar_h,
                               transform=ax.transAxes, color=CONTROL, lw=0))
    ax.text(bar_x, bar_y - 0.025, "760 sepsis", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.9, color=INK)
    ax.text(bar_x + bar_w + 0.010, bar_y - 0.025, "42 normal  (18.1:1)",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.9, color=MUTED)
    ax.text(0.98, 0.075, "fold-safe preprocessing  •  nested CV  •  audit",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.2,
            color=MUTED)

    export(fig, "Fig_ClassBalance_StudyDesign")


def panel_b() -> None:
    """Benchmark-style annotated heatmap of model performance."""
    path = FIG_DIR / "Fig_ModelPerformance_BalanceStrategies_source_data.csv"
    d = pd.read_csv(path)
    # Put all metrics on a favorable 0-1 direction without hiding the transform.
    is_brier = d["metric"].eq("Brier score")
    old_lo = d.loc[is_brier, "ci_low"].copy()
    old_hi = d.loc[is_brier, "ci_high"].copy()
    d.loc[is_brier, "estimate"] = 1 - d.loc[is_brier, "estimate"]
    d.loc[is_brier, "ci_low"] = 1 - old_hi
    d.loc[is_brier, "ci_high"] = 1 - old_lo
    d.loc[is_brier, "metric"] = "1 − Brier"

    strategy_order = ["weighted", "unweighted", "downsample"]
    strategy_labels = {
        "weighted": "Weighted CV",
        "unweighted": "Unweighted CV",
        "downsample": "Fold downsampling",
    }
    algorithm_order = ["glmBoost", "LASSO", "Ranger", "XGBoost"]
    metric_order = ["ROC AUC", "PR AUC", "Balanced accuracy", "MCC", "1 − Brier"]
    metric_labels = ["ROC AUC", "PR AUC", "Balanced acc.", "MCC", "1 − Brier"]
    fig, ax = plt.subplots(figsize=(7.20, 3.00))
    fig.patch.set_facecolor("white")
    rows = [(s, a) for s in strategy_order for a in algorithm_order]
    matrix = np.full((len(rows), len(metric_order)), np.nan)
    ci_width = np.full_like(matrix, np.nan)
    for ri, (s, a) in enumerate(rows):
        for mi, metric in enumerate(metric_order):
            z = d[(d["strategy"] == s) & (d["algorithm"] == a) & (d["metric"] == metric)]
            if len(z) != 1:
                raise ValueError(f"Expected one row for {s}/{a}/{metric}, found {len(z)}")
            matrix[ri, mi] = float(z.iloc[0]["estimate"])
            ci_width[ri, mi] = float(z.iloc[0]["ci_high"] - z.iloc[0]["ci_low"])

    cmap = LinearSegmentedColormap.from_list(
        "performance", ["#F3F5F6", "#CFDCE8", "#87A5C1", CONTROL]
    )
    im = ax.imshow(matrix, cmap=cmap, vmin=0.50, vmax=1.00, aspect="auto")
    for ri in range(matrix.shape[0]):
        for mi in range(matrix.shape[1]):
            v = matrix[ri, mi]
            color = "white" if v >= 0.90 else INK
            ax.text(mi, ri, f"{v:.3f}", ha="center", va="center",
                    fontsize=6.5, color=color,
                    fontweight="bold" if v >= 0.98 else "normal")

    ax.set_xticks(range(len(metric_order)), metric_labels, fontsize=7.7,
                  fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)), [a for _, a in rows], fontsize=7.0)
    ax.tick_params(length=0, pad=5)
    ax.set_xticks(np.arange(-0.5, len(metric_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i, strategy in enumerate(strategy_order):
        center = i * len(algorithm_order) + 1.5
        ax.text(-0.92, center, strategy_labels[strategy], ha="right", va="center",
                fontsize=7.2, color=INK, fontweight="bold", clip_on=False)
        if i:
            ax.axhline(i * len(algorithm_order) - 0.5, color="#455560", lw=0.9)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.022)
    cbar.set_label("Fold-level mean", fontsize=6.9)
    cbar.ax.tick_params(labelsize=6.5, length=2)
    cbar.outline.set_linewidth(0.5)
    ax.text(0.0, -0.095,
            "Higher is favorable; Brier score is shown as 1 − Brier.",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.4, color=MUTED)
    fig.subplots_adjust(left=0.29, right=0.93, top=0.87, bottom=0.12)
    export(fig, "Fig_ModelPerformance_BalanceStrategies")


def panel_c() -> None:
    """Direct dot plot of S100A8 selection frequency across audit regimes."""
    path = FIG_DIR / "Fig_S100A8_SelectionFrequency_source_data.csv"
    d = pd.read_csv(path)
    regimes = [
        "Weighted nested CV",
        "Unweighted nested CV",
        "Training-fold downsampling",
        "Repeated 1:1 downsampling",
        "Stratified bootstrap",
    ]
    regime_labels = ["Weighted CV", "Unweighted CV", "Fold downsample",
                     "Repeated 1:1", "Bootstrap"]
    algorithms = ["glmBoost", "LASSO", "Ranger", "XGBoost"]
    alg_colors = {"glmBoost": "#355F7D", "LASSO": CONTROL,
                  "Ranger": PRECURSOR, "XGBoost": "#7A7180"}
    offsets = dict(zip(algorithms, [0.24, 0.08, -0.08, -0.24]))

    fig, ax = plt.subplots(figsize=(5.00, 2.05))
    fig.patch.set_facecolor("white")
    ybase = np.arange(len(regimes))[::-1]
    for y in ybase:
        ax.axhline(y, color="#E1E6E9", lw=0.65, zorder=0)
    for algorithm in algorithms:
        z = d[d["algorithm"] == algorithm].set_index("regime").reindex(regimes)
        yy = ybase + offsets[algorithm]
        xx = z["selection_frequency"].to_numpy(float) * 100
        ax.scatter(xx, yy, s=32, color=alg_colors[algorithm],
                   edgecolor="white", linewidth=0.55, zorder=3, label=algorithm)
        for xval, yval in zip(xx, yy):
            if xval < 75:
                ax.text(xval + 2.0, yval, f"{xval:.0f}", ha="left", va="center",
                        fontsize=6.1, color=alg_colors[algorithm])

    ax.axvline(100, color="#9AA6AF", lw=0.7, ls=(0, (2, 2)))
    ax.set_xlim(0, 104)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("S100A8 selection frequency (%)", fontsize=8.3)
    ax.set_yticks(ybase, regime_labels, fontsize=7.6)
    ax.tick_params(axis="x", labelsize=7.2, length=3, width=0.6)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.47, 1.13),
              fontsize=6.8, handletextpad=0.35, columnspacing=0.8)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#6E7B84")
    fig.subplots_adjust(left=0.21, right=0.98, top=0.84, bottom=0.21)
    export(fig, "Fig_S100A8_SelectionFrequency")


def panel_d() -> None:
    """Conventional horizontal boxplots of S100A8 rank stability."""
    path = FIG_DIR / "Fig_S100A8_RankDistribution_source_data.csv"
    d = pd.read_csv(path)
    regimes = ["Weighted nested CV", "Repeated 1:1 downsampling", "Stratified bootstrap"]
    algorithms = ["glmBoost", "LASSO", "Ranger", "XGBoost"]
    regime_colors = {
        "Weighted nested CV": CONTROL,
        "Repeated 1:1 downsampling": PRECURSOR,
        "Stratified bootstrap": "#7A7180",
    }

    fig, ax = plt.subplots(figsize=(5.00, 2.05))
    fig.patch.set_facecolor("white")
    base_positions = np.arange(len(algorithms))[::-1] * 1.15
    offsets = [0.25, 0.0, -0.25]
    for ai, algorithm in enumerate(algorithms):
        for regime, offset in zip(regimes, offsets):
            z = d[(d["algorithm"] == algorithm) & (d["regime"] == regime)]["rank"].dropna().to_numpy(float)
            q025, q25, q50, q75, q975 = np.quantile(z, [0.025, 0.25, 0.5, 0.75, 0.975])
            y = base_positions[ai] + offset
            color = regime_colors[regime]
            ax.plot([q025, q975], [y, y], color=color, lw=0.85, zorder=2)
            ax.plot([q25, q75], [y, y], color=color, lw=5.2,
                    solid_capstyle="butt", zorder=3)
            ax.scatter(q50, y, s=24, facecolor="white", edgecolor=color,
                       linewidth=1.1, zorder=4)

    for y in base_positions:
        ax.axhline(y, color="#E4E8EB", lw=0.55, zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(0.85, 120)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100],
                  ["1", "2", "5", "10", "20", "50", "100"], fontsize=7.8)
    ax.set_xlabel("S100A8 importance rank  ←  more important", fontsize=8.5)
    ax.set_yticks(base_positions, algorithms,
                  fontsize=8.5, fontweight="bold")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", length=3, width=0.7)
    ax.set_ylim(base_positions[-1] - 0.55, base_positions[0] + 0.55)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#65727E")

    short_regime = {"Weighted nested CV":"Weighted CV", "Repeated 1:1 downsampling":"Repeated 1:1",
                    "Stratified bootstrap":"Bootstrap"}
    handles = [mpl.lines.Line2D([0], [0], color=regime_colors[r], lw=5,
                                alpha=0.65, label=short_regime[r]) for r in regimes]
    ax.legend(handles=handles, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.07),
              fontsize=7.0, handlelength=1.7, columnspacing=1.4)
    fig.subplots_adjust(left=0.15, right=0.985, top=0.86, bottom=0.17)
    export(fig, "Fig_S100A8_RankDistribution")


def panel_e() -> None:
    """Plain histogram of balanced-resample AUC values."""
    path = FIG_DIR / "Fig_AUC_Distribution_Downsampling_source_data.csv"
    d = pd.read_csv(path)
    values = d["AUC"].dropna().to_numpy(float)
    q025, q25, q50, q75, q975 = np.quantile(values, [0.025, 0.25, 0.5, 0.75, 0.975])

    xmin = max(0.985, values.min() - 0.0007)
    xmax = 1.00025
    fig, ax = plt.subplots(figsize=(5.00, 2.05))
    fig.patch.set_facecolor("white")
    bins = np.linspace(xmin, 1.00005, 34)
    ax.hist(values, bins=bins, color=CONTROL, edgecolor="white",
            linewidth=0.45, alpha=0.92)
    ax.axvline(q025, color="#7F8B94", lw=0.8, ls=(0, (3, 2)))
    ax.axvline(q975, color="#7F8B94", lw=0.8, ls=(0, (3, 2)))
    ax.axvline(q50, color=CONTROL, lw=1.5)
    ymax = ax.get_ylim()[1]
    ax.text(q50 - 0.00020, ymax * 0.78, f"median {q50:.4f}", ha="right", va="center",
            fontsize=7.2, color=CONTROL, fontweight="bold")
    ax.text((q025 + q975) / 2, ymax * 0.96,
            f"95% interval {q025:.4f}–{q975:.4f}",
            ha="center", va="bottom", fontsize=6.7, color=MUTED)
    ax.text(0.0, 1.04, "1,000 repeated balanced resamples",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.2, color=MUTED)

    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("S100A8 AUC in balanced 1:1 samples", fontsize=8.6)
    ax.set_ylabel("Resamples", fontsize=8.0)
    ax.tick_params(labelsize=7.4, length=3, width=0.7)
    ax.spines["left"].set_color("#65727E")
    ax.spines["bottom"].set_color("#65727E")
    fig.subplots_adjust(left=0.12, right=0.985, top=0.86, bottom=0.22)
    export(fig, "Fig_AUC_Distribution_Downsampling")


def panel_f() -> None:
    """Transparent bootstrap scatter, with intervals and joint median."""
    path = REVISION / "02_class_imbalance" / "Table_Bootstrap_Results.csv"
    d = pd.read_csv(path)
    x = d["AUC"].to_numpy(float)
    y = d["standardized_mean_difference"].to_numpy(float)
    xq = np.quantile(x, [0.025, 0.5, 0.975])
    yq = np.quantile(y, [0.025, 0.5, 0.975])
    corr = np.corrcoef(x, y)[0, 1]

    fig, ax = plt.subplots(figsize=(5.00, 2.20))
    fig.patch.set_facecolor("white")
    ax.axvspan(xq[0], xq[2], color="#E9EFF2", alpha=0.70, zorder=0)
    ax.axhspan(yq[0], yq[2], color="#F1ECE8", alpha=0.65, zorder=0)
    ax.scatter(x, y, s=9, color="#536B7A", alpha=0.18,
               edgecolor="none", rasterized=False, zorder=2)
    ax.scatter(xq[1], yq[1], marker="*", s=125, facecolor=TARGET,
               edgecolor="white", linewidth=0.9, zorder=5)
    ax.axvline(xq[1], color="#4D405B", lw=0.8, ls=(0, (3, 3)), alpha=0.75, zorder=3)
    ax.axhline(yq[1], color="#4D405B", lw=0.8, ls=(0, (3, 3)), alpha=0.75, zorder=3)
    ax.plot([xq[0], xq[2], xq[2], xq[0], xq[0]],
            [yq[0], yq[0], yq[2], yq[2], yq[0]],
            color="#65727E", lw=0.85, ls=(0, (4, 3)), zorder=4)

    ax.text(0.025, 0.965,
            f"Pearson r = {corr:.2f}\n"
            f"AUC 95% interval  {xq[0]:.4f}–{xq[2]:.4f}\n"
            f"SMD 95% interval  {yq[0]:.2f}–{yq[2]:.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.0,
            color=INK, linespacing=1.45,
            bbox={"boxstyle":"round,pad=0.18", "facecolor":"white",
                  "edgecolor":"none", "alpha":0.82})
    ax.annotate("joint median", xy=(xq[1], yq[1]), xytext=(10, 9),
                textcoords="offset points", fontsize=7.0, color="#7A5433",
                fontweight="bold")

    ax.set_xlabel("Bootstrap AUC", fontsize=8.6)
    ax.set_ylabel("Bootstrap SMD", fontsize=8.6)
    # Keep the terminal 1.000 tick safely inside the panel after assembly.
    xspan = max(float(np.ptp(x)), 1e-4)
    ax.set_xlim(float(np.min(x)) - 0.035 * xspan,
                max(1.00005, float(np.max(x)) + 0.055 * xspan))
    ax.tick_params(labelsize=7.7, length=3, width=0.7)
    ax.spines["left"].set_color("#65727E")
    ax.spines["bottom"].set_color("#65727E")
    ax.text(1.0, 1.02, "n = 1,000 bootstrap samples",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.1, color=MUTED)
    fig.subplots_adjust(left=0.16, right=0.95, top=0.91, bottom=0.20)
    export(fig, "Fig_Bootstrap_EffectSize")


PANEL_FUNCS = {
    "A": panel_a,
    "B": panel_b,
    "C": panel_c,
    "D": panel_d,
    "E": panel_e,
    "F": panel_f,
}


def _embedded_svg(stem: str, prefix: str, x: float, y: float,
                  width: float, height: float) -> str:
    """Embed one editable SVG panel without rasterization or cropping."""
    source = (FIG_DIR / f"{stem}.svg").read_text(encoding="utf-8")
    root_match = re.search(r"<svg\b[^>]*viewBox=\"([^\"]+)\"[^>]*>", source)
    if not root_match:
        raise ValueError(f"No SVG viewBox found in {stem}")
    vb = [float(v) for v in root_match.group(1).split()]
    sw, sh = vb[2], vb[3]
    inner = source[root_match.end():source.rfind("</svg>")]

    ids = re.findall(r'\bid="([^"]+)"', inner)
    id_map = {old: f"{prefix}_{old}" for old in ids}
    inner = re.sub(r'\bid="([^"]+)"',
                   lambda m: f'id="{id_map[m.group(1)]}"', inner)
    for old, new in sorted(id_map.items(), key=lambda item: -len(item[0])):
        inner = inner.replace(f"url(#{old})", f"url(#{new})")
        inner = inner.replace(f'href="#{old}"', f'href="#{new}"')
        inner = inner.replace(f"href='#{old}'", f"href='#{new}'")

    scale = min(width / sw, height / sh)
    tx = x + (width - sw * scale) / 2 - vb[0] * scale
    ty = y + (height - sh * scale) / 2 - vb[1] * scale
    return (
        f'<g id="panel_{prefix}" transform="translate({tx:.4f} {ty:.4f}) '
        f'scale({scale:.7f})">{inner}</g>'
    )


def compose_fig2() -> None:
    """Compose the six independent vector panels into final SVG and PDF."""
    stems = {
        "A": "Fig_ClassBalance_StudyDesign",
        "B": "Fig_ModelPerformance_BalanceStrategies",
        "C": "Fig_S100A8_SelectionFrequency",
        "D": "Fig_S100A8_RankDistribution",
        "E": "Fig_AUC_Distribution_Downsampling",
        "F": "Fig_Bootstrap_EffectSize",
    }
    # 7.5 x 7.64 in canvas.  Full-width audit and performance rows preserve
    # 8-12 pt final type; the four diagnostics use balanced half-width cells.
    canvas_w, canvas_h = 540.0, 550.0
    cells = {
        "A": (14, 8, 516, 104),
        "B": (14, 126, 516, 180),
        "C": (14, 320, 252, 100),
        "D": (282, 320, 248, 100),
        "E": (14, 434, 252, 103),
        "F": (282, 434, 248, 103),
    }
    label_xy = {
        "A": (1.5, 20), "B": (1.5, 138), "C": (1.5, 332),
        "D": (269, 332), "E": (1.5, 446), "F": (269, 446),
    }
    pieces = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{canvas_w}pt" '
        f'height="{canvas_h}pt" viewBox="0 0 {canvas_w} {canvas_h}" version="1.1">',
        '<rect width="540" height="550" fill="#ffffff"/>',
    ]
    for label in ["A", "B", "C", "E", "D", "F"]:
        pieces.append(_embedded_svg(stems[label], label, *cells[label]))
        lx, ly = label_xy[label]
        pieces.append(
            f'<text x="{lx}" y="{ly}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="12" font-weight="700" fill="#17232d">{label}</text>'
        )
    pieces.append("</svg>")
    output = FINAL_DIR / "Fig2.svg"
    output.write_text("\n".join(pieces), encoding="utf-8")

    # The PDF compositor is the Illustrator-safe master: each source PDF is
    # placed as vector content, avoiding Adobe's inconsistent interpretation of
    # nested SVG clip paths.  No panel is rasterized or cropped.
    import fitz
    pdf_master = FINAL_DIR / "Fig2_python_master.pdf"
    final_doc = fitz.open()
    page = final_doc.new_page(width=canvas_w, height=canvas_h)
    arial_bold_file = Path(r"C:\Windows\Fonts\arialbd.ttf")
    panel_font = "hebo"
    if arial_bold_file.exists():
        page.insert_font(fontname="ArialBold", fontfile=str(arial_bold_file))
        panel_font = "ArialBold"
    for label in ["A", "B", "C", "E", "D", "F"]:
        x, y, w, h = cells[label]
        source_doc = fitz.open(FIG_DIR / f"{stems[label]}.pdf")
        page.show_pdf_page(fitz.Rect(x, y, x + w, y + h), source_doc, 0,
                           keep_proportion=True, overlay=True)
        source_doc.close()
        lx, ly = label_xy[label]
        page.insert_text((lx, ly), label, fontname=panel_font, fontsize=12,
                         color=(0.09, 0.14, 0.18), overlay=True)
    final_doc.save(pdf_master, garbage=4, deflate=True, clean=True)
    final_doc.close()

    preview_doc = fitz.open(pdf_master)
    preview_page = preview_doc[0]
    zoom = 2250 / preview_page.rect.width
    preview_page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False).save(
        FINAL_DIR / "Fig2_preview.png"
    )
    preview_doc.close()
    print(f"Composed editable Fig2 masters: {output} and {pdf_master}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel", default="A",
        choices=["A", "B", "C", "D", "E", "F", "all", "compose"],
    )
    args = parser.parse_args()
    if args.panel == "compose":
        compose_fig2()
        return
    targets = list(PANEL_FUNCS) if args.panel == "all" else [args.panel]
    missing = [p for p in targets if p not in PANEL_FUNCS]
    if missing:
        raise SystemExit(f"Panel(s) not implemented yet: {', '.join(missing)}")
    for panel in targets:
        PANEL_FUNCS[panel]()
        print(f"Rendered Fig2{panel} independent source panel")


if __name__ == "__main__":
    main()
