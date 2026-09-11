"""Redraw Fig. 3 panels as six independent, single-chart Python figures."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap


ROOT = Path(r"D:\桌面\sepsis")
Q = ROOT / "review_revision" / "09_figures"
CLIN = ROOT / "review_revision" / "03_clinical_validation"
EXT = ROOT / "sepsis-s100a8-multiomics" / "data" / "08_External Validation" / "output"

INK = "#26323A"
MUTED = "#687780"
GRID = "#E8ECEE"
SEPSIS = "#D66A5E"
HEALTHY = "#4C78A8"
SIMILAR = "#D49A45"
ROC_COLORS = ["#355F7D", "#4C78A8", "#D49A45", "#D66A5E"]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
})


def save(fig, stem, directory=Q):
    for item in fig.findobj(match=mpl.text.Text):
        if item.get_text().strip() and item.get_fontsize() < 14.5:
            item.set_fontsize(14.5)
    try:
        fig.tight_layout(pad=0.45)
    except Exception:
        pass
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf"):
        fig.savefig(directory / f"{stem}.{ext}")
    fig.savefig(directory / f"{stem}.tiff", dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(directory / f"{stem}.png", dpi=300)
    plt.close(fig)


def clean(ax, xgrid=False):
    ax.tick_params(length=3, width=0.7)
    if xgrid:
        ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def panel_a():
    d = pd.read_csv(Q / "Fig3A_ExternalROC_source_data.csv")
    order = ["GSE9692", "GSE26440", "GSE28750", "GSE69528"]
    fig, ax = plt.subplots(figsize=(4.25, 3.55))
    for color, name in zip(ROC_COLORS, order):
        z = d[d.dataset == name].sort_values("specificity", ascending=False)
        auc = z.AUC.iloc[0]
        ax.step(1 - z.specificity, z.sensitivity, where="post", lw=1.8,
                color=color, label=f"{name}   AUC {auc:.2f}")
    ax.plot([0, 1], [0, 1], color="#B4BDC2", lw=0.8, ls=(0, (3, 3)))
    ax.set(xlim=(-0.01, 1.01), ylim=(-0.01, 1.01),
           xlabel="1 − specificity", ylabel="Sensitivity")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", fontsize=7.4, handlelength=2.2, labelspacing=0.45)
    clean(ax)
    save(fig, "Combined_ROC_Validation_Final", EXT)


def panel_b():
    d = pd.read_csv(Q / "Fig_ClinicalComparator_EffectForest_source_data.csv")
    d = d.sort_values(["family", "estimate"], ascending=[True, False]).reset_index(drop=True)
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(6.25, 4.25))
    for yi, row in zip(y, d.itertuples()):
        col = HEALTHY if row.family == "Healthy case-control" else SIMILAR
        ax.plot([row.ci_low, row.ci_high], [yi, yi], color=col, lw=1.5, solid_capstyle="round")
        ax.scatter(row.estimate, yi, s=31, color=col, edgecolor="white", linewidth=0.55, zorder=3)
    ax.axvline(0, color="#8E979C", lw=0.9, ls=(0, (3, 3)))
    ax.set_yticks(y, d.label.str.replace(" vs ", "  |  ", regex=False))
    ax.set_xlabel("Standardized mean difference (95% CI)")
    ax.set_xlim(min(-1.55, d.ci_low.min() - .15), d.ci_high.max() + .95)
    ax.scatter([], [], color=HEALTHY, s=27, label="Healthy")
    ax.scatter([], [], color=SIMILAR, s=27, label="Similar")
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.015), ncol=2, fontsize=7.1)
    clean(ax, xgrid=True)
    save(fig, "Fig_ClinicalComparator_EffectForest")


COMPARISONS = [
    ("GSE131411", "Cardiogenic shock", "GSE131411 | cardiogenic shock", "similar"),
    ("GSE134347", "Non-infectious ICU", "GSE134347 | non-infectious ICU", "similar"),
    ("GSE28750", "Postoperative sterile inflammation", "GSE28750 | postoperative", "similar"),
    ("GSE69528", "Uninfected chronic disease control", "GSE69528 | chronic disease", "similar"),
    ("GSE134347", "Healthy", "GSE134347 | healthy", "healthy"),
    ("GSE26440", "Healthy", "GSE26440 | healthy", "healthy"),
    ("GSE28750", "Healthy", "GSE28750 | healthy", "healthy"),
    ("GSE69528", "Healthy", "GSE69528 | healthy", "healthy"),
    ("GSE9692", "Healthy", "GSE9692 | healthy", "healthy"),
]


def panel_c():
    raw = pd.read_csv(CLIN / "clinical_validation_S100A8_samples.csv")
    rows = []
    for dataset, comparator, label, family in COMPARISONS:
        z = raw[(raw.dataset == dataset) & raw.group.isin(["Sepsis", comparator])].copy()
        z["z"] = (z.expression - z.expression.mean()) / z.expression.std(ddof=1)
        z["label"], z["family"] = label, family
        rows.append(z)
    d = pd.concat(rows, ignore_index=True)
    pd.DataFrame(d).to_csv(Q / "Fig_ClinicalComparator_Expression_source_data.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    rng = np.random.default_rng(20260824)
    labels = [x[2] for x in COMPARISONS]
    ybase = np.arange(len(labels))[::-1]
    for yi, (_, comparator, label, family) in zip(ybase, COMPARISONS):
        z = d[d.label == label]
        for group, off, col in [(comparator, -0.14, HEALTHY if family == "healthy" else SIMILAR),
                                ("Sepsis", 0.14, SEPSIS)]:
            vals = z.loc[z.group == group, "z"].to_numpy(float)
            jitter = rng.normal(0, 0.025, len(vals))
            ax.scatter(vals, np.full(len(vals), yi + off) + jitter, s=8, color=col,
                       alpha=0.22, edgecolors="none", rasterized=False)
            q1, med, q3 = np.quantile(vals, [.25, .5, .75])
            ax.plot([q1, q3], [yi + off, yi + off], color=col, lw=3.4, solid_capstyle="round")
            ax.scatter(med, yi + off, s=24, color=col, edgecolor="white", linewidth=.5, zorder=4)
    ax.axvline(0, color="#A8B0B5", lw=.8, ls=(0, (3, 3)))
    ax.set_yticks(ybase, labels)
    ax.set_xlabel("S100A8 expression z score")
    ax.scatter([], [], color=SEPSIS, s=25, label="Sepsis")
    ax.scatter([], [], color=HEALTHY, s=25, label="Healthy")
    ax.scatter([], [], color=SIMILAR, s=25, label="Clin. similar")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(.5, 1.01), fontsize=7.1)
    clean(ax, xgrid=True)
    save(fig, "Fig_ClinicalComparator_Expression")


def panel_d():
    table = pd.read_csv(CLIN / "Table_ClinicalValidation_long.csv")
    d = table[table.metric == "AUC"].copy()
    d["family"] = np.where(d.comparator.eq("Healthy"), "Healthy comparator", "Clinically similar comparator")
    d["label"] = d.dataset + "  |  " + d.comparator.str.replace("Uninfected chronic disease control", "chronic disease", regex=False).str.replace("Postoperative sterile inflammation", "postoperative", regex=False)
    d = d.sort_values(["family", "estimate"], ascending=[True, False]).reset_index(drop=True)
    d.to_csv(Q / "Fig_ClinicalComparator_ROC_source_data.csv", index=False)
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(6.15, 4.25))
    ax.axvspan(.5, .7, color="#F3F5F6", zorder=0)
    ax.axvline(.5, color="#90999E", lw=.9, ls=(0, (3, 3)))
    ax.axvline(.8, color="#C1C8CC", lw=.7, ls=(0, (2, 3)))
    for yi, row in zip(y, d.itertuples()):
        col = HEALTHY if row.family == "Healthy comparator" else SIMILAR
        lo, hi = row.ci_low, min(1.02, row.ci_high)
        ax.plot([lo, hi], [yi, yi], color=col, lw=1.5)
        ax.scatter(row.estimate, yi, color=col, s=31, edgecolor="white", linewidth=.5, zorder=3)
    ax.set_yticks(y, d.label)
    ax.set_xlim(max(.18, d.ci_low.min() - .03), 1.075)
    ax.set_xlabel("ROC AUC (95% CI)")
    ax.text(.5, len(d)-.15, "chance", ha="center", va="bottom", fontsize=6.6, color=MUTED)
    clean(ax, xgrid=True)
    save(fig, "Fig_ClinicalComparator_ROC")


def panel_e():
    d = pd.read_csv(Q / "Fig_DecisionCurve_source_data.csv")
    labels = {
        "GSE131411_CardiogenicShock": "Cardiogenic shock",
        "GSE134347_NoninfectiousICU": "Non-infectious ICU",
        "GSE28750_Postoperative": "Postoperative sterile inflammation",
    }
    rows = []
    for cid, label in labels.items():
        z = d[d.comparison_id == cid].pivot(index="threshold", columns="strategy", values="net_benefit").sort_index()
        biom = z["S100A8 (cross-validated)"]
        ref = np.maximum(z["Treat all"], z["Treat none"])
        rows.append(pd.DataFrame({"comparison_id": cid, "label": label,
                                  "threshold": z.index, "delta_net_benefit": biom - ref}))
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(Q / "Fig3E_NetBenefitAdvantage_source_data.csv", index=False)
    heat = out.pivot(index="label", columns="threshold", values="delta_net_benefit").loc[list(labels.values())]
    vals = heat.to_numpy(float)
    lim = max(.08, np.nanquantile(np.abs(vals), .97))
    cmap = LinearSegmentedColormap.from_list("clinical_delta", [SEPSIS, "#F7F7F5", HEALTHY])
    cmap.set_bad("#ECEFF1")
    fig, ax = plt.subplots(figsize=(6.1, 3.15))
    im = ax.imshow(vals, aspect="auto", interpolation="nearest", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
                   extent=[heat.columns.min(), heat.columns.max(), len(heat)-.5, -.5])
    ax.set_yticks(range(len(heat)), heat.index)
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("")
    cbar = fig.colorbar(im, ax=ax, fraction=.035, pad=.025)
    # A horizontal colorbar title avoids collision with the vertical tick labels.
    cbar.ax.set_title("ΔNB", fontsize=7.5, pad=10)
    cbar.ax.tick_params(labelsize=6.8, length=2)
    ax.spines.top.set_visible(False); ax.spines.right.set_visible(False)
    clean(ax)
    save(fig, "Fig_DecisionCurve")


def panel_f():
    d = pd.read_csv(Q / "Fig_S100A8_Prognosis_Extended_source_data.csv")
    d = d.sort_values(["dataset", "endpoint"]).reset_index(drop=True)
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(6.15, 3.45))
    ax.axvline(1, color="#8F989D", lw=.9, ls=(0, (3, 3)))
    for yi, row in zip(y, d.itertuples()):
        ax.plot([row.ci_low, row.ci_high], [yi, yi], color=MUTED, lw=1.5)
        ax.scatter(row.estimate, yi, s=34, color=SEPSIS, edgecolor="white", linewidth=.55, zorder=3)
    labels = [f"{r.dataset}  |  {r.endpoint}\n{int(r.events)}/{int(r.n)} events" for r in d.itertuples()]
    ax.set_yticks(y, labels)
    ax.set_xscale("log")
    ax.set_xlim(.24, 2.65)
    ax.set_xticks([.25, .5, 1, 2], ["0.25", "0.5", "1", "2"])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("OR per 1-SD S100A8")
    for yi, row in zip(y, d.itertuples()):
        ax.text(2.55, yi, f"P={row.p_value:.3f}", ha="right", va="center", fontsize=7, color=MUTED)
    clean(ax, xgrid=True)
    save(fig, "Fig_S100A8_Prognosis_Extended")


PANELS = {"A": panel_a, "B": panel_b, "C": panel_c, "D": panel_d, "E": panel_e, "F": panel_f}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=[*PANELS, "all"], default="all")
    args = ap.parse_args()
    if args.panel == "all":
        for fn in PANELS.values():
            fn()
    else:
        PANELS[args.panel]()
