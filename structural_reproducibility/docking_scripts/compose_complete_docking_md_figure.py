from __future__ import annotations

from pathlib import Path
import json
import runpy
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import gaussian_kde


# Reuse the already audited inputs, colors and helpers from the compact figure script.
BASE_SCRIPT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal\scripts\compose_docking_md_figure.py")
g = runpy.run_path(str(BASE_SCRIPT))

ROOT = g["ROOT"]
MD_ROOT = g["MD_ROOT"]
ASSETS = g["ASSETS"]
trim_white = g["trim_white"]
boxed = g["boxed"]
letter = g["letter"]
median_range = g["median_range"]

TPL = g["TPL"]
OH16 = g["OH16"]
COMPLEX = g["COMPLEX"]
APO = g["APO"]
REP_COLORS = g["REP_COLORS"]
CHAIN_B = g["CHAIN_B"]
CHAIN_D = g["CHAIN_D"]
NEUTRAL = g["NEUTRAL"]

OUT = ROOT / "publication_figure_v2_complete"
SOURCE = OUT / "source_data"
PREVIEW = OUT / "preview"
ASSET_OUT = OUT / "assets"
for folder in (OUT, SOURCE, PREVIEW, ASSET_OUT):
    folder.mkdir(parents=True, exist_ok=True)

for asset_name in ("compound_structures_vertical_v2.png",
                   "docking_overlay_overview_v3.png",
                   "docking_overlay_closeup_v3.png"):
    shutil.copy2(ASSETS / asset_name, ASSET_OUT / asset_name)

mpl.rcParams.update({
    "font.size": 6.2,
    "axes.labelsize": 6.2,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "legend.fontsize": 5.3,
    "axes.linewidth": 0.70,
})


def md_csv(name: str) -> pd.DataFrame:
    path = MD_ROOT / name
    shutil.copy2(path, SOURCE / name)
    return pd.read_csv(path)


def ligand_lines(ax: plt.Axes, frame: pd.DataFrame, ylabel: str,
                 xlim: tuple[float, float] = (0, 100), legend: bool = False) -> None:
    for rep, color in zip((1, 2, 3), REP_COLORS):
        dat = frame[frame.replicate == rep].sort_values("time_ns").iloc[::5]
        ax.plot(dat.time_ns, dat.display_value, color=color, lw=0.70, label=f"Rep {rep}")
    ax.set_xlim(*xlim)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    if legend:
        ax.legend(loc="upper left", ncol=3, handlelength=1.1,
                  columnspacing=0.5, borderaxespad=0.25)
    boxed(ax)


# Docking tables frozen into the output package.
dock_reps_path = ROOT / "results" / "comparative_analysis" / "matched_local_replicates_with_rmsd.csv"
dock_pose_path = ROOT / "results" / "comparative_analysis" / "representative_pose_cross_compound_comparison.csv"
similarity_path = ROOT / "results" / "comparative_analysis" / "compound_similarity.csv"
for path in (dock_reps_path, dock_pose_path, similarity_path):
    shutil.copy2(path, SOURCE / path.name)
dock_reps = pd.read_csv(dock_reps_path)
dock_pose = pd.read_csv(dock_pose_path)

# Complete accepted v3 MD source-data set used by this figure.
rmsd = md_csv("system_rmsd_replicate_median_range.csv")
distance = md_csv("system_distance_replicate_median_range.csv")
rmsf = md_csv("chain_local_backbone_rmsf_20_100ns.csv")
rg = md_csv("system_rg_replicate_median_range.csv")
sasa = md_csv("system_sasa_replicate_median_range.csv")
inter_hbond = md_csv("system_hbond_replicate_median_range.csv")
lig_rmsd = md_csv("ligand_rmsd.csv")
lig_pocket = md_csv("ligand_pocket.csv")
lig_mindist = md_csv("ligand_mindist.csv")
lig_contacts = md_csv("ligand_contacts.csv")
lig_hbond = md_csv("ligand_hbond.csv")
contact_occ = md_csv("ligand_contact_occupancy_top14.csv")
pca = md_csv("pooled_basis_pc1_pc2_projections.csv")
fel = md_csv("condition_free_energy_landscapes_300K.csv")

fig = plt.figure(figsize=(7.20, 9.55), constrained_layout=False)
gs = fig.add_gridspec(
    5, 4, height_ratios=[1.42, 1.0, 1.0, 1.0, 1.18],
    left=0.055, right=0.982, bottom=0.042, top=0.985,
    wspace=0.43, hspace=0.66,
)

# A–C: equal-width high-resolution image panels.
top = gs[0, :].subgridspec(1, 3, width_ratios=[1.10, 1.00, 1.10], wspace=0.13)
ax_a = fig.add_subplot(top[0, 0])
ax_a.imshow(trim_white(Image.open(ASSETS / "compound_structures_vertical_v2.png"), pad=8))
ax_a.axis("off")
ax_a.text(0.5, -0.035, "Morgan  0.831     MACCS  0.885",
          transform=ax_a.transAxes, ha="center", va="top", fontsize=5.8)

ax_b = fig.add_subplot(top[0, 1])
ax_b.imshow(trim_white(Image.open(ASSETS / "docking_overlay_overview_v3.png"), pad=4))
ax_b.axis("off")
ax_b.legend(handles=[
    Line2D([0], [0], color=CHAIN_B, lw=3, label="Chain B"),
    Line2D([0], [0], color=CHAIN_D, lw=3, label="Chain D"),
    Line2D([0], [0], color=TPL, lw=3, label="Triptolide"),
    Line2D([0], [0], color=OH16, lw=3, label="16-OH-TPL"),
], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.12),
    handlelength=1.2, columnspacing=0.7, labelspacing=0.25,
    borderaxespad=0)

ax_c = fig.add_subplot(top[0, 2])
ax_c.imshow(trim_white(Image.open(ASSETS / "docking_overlay_closeup_v3.png"), pad=4))
ax_c.axis("off")

# Figure-coordinate labels guarantee exact A/B/C horizontal alignment.
top_label_y = 0.982
for ax, label_text in ((ax_a, "A"), (ax_b, "B"), (ax_c, "C")):
    bbox = ax.get_position()
    fig.text(bbox.x0 - 0.018, top_label_y, label_text, fontsize=9,
             fontweight="bold", va="top", ha="left")

# D: seed-matched docking-score difference.
ax = fig.add_subplot(gs[1, 0])
pivot = dock_reps.pivot_table(index=["receptor", "scoring", "seed"],
                              columns="compound", values="best_score_kcal_mol").reset_index()
pivot["difference"] = pivot["16-hydroxytriptolide"] - pivot["triptolide"]
conditions = [("AC", "vina"), ("BD", "vina"), ("AC", "vinardo"), ("BD", "vinardo")]
rng_jitter = np.random.default_rng(20260821)
for idx, (receptor, scoring) in enumerate(conditions):
    values = pivot[(pivot.receptor == receptor) & (pivot.scoring == scoring)].difference.to_numpy()
    ax.scatter(idx + rng_jitter.uniform(-0.08, 0.08, len(values)), values,
               s=8, color=OH16, edgecolor="white", linewidth=0.2, zorder=3)
    ax.plot([idx - 0.15, idx + 0.15], [np.median(values)] * 2, color="black", lw=0.9)
ax.axhline(0, color="#777777", lw=0.6, ls="--")
ax.set_xticks(range(4), ["AC\nVina", "BD\nVina", "AC\nVinardo", "BD\nVinardo"])
ax.set_ylabel(r"Score difference (kcal mol$^{-1}$)")
ax.text(0.02, 0.98, "16-OH-TPL − TPL", transform=ax.transAxes,
        va="top", ha="left", fontsize=5.1, color=NEUTRAL)
boxed(ax); letter(ax, "D")

# E: representative-pose contact overlap and common-scaffold RMSD.
ax = fig.add_subplot(gs[1, 1])
dock_pose = dock_pose.copy()
dock_pose["key"] = dock_pose.receptor + "\n" + dock_pose.scoring.str.capitalize()
order = ["AC\nVina", "BD\nVina", "AC\nVinardo", "BD\nVinardo"]
dock_pose["key"] = pd.Categorical(dock_pose.key, order, ordered=True)
dock_pose = dock_pose.sort_values("key")
bars = ax.bar(np.arange(4), dock_pose.contact_jaccard, width=0.62,
              color="#7FAFA4", edgecolor="black", linewidth=0.4)
for bar, value in zip(bars, dock_pose.common_scaffold_direct_rmsd_A):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.025,
            f"{value:.2f} Å", ha="center", va="bottom", fontsize=4.9)
ax.set_ylim(0, 1.12)
ax.set_xticks(range(4), order)
ax.set_ylabel("Contact Jaccard index")
boxed(ax); letter(ax, "E")

# F–K: all system-level stability and interface metrics.
ax = fig.add_subplot(gs[1, 2])
median_range(ax, rmsd, "Backbone RMSD (nm)", (0, 100), show_legend=True)
ax.axvspan(0, 20, color="#D9D9D9", alpha=0.22, linewidth=0)
ax.axvline(20, color="#888888", lw=0.55, ls="--")
letter(ax, "F")

ax = fig.add_subplot(gs[1, 3])
median_range(ax, distance, "Chain B–D distance (nm)", (20, 100)); letter(ax, "G")

ax = fig.add_subplot(gs[2, 0])
for system, color, label_text in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
    dat = rmsf[rmsf.system == system].sort_values("protein_position")
    x = dat.protein_position.to_numpy()
    ax.fill_between(x, dat.minimum_nm.to_numpy(), dat.maximum_nm.to_numpy(),
                    color=color, alpha=0.15, linewidth=0)
    ax.plot(x, dat.median_nm, color=color, lw=0.78, label=label_text)
ax.axvline(93.5, color="#888888", lw=0.55, ls="--")
ax.set_xlim(1, 186)
ax.set_xlabel("Protein position")
ax.set_ylabel("Backbone RMSF (nm)")
ax.text(47, 0.96, "Chain B", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.0)
ax.text(140, 0.96, "Chain D", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.0)
boxed(ax); letter(ax, "H")

ax = fig.add_subplot(gs[2, 1]); median_range(ax, rg, r"$R_g$ (nm)", (20, 100)); letter(ax, "I")
ax = fig.add_subplot(gs[2, 2]); median_range(ax, sasa, r"SASA (nm$^2$)", (20, 100)); letter(ax, "J")
ax = fig.add_subplot(gs[2, 3]); median_range(ax, inter_hbond, "Inter-chain H-bonds", (20, 100)); letter(ax, "K")

# L–P: complete ligand-behaviour metrics for all three TPL replicas.
ax = fig.add_subplot(gs[3, 0]); ligand_lines(ax, lig_rmsd, "Unwrapped ligand RMSD (nm)", legend=True); letter(ax, "L")
ax = fig.add_subplot(gs[3, 1]); ligand_lines(ax, lig_pocket, "Pocket displacement (nm)"); letter(ax, "M")
ax = fig.add_subplot(gs[3, 2]); ligand_lines(ax, lig_mindist, "Protein minimum distance (nm)"); letter(ax, "N")
ax = fig.add_subplot(gs[3, 3]); ligand_lines(ax, lig_contacts, "Protein contacts"); letter(ax, "O")

ax = fig.add_subplot(gs[4, 0]); ligand_lines(ax, lig_hbond, "Ligand–protein H-bonds", (20, 100)); letter(ax, "P")

# Q: per-residue contact occupancy for each independent complex replica.
ax = fig.add_subplot(gs[4, 1])
matrix = contact_occ[["occupancy_rep1", "occupancy_rep2", "occupancy_rep3"]].to_numpy()
im = ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1,
               interpolation="nearest")
ax.set_xticks(range(3), ["Rep 1", "Rep 2", "Rep 3"])
ax.set_yticks(range(len(contact_occ)), contact_occ.label, fontsize=4.35)
ax.tick_params(axis="both", direction="in", top=True, right=True, length=1.8, pad=1.2)
for row in range(matrix.shape[0]):
    for col in range(matrix.shape[1]):
        value = matrix[row, col]
        ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=3.8,
                color="white" if value > 0.58 else "#222222")
for spine in ax.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.7)
letter(ax, "Q")

# R: pooled-basis PCA projection as crisp probability-density contours.
ax = fig.add_subplot(gs[4, 2])
x_grid = np.linspace(pca.pc1_nm.min() - 0.35, pca.pc1_nm.max() + 0.35, 150)
y_grid = np.linspace(pca.pc2_nm.min() - 0.35, pca.pc2_nm.max() + 0.35, 150)
xx, yy = np.meshgrid(x_grid, y_grid)
density_rows = []
condition_handles = []
for system, color, label_text in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
    dat = pca[pca.system == system]
    values = np.vstack([dat.pc1_nm.to_numpy(), dat.pc2_nm.to_numpy()])
    density = gaussian_kde(values)(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    ordered = np.sort(density.ravel())[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    thresholds = [ordered[np.searchsorted(cumulative, probability)]
                  for probability in (0.95, 0.75, 0.50)]
    ax.contour(xx, yy, density, levels=thresholds,
               colors=[color] * 3, linestyles=[":", "--", "-"],
               linewidths=[0.75, 0.95, 1.25], zorder=2)
    for rep, marker in zip((1, 2, 3), ("o", "s", "^")):
        rep_dat = dat[dat.replicate == rep]
        ax.scatter(rep_dat.pc1_nm.mean(), rep_dat.pc2_nm.mean(), s=19,
                   marker=marker, color=color, edgecolor="white", linewidth=0.45,
                   zorder=4)
    condition_handles.append(Line2D([0], [0], color=color, lw=1.2, label=label_text))
    for x_value, y_value, density_value in zip(xx.ravel(), yy.ravel(), density.ravel()):
        density_rows.append({"system": system, "pc1_nm": x_value,
                             "pc2_nm": y_value, "density": density_value})
pd.DataFrame(density_rows).to_csv(SOURCE / "pooled_basis_pca_density_grid.csv",
                                  index=False, encoding="utf-8-sig")
ax.set_xlabel("PC1 (36.6%)")
ax.set_ylabel("PC2 (15.1%)")
ax.legend(handles=condition_handles, loc="upper right", ncol=1,
          handlelength=1.3, handletextpad=0.35, borderaxespad=0.2)
boxed(ax); letter(ax, "R")

# S: condition-specific free-energy landscapes on the same pooled PCA basis.
sub = gs[4, 3].subgridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.18)
fel_axes = [fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])]
groups = ["complex_all_replicates", "apo_all_replicates"]
names = ["Complex", "Apo"]
mesh = None
for sub_ax, group, name in zip(fel_axes, groups, names):
    dat = fel[fel.group == group]
    grid = dat.pivot(index="pc2_center_nm", columns="pc1_center_nm",
                     values="free_energy_kj_mol").sort_index()
    values = np.ma.masked_invalid(grid.to_numpy())
    mesh = sub_ax.pcolormesh(grid.columns.to_numpy(), grid.index.to_numpy(), values,
                             cmap="viridis_r", vmin=0, vmax=12, shading="auto")
    sub_ax.text(0.04, 0.96, name, transform=sub_ax.transAxes, va="top", ha="left", fontsize=4.9)
    sub_ax.set_xlabel("PC1")
    boxed(sub_ax)
fel_axes[0].set_ylabel("PC2")
fel_axes[1].set_yticklabels([])
cax = fig.add_subplot(sub[0, 2])
cb = fig.colorbar(mesh, cax=cax)
cb.set_ticks([0, 6, 12])
cb.ax.set_title(r"$\Delta G$", fontsize=4.2, pad=2)
cb.ax.tick_params(labelsize=4.2, length=1.5)
bbox_s = fel_axes[0].get_position()
fig.text(bbox_s.x0 - 0.020, bbox_s.y1 + 0.012, "S", fontsize=9,
         fontweight="bold", va="bottom", ha="left")

stem = "FINAL_Figure_S100A8_docking_and_triptolide_MD"
fig.savefig(OUT / f"{stem}.svg", dpi=600, facecolor="white")
fig.savefig(OUT / f"{stem}.pdf", dpi=600, facecolor="white")
fig.savefig(OUT / f"{stem}.tiff", dpi=600, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white")
fig.savefig(PREVIEW / f"{stem}.png", dpi=230, facecolor="white")
plt.close(fig)

manifest = {
    "figure": stem,
    "dimensions_mm": {"width": 182.88, "height": 242.57},
    "tiff_dpi": 600,
    "panel_labels": list("ABCDEFGHIJKLMNOPQRS"),
    "complete_md_panels": [
        "backbone RMSD", "chain B-D distance", "chain-local backbone RMSF",
        "radius of gyration", "protein SASA", "inter-chain hydrogen bonds",
        "unwrapped ligand RMSD", "original-pocket displacement",
        "protein minimum distance", "protein contacts", "ligand-protein hydrogen bonds",
        "per-residue contact occupancy", "PCA", "free-energy landscapes",
    ],
    "md_replicates": "complex 3 x 100 ns and apo 3 x 100 ns",
    "interpretive_boundary": "All MD panels concern triptolide only; no 16-hydroxytriptolide MD is implied.",
}
(OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(OUT)
