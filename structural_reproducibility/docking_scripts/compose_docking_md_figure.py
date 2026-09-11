from __future__ import annotations

from pathlib import Path
import json
import shutil
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
MD_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_MD_topjournal\comparative_analysis_v3\publication_figures_python_v2\source_data")
OUT = ROOT / "publication_figure_v1"
ASSETS = OUT / "assets"
SOURCE = OUT / "source_data"
PREVIEW = OUT / "preview"
for folder in (OUT, ASSETS, SOURCE, PREVIEW):
    folder.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 6.5,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 5.8,
    "ytick.labelsize": 5.8,
    "legend.fontsize": 5.7,
    "axes.linewidth": 0.72,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

TPL = "#E8892F"
OH16 = "#A94B8C"
COMPLEX = "#B64342"
APO = "#3775BA"
REP_COLORS = ["#4B4B7A", "#36939A", "#B64342"]
CHAIN_B = "#6FA6C1"
CHAIN_D = "#8DC7B8"
NEUTRAL = "#666666"


def trim_white(image: Image.Image, pad: int = 10) -> Image.Image:
    rgb = np.asarray(image.convert("RGB"))
    mask = np.any(rgb < 248, axis=2)
    if not mask.any():
        return image
    yy, xx = np.where(mask)
    left = max(0, int(xx.min()) - pad)
    top = max(0, int(yy.min()) - pad)
    right = min(image.width, int(xx.max()) + pad + 1)
    bottom = min(image.height, int(yy.max()) + pad + 1)
    return image.crop((left, top, right, bottom))


def boxed(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.72)
        spine.set_color("black")
    ax.tick_params(direction="in", length=2.5, width=0.65, pad=1.8,
                   top=True, right=True)
    ax.grid(False)


def letter(ax: plt.Axes, value: str, x: float = -0.16, y: float = 1.04) -> None:
    ax.text(x, y, value, transform=ax.transAxes, fontsize=9.0,
            fontweight="bold", va="bottom", ha="left", clip_on=False)


def load_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(MD_ROOT / name)


def median_range(ax: plt.Axes, frame: pd.DataFrame, ylabel: str,
                 xlim: tuple[float, float], show_legend: bool = False) -> None:
    for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
        dat = frame[frame.system == system].sort_values("time_ns").iloc[::10]
        x = dat.time_ns.to_numpy()
        ax.fill_between(x, dat.minimum.to_numpy(), dat.maximum.to_numpy(),
                        color=color, alpha=0.15, linewidth=0)
        ax.plot(x, dat["median"], color=color, lw=0.85, label=label)
    ax.set_xlim(*xlim)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    if show_legend:
        ax.legend(loc="upper left", ncol=2, handlelength=1.5,
                  columnspacing=0.8, borderaxespad=0.3)
    boxed(ax)


chem_asset = ASSETS / "compound_structures.png"
vertical_asset = ASSETS / "compound_structures_vertical_v2.png"
if not chem_asset.exists() and vertical_asset.exists():
    shutil.copy2(vertical_asset, chem_asset)
if not chem_asset.exists():
    from rdkit import Chem
    from rdkit.Chem import Draw
    # RDKit on this Windows build cannot open Unicode paths, so use read-only ASCII-path copies.
    chem_tmp = Path(tempfile.mkdtemp(prefix="s100a8_compounds_"))
    tpl_tmp = chem_tmp / "triptolide.sdf"
    oh_tmp = chem_tmp / "hydroxytriptolide.sdf"
    shutil.copy2(TPL_ROOT / "input" / "prepared" / "triptolide_primary.sdf", tpl_tmp)
    shutil.copy2(ROOT / "input" / "prepared" / "16-hydroxytriptolide_primary.sdf", oh_tmp)
    tpl_mol = Chem.SDMolSupplier(str(tpl_tmp), removeHs=True)[0]
    oh_mol = Chem.SDMolSupplier(str(oh_tmp), removeHs=True)[0]
    if tpl_mol is None or oh_mol is None:
        raise RuntimeError("Could not read one or both primary ligand SDF files")
    chem_img = Draw.MolsToGridImage(
        [tpl_mol, oh_mol], molsPerRow=1, subImgSize=(700, 310),
        legends=["Triptolide", "16-hydroxytriptolide"], useSVG=False,
    )
    chem_img.save(chem_asset, dpi=(600, 600))
    shutil.rmtree(chem_tmp, ignore_errors=True)

# Freeze/copy exactly the tabular inputs used by this composite.
dock_reps_src = ROOT / "results" / "comparative_analysis" / "matched_local_replicates_with_rmsd.csv"
dock_pose_src = ROOT / "results" / "comparative_analysis" / "representative_pose_cross_compound_comparison.csv"
similarity_src = ROOT / "results" / "comparative_analysis" / "compound_similarity.csv"
for path in (dock_reps_src, dock_pose_src, similarity_src):
    shutil.copy2(path, SOURCE / path.name)

md_files = {
    "rmsd": "system_rmsd_replicate_median_range.csv",
    "rmsf": "chain_local_backbone_rmsf_20_100ns.csv",
    "rg": "system_rg_replicate_median_range.csv",
    "sasa": "system_sasa_replicate_median_range.csv",
    "hbond": "system_hbond_replicate_median_range.csv",
    "pocket": "ligand_pocket.csv",
}
for filename in md_files.values():
    shutil.copy2(MD_ROOT / filename, SOURCE / filename)

dock_reps = pd.read_csv(dock_reps_src)
dock_pose = pd.read_csv(dock_pose_src)
rmsd = load_csv(md_files["rmsd"])
rmsf = load_csv(md_files["rmsf"])
rg = load_csv(md_files["rg"])
sasa = load_csv(md_files["sasa"])
hbond = load_csv(md_files["hbond"])
pocket = load_csv(md_files["pocket"])

fig = plt.figure(figsize=(7.20, 6.62), constrained_layout=False)
gs = fig.add_gridspec(3, 4, height_ratios=[1.23, 1.0, 1.0],
                      left=0.055, right=0.988, bottom=0.065, top=0.988,
                      wspace=0.42, hspace=0.60)

# A: two compounds and numerical similarity.
ax = fig.add_subplot(gs[0, 0])
ax.imshow(trim_white(Image.open(chem_asset), pad=15))
ax.axis("off")
ax.text(0.5, -0.025, "Morgan similarity  0.831\nMACCS similarity  0.885",
        transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
        linespacing=1.25)
letter(ax, "A", x=-0.08, y=1.0)

# B: whole dimer/interface view.
ax = fig.add_subplot(gs[0, 1:3])
ax.imshow(trim_white(Image.open(ASSETS / "docking_overlay_overview_v2.png"), pad=5))
ax.axis("off")
ax.legend(handles=[
    Line2D([0], [0], color=CHAIN_B, lw=3, label="Chain B"),
    Line2D([0], [0], color=CHAIN_D, lw=3, label="Chain D"),
    Line2D([0], [0], color=TPL, lw=3, label="Triptolide"),
    Line2D([0], [0], color=OH16, lw=3, label="16-OH-TPL"),
], loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.04),
          handlelength=1.3, columnspacing=0.8)
letter(ax, "B", x=-0.04, y=0.94)

# C: shared interface pocket close-up.
ax = fig.add_subplot(gs[0, 3])
ax.imshow(trim_white(Image.open(ASSETS / "docking_overlay_closeup_v2.png"), pad=5))
ax.axis("off")
letter(ax, "C", x=-0.08, y=1.0)

# D: seed-matched score differences; positive values mean less-negative score for 16-OH-TPL.
ax = fig.add_subplot(gs[1, 0])
pivot = dock_reps.pivot_table(index=["receptor", "scoring", "seed"],
                              columns="compound", values="best_score_kcal_mol").reset_index()
pivot["difference"] = pivot["16-hydroxytriptolide"] - pivot["triptolide"]
conditions = [("AC", "vina"), ("BD", "vina"), ("AC", "vinardo"), ("BD", "vinardo")]
rng = np.random.default_rng(20260821)
for idx, (receptor, scoring) in enumerate(conditions):
    values = pivot[(pivot.receptor == receptor) & (pivot.scoring == scoring)].difference.to_numpy()
    jitter = rng.uniform(-0.085, 0.085, size=len(values))
    ax.scatter(np.full(len(values), idx) + jitter, values, s=10,
               color=OH16, edgecolor="white", linewidth=0.25, zorder=3)
    ax.plot([idx - 0.16, idx + 0.16], [np.median(values)] * 2, color="black", lw=1.0, zorder=4)
ax.axhline(0, color="#777777", lw=0.65, ls="--")
ax.set_xticks(range(4), ["AC\nVina", "BD\nVina", "AC\nVinardo", "BD\nVinardo"])
ax.set_ylabel(r"Score difference (kcal mol$^{-1}$)")
ax.text(0.02, 0.98, "16-OH-TPL − TPL", transform=ax.transAxes,
        ha="left", va="top", fontsize=5.5, color=NEUTRAL)
boxed(ax); letter(ax, "D")

# E: representative-pose contact overlap, with common-scaffold RMSD annotated.
ax = fig.add_subplot(gs[1, 1])
dock_pose = dock_pose.copy()
dock_pose["key"] = dock_pose.receptor + "\n" + dock_pose.scoring.str.capitalize()
order = ["AC\nVina", "BD\nVina", "AC\nVinardo", "BD\nVinardo"]
dock_pose["key"] = pd.Categorical(dock_pose.key, order, ordered=True)
dock_pose = dock_pose.sort_values("key")
bars = ax.bar(np.arange(4), dock_pose.contact_jaccard, width=0.62,
              color="#7FAFA4", edgecolor="black", linewidth=0.45)
for bar, rmsd_value in zip(bars, dock_pose.common_scaffold_direct_rmsd_A):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.025,
            f"{rmsd_value:.2f} Å", ha="center", va="bottom", fontsize=5.2)
ax.set_ylim(0, 1.12)
ax.set_xticks(range(4), order)
ax.set_ylabel("Contact Jaccard index")
boxed(ax); letter(ax, "E")

# F: whole-protein backbone RMSD, full range across three independent trajectories.
ax = fig.add_subplot(gs[1, 2])
median_range(ax, rmsd, "Backbone RMSD (nm)", (0, 100), show_legend=True)
ax.axvspan(0, 20, color="#D9D9D9", alpha=0.22, linewidth=0)
ax.axvline(20, color="#888888", lw=0.6, ls="--")
letter(ax, "F")

# G: chain-local RMSF (20–100 ns).
ax = fig.add_subplot(gs[1, 3])
for system, color, label_text in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
    dat = rmsf[rmsf.system == system].sort_values("protein_position")
    x = dat.protein_position.to_numpy()
    ax.fill_between(x, dat.minimum_nm.to_numpy(), dat.maximum_nm.to_numpy(),
                    color=color, alpha=0.15, linewidth=0)
    ax.plot(x, dat.median_nm, color=color, lw=0.85, label=label_text)
ax.axvline(93.5, color="#888888", lw=0.6, ls="--")
ax.set_xlim(1, 186)
ax.set_xlabel("Protein position")
ax.set_ylabel("Backbone RMSF (nm)")
ax.text(47, 0.96, "Chain B", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.3)
ax.text(140, 0.96, "Chain D", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.3)
boxed(ax); letter(ax, "G")

# H–J: compactness, solvent exposure and inter-chain hydrogen bonds.
ax = fig.add_subplot(gs[2, 0])
median_range(ax, rg, r"$R_g$ (nm)", (20, 100)); letter(ax, "H")

ax = fig.add_subplot(gs[2, 1])
median_range(ax, sasa, r"SASA (nm$^2$)", (20, 100)); letter(ax, "I")

ax = fig.add_subplot(gs[2, 2])
median_range(ax, hbond, "Inter-chain H-bonds", (20, 100)); letter(ax, "J")

# K: ligand displacement from the original docking pocket; all three TPL replicas are shown.
ax = fig.add_subplot(gs[2, 3])
for rep, color in zip((1, 2, 3), REP_COLORS):
    dat = pocket[pocket.replicate == rep].sort_values("time_ns").iloc[::5]
    ax.plot(dat.time_ns, dat.display_value, color=color, lw=0.78, label=f"Rep {rep}")
ax.set_xlim(0, 100)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Pocket displacement (nm)")
ax.legend(loc="upper left", ncol=3, handlelength=1.2, columnspacing=0.65,
          borderaxespad=0.3)
boxed(ax); letter(ax, "K")

stem = "Figure_docking_comparison_and_triptolide_MD"
fig.savefig(OUT / f"{stem}.svg", dpi=600, facecolor="white")
fig.savefig(OUT / f"{stem}.pdf", dpi=600, facecolor="white")
fig.savefig(OUT / f"{stem}.tiff", dpi=600, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(PREVIEW / f"{stem}.png", dpi=250, facecolor="white")
plt.close(fig)

manifest = {
    "figure": stem,
    "dimensions_mm": {"width": 182.88, "height": 168.15},
    "tiff_dpi": 600,
    "panel_labels": list("ABCDEFGHIJK"),
    "docking_receptor": "PDB 5HLO S100A8 homodimer interfaces AC and BD",
    "docking_compounds": ["triptolide", "16-hydroxytriptolide"],
    "md_scope": "triptolide-S100A8 complex only; 3 x 100 ns plus apo 3 x 100 ns",
    "md_pbc": "cluster-to-nojump corrected v3 trajectories",
    "plotting_backend": "Python/matplotlib; RDKit for 2D structures; PyMOL for 3D assets",
    "interpretive_boundary": "No 16-hydroxytriptolide MD or wet-lab evidence is implied.",
}
(OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(OUT)
