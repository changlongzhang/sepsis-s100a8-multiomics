from __future__ import annotations

from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "comparative_analysis_v3"
OUT = BASE / "publication_figures_python_v2"
SOURCE = OUT / "source_data"
PREVIEW = OUT / "preview"
for folder in (OUT, SOURCE, PREVIEW):
    folder.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

COMPLEX = "#D56A54"
APO = "#315B78"
REP_COLORS = ["#315B78", "#5B958C", "#D56A54"]
NEUTRAL = "#6E6E6E"


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip()
        if not text or text[0] in "#@" or text == "&":
            continue
        try:
            rows.append([float(value) for value in text.split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric data: {path}")
    return np.asarray(rows, dtype=float)


def load_series(path: Path, start: float | None = None, end: float | None = None) -> pd.DataFrame:
    arr = read_xvg(path)
    time = arr[:, 0] / 1000.0 if np.nanmax(arr[:, 0]) > 1000 else arr[:, 0]
    frame = pd.DataFrame({"time_ns": time, "value": arr[:, 1]})
    if start is not None:
        frame = frame[frame.time_ns >= start]
    if end is not None:
        frame = frame[frame.time_ns <= end + 1e-6]
    return frame.reset_index(drop=True)


def run_path(system: str, rep: int, filename: str) -> Path:
    return BASE / f"{system}_rep{rep}" / filename


def aggregate(filename: str, start: float, end: float, smooth: int = 1) -> dict[str, pd.DataFrame]:
    result = {}
    for system in ("complex", "apo"):
        frames = [load_series(run_path(system, rep, filename), start, end) for rep in range(1, 4)]
        time = frames[0].time_ns.to_numpy()
        values = []
        for frame in frames:
            if len(frame) != len(time) or not np.allclose(frame.time_ns, time, atol=1e-4):
                values.append(np.interp(time, frame.time_ns, frame.value))
            else:
                values.append(frame.value.to_numpy())
        matrix = pd.DataFrame(np.vstack(values).T)
        if smooth > 1:
            matrix = matrix.rolling(smooth, center=True, min_periods=max(1, smooth // 4)).mean()
        result[system] = pd.DataFrame({
            "system": system,
            "time_ns": time,
            "median": matrix.median(axis=1),
            "minimum": matrix.min(axis=1),
            "maximum": matrix.max(axis=1),
            "rep1": matrix.iloc[:, 0],
            "rep2": matrix.iloc[:, 1],
            "rep3": matrix.iloc[:, 2],
        })
    return result


def write_source(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(SOURCE / filename, index=False, encoding="utf-8-sig")


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.16, 1.04, letter.upper(), transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom")
    ax.tick_params(direction="out", length=2.5, width=0.6, pad=2)


def median_range(ax: plt.Axes, data: dict[str, pd.DataFrame], ylabel: str,
                 xlim: tuple[float, float], letter: str | None = None) -> None:
    for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
        frame = data[system]
        x = frame.time_ns.to_numpy()
        ax.fill_between(x, frame.minimum.to_numpy(), frame.maximum.to_numpy(),
                        color=color, alpha=0.14, linewidth=0)
        ax.plot(x, frame["median"], color=color, lw=1.1, label=label)
    ax.set(xlim=xlim, xlabel="Time (ns)", ylabel=ylabel)
    if letter:
        panel(ax, letter)


def save_bundle(fig: plt.Figure, stem: str, width: float, height: float) -> None:
    fig.set_size_inches(width, height)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.04,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(PREVIEW / f"{stem}.png", dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# System-level metrics. Shading is the full range across three independent trajectories.
metric_specs = {
    "rmsd": ("whole_backbone_rmsd_0_100ns.xvg", 0, 100, 21),
    "distance": ("chainB_chainD_com_distance_pbc_nojump_0_100ns.xvg", 20, 100, 21),
    "rg": ("protein_rg_20_100ns.xvg", 20, 100, 21),
    "sasa": ("protein_sasa_20_100ns.xvg", 20, 100, 5),
    "hbond": ("interchain_hbonds_pbc_nojump_20_100ns.xvg", 20, 100, 5),
}
system_data = {}
for key, (filename, start, end, smooth) in metric_specs.items():
    system_data[key] = aggregate(filename, start, end, smooth)
    write_source(pd.concat(system_data[key].values(), ignore_index=True),
                 f"system_{key}_replicate_median_range.csv")


# Chain-local RMSF avoids artificial motion between protein chains.
residue = pd.read_csv(BASE / "summaries" / "per_residue_rmsf_and_contact_occupancy.csv")
rmsf_rows = []
for system in ("complex", "apo"):
    subset = residue[residue.system == system]
    for position, group in subset.groupby("protein_position", sort=True):
        rmsf_rows.append({
            "system": system,
            "protein_position": int(position),
            "source_chain": group.source_chain.iloc[0],
            "residue_number_in_chain": int(group.residue_number_in_chain.iloc[0]),
            "residue_name": group.residue_name.iloc[0],
            "median_nm": group.backbone_localfit_rmsf_nm.median(),
            "minimum_nm": group.backbone_localfit_rmsf_nm.min(),
            "maximum_nm": group.backbone_localfit_rmsf_nm.max(),
        })
rmsf = pd.DataFrame(rmsf_rows)
write_source(rmsf, "chain_local_backbone_rmsf_20_100ns.csv")


def draw_rmsf(ax: plt.Axes, letter: str | None = None) -> None:
    for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
        frame = rmsf[rmsf.system == system].sort_values("protein_position")
        x = frame.protein_position.to_numpy()
        ax.fill_between(x, frame.minimum_nm.to_numpy(), frame.maximum_nm.to_numpy(),
                        color=color, alpha=0.14, linewidth=0)
        ax.plot(x, frame.median_nm, color=color, lw=1.05, label=label)
    ax.axvline(93.5, color="#BDBDBD", lw=0.7, ls="--")
    ax.text(46.5, 0.97, "Chain B", transform=ax.get_xaxis_transform(), ha="center",
            va="top", fontsize=8, color=NEUTRAL)
    ax.text(139.5, 0.97, "Chain D", transform=ax.get_xaxis_transform(), ha="center",
            va="top", fontsize=8, color=NEUTRAL)
    ax.set(xlim=(1, 186), xlabel="Protein position", ylabel="Backbone RMSF (nm)")
    if letter:
        panel(ax, letter)


fig, axes = plt.subplots(2, 3, constrained_layout=True)
median_range(axes[0, 0], system_data["rmsd"], "Backbone RMSD (nm)", (0, 100), "a")
axes[0, 0].axvspan(0, 20, color="#D9D9D9", alpha=0.25, linewidth=0)
axes[0, 0].axvline(20, color="#999999", lw=0.6, ls="--")
median_range(axes[0, 1], system_data["distance"], "Chain B-D distance (nm)", (20, 100), "b")
draw_rmsf(axes[0, 2], "c")
median_range(axes[1, 0], system_data["rg"], r"$R_g$ (nm)", (20, 100), "d")
median_range(axes[1, 1], system_data["sasa"], r"SASA (nm$^2$)", (20, 100), "e")
median_range(axes[1, 2], system_data["hbond"], "Inter-chain H-bonds\n(0.5 ns mean)", (20, 100), "f")
fig.legend(handles=[Line2D([0], [0], color=COMPLEX, lw=1.5, label="Complex"),
                    Line2D([0], [0], color=APO, lw=1.5, label="Apo")],
           loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.015), handlelength=2.3)
save_bundle(fig, "Figure_1_MD_system_metrics", 7.20, 4.85)


# Direct replacements for the three legacy figures.
fig, ax = plt.subplots(constrained_layout=True)
median_range(ax, system_data["rmsd"], "Backbone RMSD (nm)", (0, 100))
ax.axvspan(0, 20, color="#D9D9D9", alpha=0.25, linewidth=0)
ax.axvline(20, color="#999999", lw=0.6, ls="--")
ax.legend(loc="upper left", ncol=2)
save_bundle(fig, "RMSD_latest", 3.50, 2.55)

fig, ax = plt.subplots(constrained_layout=True)
draw_rmsf(ax)
ax.legend(loc="upper center", ncol=2)
save_bundle(fig, "RMSF_latest", 3.50, 2.55)


# Ligand metrics distinguish unwrapped departure from periodic protein association.
ligand_specs = {
    "rmsd": ("ligand_rmsd_unwrapped_after_proteinfit_0_100ns.xvg", "Unwrapped ligand RMSD (nm)", 51),
    "pocket": ("ligand_original_pocket_com_distance_unwrapped_0_100ns.xvg", "Original-pocket displacement (nm)", 51),
    "mindist": ("ligand_protein_mindist_pbc_nojump_0_100ns.xvg", "Protein minimum distance (nm)", 101),
    "contacts": ("ligand_protein_contacts_pbc_nojump_0_100ns.xvg", "Protein contacts (1 ns mean)", 101),
    "hbond": ("ligand_protein_hbonds_pbc_nojump_20_100ns.xvg", "Ligand-protein H-bonds\n(0.5 ns mean)", 5),
}
ligand_data: dict[str, list[pd.DataFrame]] = {}
for key, (filename, _, smooth) in ligand_specs.items():
    start = 20 if key == "hbond" else 0
    frames = []
    for rep in range(1, 4):
        frame = load_series(run_path("complex", rep, filename), start, 100)
        frame["display_value"] = frame.value.rolling(
            smooth, center=True, min_periods=max(1, smooth // 4)).mean()
        frame.insert(0, "replicate", rep)
        frames.append(frame)
    ligand_data[key] = frames
    write_source(pd.concat(frames, ignore_index=True), f"ligand_{key}.csv")

contact = pd.read_csv(BASE / "summaries" / "ligand_contact_occupancy_replicate_consensus.csv")
contact = contact.sort_values(["occupancy_median", "occupancy_max"], ascending=False).head(14).copy()
contact["label"] = contact.source_chain + ":" + contact.residue_name + contact.residue_number_in_chain.astype(str)
write_source(contact, "ligand_contact_occupancy_top14.csv")

fig, axes = plt.subplots(2, 3, constrained_layout=True)
for ax, key, letter in zip(axes.flat[:5], ligand_specs, "abcde"):
    _, ylabel, _ = ligand_specs[key]
    for rep, frame in enumerate(ligand_data[key], start=1):
        ax.plot(frame.time_ns, frame.display_value, color=REP_COLORS[rep - 1],
                lw=0.82, alpha=0.94, label=f"Rep {rep}")
    ax.set(xlim=(20 if key == "hbond" else 0, 100), xlabel="Time (ns)", ylabel=ylabel)
    panel(ax, letter)
axes[0, 0].legend(loc="upper left", ncol=3, handlelength=1.4)

ax = axes[1, 2]
matrix = contact[["occupancy_rep1", "occupancy_rep2", "occupancy_rep3"]].to_numpy()
im = ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(3), ["Rep 1", "Rep 2", "Rep 3"])
ax.set_yticks(range(len(contact)), contact.label)
ax.tick_params(axis="both", length=0)
for row in range(matrix.shape[0]):
    for col in range(matrix.shape[1]):
        value = matrix[row, col]
        ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=8,
                color="white" if value > 0.55 else "#252525")
panel(ax, "f")
cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.025)
cbar.set_label("Contact occupancy", fontsize=8)
cbar.ax.tick_params(labelsize=8, length=2)
save_bundle(fig, "Figure_2_ligand_behavior", 7.20, 5.10)


fig, axes = plt.subplots(1, 2, constrained_layout=True)
median_range(axes[0], system_data["hbond"], "Inter-chain H-bonds\n(0.5 ns mean)", (20, 100), "a")
axes[0].legend(loc="upper center", ncol=2)
for rep, frame in enumerate(ligand_data["hbond"], start=1):
    axes[1].plot(frame.time_ns, frame.display_value, color=REP_COLORS[rep - 1],
                 lw=0.9, label=f"Rep {rep}")
axes[1].set(xlim=(20, 100), xlabel="Time (ns)", ylabel="Ligand-protein H-bonds\n(0.5 ns mean)")
axes[1].legend(loc="upper center", ncol=3, handlelength=1.5)
panel(axes[1], "b")
save_bundle(fig, "Hydrogen_bonds_latest", 7.20, 2.65)


# PCA and FEL from the corrected pooled C-alpha basis.
fel_root = BASE / "pooled_ca_pca_v3" / "fel_data"
projection = pd.read_csv(fel_root / "pc1_pc2_projection_all_replicates.csv")
variance = pd.read_csv(fel_root / "pooled_pca_explained_variance.csv")
fel = pd.read_csv(fel_root / "free_energy_landscape_common_grid_300K.csv")
fel_condition = fel[fel.group.isin(["complex_all_replicates", "apo_all_replicates"])].copy()
write_source(projection, "pooled_basis_pc1_pc2_projections.csv")
write_source(variance, "pooled_pca_explained_variance.csv")
write_source(fel_condition, "condition_free_energy_landscapes_300K.csv")
pc1 = variance.loc[variance["mode"] == 1, "variance_percent"].iloc[0]
pc2 = variance.loc[variance["mode"] == 2, "variance_percent"].iloc[0]

fig, axes = plt.subplots(1, 3, constrained_layout=True)
ax = axes[0]
for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
    frame = projection[projection.system == system]
    ax.scatter(frame.pc1_nm, frame.pc2_nm, s=4, color=color, alpha=0.14,
               linewidths=0, label=label)
ax.set(xlabel=f"PC1 ({pc1:.1f}%)", ylabel=f"PC2 ({pc2:.1f}%)")
ax.legend(loc="upper right", markerscale=2.3)
panel(ax, "a")

norm = Normalize(vmin=0, vmax=12)
mesh = None
for ax, group, letter, label in zip(
        axes[1:], ["complex_all_replicates", "apo_all_replicates"], "bc", ["Complex", "Apo"]):
    frame = fel_condition[fel_condition.group == group]
    x = np.sort(frame.pc1_center_nm.unique())
    y = np.sort(frame.pc2_center_nm.unique())
    grid = frame.pivot(index="pc2_center_nm", columns="pc1_center_nm",
                       values="free_energy_kj_mol").reindex(index=y, columns=x).to_numpy()
    mesh = ax.pcolormesh(x, y, np.ma.masked_invalid(grid), shading="nearest",
                         cmap="viridis_r", norm=norm)
    ax.text(0.04, 0.95, label, transform=ax.transAxes, va="top", fontsize=8)
    ax.set(xlabel=f"PC1 ({pc1:.1f}%)", ylabel=f"PC2 ({pc2:.1f}%)")
    panel(ax, letter)
cbar = fig.colorbar(mesh, ax=axes[1:], fraction=0.035, pad=0.025)
cbar.set_label(r"Free energy (kJ mol$^{-1}$)", fontsize=8)
cbar.ax.tick_params(labelsize=8, length=2)
save_bundle(fig, "Figure_3_PCA_FEL", 7.20, 2.65)


manifest = {
    "backend": "Python/matplotlib",
    "source": "comparative_analysis_v3 only; PBC cluster then nojump",
    "independent_unit": "trajectory replicate",
    "replicates_per_condition": 3,
    "production_length": "100 ns per trajectory; six trajectories",
    "analysis_window": "20-100 ns; selected time-series displayed from 0 ns",
    "summary": "replicate median with full replicate range; no frame-level inference",
    "ligand_departure": "unwrapped protein-fit coordinates",
    "ligand_association": "PBC-aware unrotated nojump coordinates",
    "PCA": f"pooled 186 C-alpha basis; 4806 frames; PC1={pc1:.3f}%; PC2={pc2:.3f}%",
    "exports": "editable SVG/PDF, 600 dpi LZW TIFF, PNG preview",
    "titles": "none; lowercase panel letters and axis labels only",
}
with (OUT / "FIGURE_TECHNICAL_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.writer(handle)
    writer.writerow(["field", "value"])
    writer.writerows(manifest.items())

print(OUT)
