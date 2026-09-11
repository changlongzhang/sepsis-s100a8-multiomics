from __future__ import annotations

from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "comparative_analysis_v1"
OUT = BASE / "publication_figures_python_v1"
SOURCE = OUT / "source_data"
PREVIEW = OUT / "preview"
for folder in (OUT, SOURCE, PREVIEW):
    folder.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.2,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

COMPLEX = "#B64342"
APO = "#3775BA"
REP_COLORS = ["#484878", "#42949E", "#B64342"]
NEUTRAL = "#767676"


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@" or s == "&":
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            pass
    if not rows:
        raise ValueError(f"No numeric rows in {path}")
    return np.asarray(rows, dtype=float)


def normalize_time(time: np.ndarray) -> np.ndarray:
    return time / 1000.0 if float(np.nanmax(time)) > 1000.0 else time


def load_series(path: Path, start: float | None = None, end: float | None = None) -> pd.DataFrame:
    arr = read_xvg(path)
    df = pd.DataFrame({"time_ns": normalize_time(arr[:, 0]), "value": arr[:, 1]})
    if start is not None:
        df = df[df.time_ns >= start]
    if end is not None:
        df = df[df.time_ns <= end + 1e-3]
    return df.reset_index(drop=True)


def analysis_path(system: str, rep: int, complex_name: str, apo_name: str | None = None) -> Path:
    name = complex_name if system == "complex" else (apo_name or complex_name)
    return ROOT / "systems" / system / f"rep{rep}" / "analysis" / name


def extended_path(system: str, rep: int, name: str) -> Path:
    return BASE / "per_trajectory" / f"{system}_rep{rep}" / name


def aggregate(paths: list[Path], start: float, end: float, smooth_points: int = 1) -> pd.DataFrame:
    frames = [load_series(p, start, end) for p in paths]
    time = frames[0].time_ns.to_numpy()
    values = []
    for frame in frames:
        if len(frame) != len(time) or not np.allclose(frame.time_ns, time, atol=1e-4):
            values.append(np.interp(time, frame.time_ns, frame.value))
        else:
            values.append(frame.value.to_numpy())
    mat = pd.DataFrame(np.vstack(values).T)
    if smooth_points > 1:
        mat = mat.rolling(smooth_points, center=True, min_periods=max(1, smooth_points // 4)).mean()
    return pd.DataFrame({
        "time_ns": time,
        "median": mat.median(axis=1).to_numpy(),
        "minimum": mat.min(axis=1).to_numpy(),
        "maximum": mat.max(axis=1).to_numpy(),
    })


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
    ax.tick_params(direction="out", length=2.5, width=0.6, pad=2)


def median_range_plot(ax: plt.Axes, complex_df: pd.DataFrame, apo_df: pd.DataFrame,
                      ylabel: str, xmin: float, xmax: float, letter: str | None = None) -> None:
    for df, color, label in ((complex_df, COMPLEX, "Complex"), (apo_df, APO, "Apo")):
        x = df.time_ns.to_numpy()
        ax.fill_between(x, df.minimum, df.maximum, color=color, alpha=0.14, linewidth=0)
        ax.plot(x, df["median"], color=color, lw=1.15, label=label)
    ax.set(xlim=(xmin, xmax), xlabel="Time (ns)", ylabel=ylabel)
    if letter:
        panel(ax, letter)


def save_bundle(fig: plt.Figure, stem: str, width_in: float, height_in: float) -> None:
    fig.set_size_inches(width_in, height_in)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.04,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(PREVIEW / f"{stem}.png", dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def write_source(df: pd.DataFrame, name: str) -> None:
    df.to_csv(SOURCE / name, index=False, encoding="utf-8-sig")


system_series: dict[str, dict[str, pd.DataFrame]] = {}
metric_specs = {
    "rmsd": (lambda s, r: analysis_path(s, r, "protein_backbone_rmsd.xvg", "rmsd_backbone.xvg"), 0.0, 100.0, 1),
    "com": (lambda s, r: analysis_path(s, r, "chainB_chainD_com_distance.xvg", "chainB_chainD_COM_distance.xvg"), 20.0, 100.0, 1),
    "rg": (lambda s, r: extended_path(s, r, "protein_rg_20_100ns.xvg"), 20.0, 100.0, 1),
    "sasa": (lambda s, r: extended_path(s, r, "protein_sasa_20_100ns.xvg"), 20.0, 100.0, 21),
    "hbond": (lambda s, r: extended_path(s, r, "interchain_hbonds_20_100ns.xvg"), 20.0, 100.0, 101),
}
for metric, (locator, start, end, smooth) in metric_specs.items():
    system_series[metric] = {}
    for system in ("complex", "apo"):
        df = aggregate([locator(system, rep) for rep in range(1, 4)], start, end, smooth)
        df.insert(0, "system", system)
        system_series[metric][system] = df
    write_source(pd.concat(system_series[metric].values(), ignore_index=True), f"system_{metric}_median_range.csv")


def local_rmsf_summary() -> pd.DataFrame:
    rows = []
    folder = BASE / "chain_local_rmsf"
    for system in ("complex", "apo"):
        for chain_i, chain in enumerate(("B", "D")):
            vals = []
            for rep in range(1, 4):
                arr = read_xvg(folder / f"{system}_rep{rep}_chain{chain}_backbone_localfit_rmsf_20_100ns.xvg")
                if len(arr) != 93:
                    raise ValueError(f"Expected 93 RMSF rows for {system} rep{rep} chain {chain}")
                vals.append(arr[:, 1])
            mat = np.vstack(vals)
            for i in range(93):
                rows.append({
                    "system": system,
                    "source_chain": chain,
                    "residue_number": i + 1,
                    "protein_position": chain_i * 93 + i + 1,
                    "median_nm": float(np.median(mat[:, i])),
                    "minimum_nm": float(np.min(mat[:, i])),
                    "maximum_nm": float(np.max(mat[:, i])),
                    "rep1_nm": float(mat[0, i]),
                    "rep2_nm": float(mat[1, i]),
                    "rep3_nm": float(mat[2, i]),
                })
    return pd.DataFrame(rows)


rmsf = local_rmsf_summary()
write_source(rmsf, "chain_local_backbone_rmsf_20_100ns.csv")


def draw_rmsf(ax: plt.Axes, letter: str | None = None) -> None:
    for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
        d = rmsf[rmsf.system == system].sort_values("protein_position")
        x = d.protein_position.to_numpy()
        ax.fill_between(x, d.minimum_nm, d.maximum_nm, color=color, alpha=0.14, linewidth=0)
        ax.plot(x, d.median_nm, color=color, lw=1.1, label=label)
    ax.axvline(93.5, color="#BDBDBD", lw=0.7, ls="--")
    ax.text(46.5, 0.97, "Chain B", transform=ax.get_xaxis_transform(), ha="center", va="top", color=NEUTRAL, fontsize=6)
    ax.text(139.5, 0.97, "Chain D", transform=ax.get_xaxis_transform(), ha="center", va="top", color=NEUTRAL, fontsize=6)
    ax.set(xlim=(1, 186), xlabel="Protein position", ylabel="Backbone RMSF (nm)")
    if letter:
        panel(ax, letter)


# Figure 1: system-level metrics
fig, axes = plt.subplots(2, 3, constrained_layout=True)
median_range_plot(axes[0, 0], system_series["rmsd"]["complex"], system_series["rmsd"]["apo"], "Backbone RMSD (nm)", 0, 100, "a")
axes[0, 0].axvspan(0, 20, color="#D9D9D9", alpha=0.25, lw=0)
axes[0, 0].axvline(20, color="#A0A0A0", lw=0.6, ls="--")
median_range_plot(axes[0, 1], system_series["com"]["complex"], system_series["com"]["apo"], "Chain B-D distance (nm)", 20, 100, "b")
draw_rmsf(axes[0, 2], "c")
median_range_plot(axes[1, 0], system_series["rg"]["complex"], system_series["rg"]["apo"], r"$R_g$ (nm)", 20, 100, "d")
median_range_plot(axes[1, 1], system_series["sasa"]["complex"], system_series["sasa"]["apo"], r"SASA (nm$^2$)", 20, 100, "e")
median_range_plot(axes[1, 2], system_series["hbond"]["complex"], system_series["hbond"]["apo"], "Inter-chain H-bonds\n(1 ns mean)", 20, 100, "f")
handles = [Line2D([0], [0], color=COMPLEX, lw=1.5, label="Complex"), Line2D([0], [0], color=APO, lw=1.5, label="Apo")]
fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.015), handlelength=2.3)
save_bundle(fig, "Figure_1_MD_system_metrics", 7.20, 4.85)


# Standalone replacements for the three legacy plots
fig, ax = plt.subplots(constrained_layout=True)
median_range_plot(ax, system_series["rmsd"]["complex"], system_series["rmsd"]["apo"], "Backbone RMSD (nm)", 0, 100)
ax.axvspan(0, 20, color="#D9D9D9", alpha=0.25, lw=0)
ax.axvline(20, color="#A0A0A0", lw=0.6, ls="--")
ax.legend(loc="upper left", ncol=2)
save_bundle(fig, "RMSD_latest", 3.50, 2.55)

fig, ax = plt.subplots(constrained_layout=True)
draw_rmsf(ax)
ax.legend(loc="upper center", ncol=2)
save_bundle(fig, "RMSF_latest", 3.50, 2.55)

fig, ax = plt.subplots(constrained_layout=True)
median_range_plot(ax, system_series["hbond"]["complex"], system_series["hbond"]["apo"], "Inter-chain H-bonds\n(1 ns mean)", 20, 100)
ax.legend(loc="upper center", ncol=2)
save_bundle(fig, "Hydrogen_bonds_latest", 3.50, 2.55)


# Figure 2: ligand behaviour across complex replicates
ligand_specs = {
    "rmsd": ("ligand_rmsd_after_proteinfit.xvg", "Ligand RMSD (nm)", 1),
    "pocket": ("ligand_pocket_com_distance.xvg", "Original-pocket distance (nm)", 1),
    "contacts": ("ligand_protein_contacts.xvg", "Protein contacts\n(1 ns mean)", 101),
}
ligand_data: dict[str, list[pd.DataFrame]] = {}
for key, (name, _, smooth) in ligand_specs.items():
    ligand_data[key] = []
    source_rows = []
    for rep in range(1, 4):
        d = load_series(analysis_path("complex", rep, name), 20, 100)
        if smooth > 1:
            d["display_value"] = d.value.rolling(smooth, center=True, min_periods=smooth // 4).mean()
        else:
            d["display_value"] = d.value
        d.insert(0, "replicate", rep)
        ligand_data[key].append(d)
        source_rows.append(d)
    write_source(pd.concat(source_rows, ignore_index=True), f"ligand_{key}_20_100ns.csv")

contact = pd.read_csv(BASE / "summaries" / "ligand_contact_occupancy_replicate_consensus.csv")
contact = contact.sort_values(["occupancy_median", "occupancy_max"], ascending=False).head(15).copy()
contact["label"] = contact.source_chain + ":" + contact.residue_name + contact.residue_number_in_chain.astype(str)
write_source(contact, "ligand_contact_occupancy_top15.csv")

fig, axes = plt.subplots(2, 2, constrained_layout=True)
for ax, key, letter in zip(axes.flat[:3], ("rmsd", "pocket", "contacts"), ("a", "b", "c")):
    ylabel = ligand_specs[key][1]
    for rep, d in enumerate(ligand_data[key], start=1):
        ax.plot(d.time_ns, d.display_value, color=REP_COLORS[rep - 1], lw=0.8, alpha=0.92, label=f"Rep {rep}")
    ax.set(xlim=(20, 100), xlabel="Time (ns)", ylabel=ylabel)
    panel(ax, letter)
axes[0, 0].legend(loc="upper center", ncol=3, handlelength=1.5)

ax = axes[1, 1]
mat = contact[["occupancy_rep1", "occupancy_rep2", "occupancy_rep3"]].to_numpy()
im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=max(0.65, float(np.nanmax(mat))))
ax.set_xticks(range(3), ["Rep 1", "Rep 2", "Rep 3"])
ax.set_yticks(range(len(contact)), contact.label)
ax.set_xlabel("Ligand contact occupancy")
ax.tick_params(axis="both", length=0)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                color="white" if mat[i, j] > 0.38 else "#272727")
panel(ax, "d")
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
cbar.set_label("Fraction of frames", fontsize=6)
cbar.ax.tick_params(labelsize=5.5, length=2)
save_bundle(fig, "Figure_2_ligand_behavior", 7.20, 5.35)


# Figure 3: pooled-basis PCA and condition-level FEL on the same grid
proj = pd.read_csv(BASE / "pooled_ca_pca_v2" / "fel_data" / "pc1_pc2_projection_all_replicates.csv")
write_source(proj, "pooled_basis_pc1_pc2_projections.csv")
fel = pd.read_csv(BASE / "pooled_ca_pca_v2" / "fel_data" / "free_energy_landscape_common_grid_300K.csv")
fel_sel = fel[fel.group.isin(["complex_all_replicates", "apo_all_replicates"])].copy()
write_source(fel_sel, "condition_free_energy_landscapes_300K.csv")

fig, axes = plt.subplots(1, 3, constrained_layout=True)
ax = axes[0]
for system, color, label in (("complex", COMPLEX, "Complex"), ("apo", APO, "Apo")):
    d = proj[proj.system == system]
    ax.scatter(d.pc1_nm, d.pc2_nm, s=4, color=color, alpha=0.12, linewidths=0, label=label)
ax.set(xlabel="PC1 (nm)", ylabel="PC2 (nm)")
ax.legend(loc="upper right", markerscale=2.5)
panel(ax, "a")

norm = Normalize(vmin=0, vmax=12)
mesh = None
for ax, group, letter, label in zip(axes[1:], ("complex_all_replicates", "apo_all_replicates"), ("b", "c"), ("Complex", "Apo")):
    d = fel_sel[fel_sel.group == group]
    x = np.sort(d.pc1_center_nm.unique())
    y = np.sort(d.pc2_center_nm.unique())
    grid = d.pivot(index="pc2_center_nm", columns="pc1_center_nm", values="free_energy_kj_mol").reindex(index=y, columns=x).to_numpy()
    masked = np.ma.masked_invalid(grid)
    mesh = ax.pcolormesh(x, y, masked, shading="nearest", cmap="viridis_r", norm=norm)
    ax.text(0.04, 0.95, label, transform=ax.transAxes, va="top", ha="left", fontsize=6.5, color="#272727")
    ax.set(xlabel="PC1 (nm)", ylabel="PC2 (nm)")
    panel(ax, letter)
cbar = fig.colorbar(mesh, ax=axes[1:], fraction=0.035, pad=0.025)
cbar.set_label(r"Free energy (kJ mol$^{-1}$)", fontsize=6.5)
cbar.ax.tick_params(labelsize=5.5, length=2)
save_bundle(fig, "Figure_3_PCA_FEL", 7.20, 2.65)


# Export a concise technical manifest.
manifest = {
    "backend": "Python/matplotlib",
    "independent_unit": "trajectory replicate",
    "replicates_per_condition": 3,
    "analysis_window_ns": "20-100 (RMSD display includes 0-100)",
    "time_trace_summary": "replicate median and full replicate range",
    "rmsf": "chain-backbone local fit; 93 residues per source chain",
    "fel": "accepted pooled_ca_pca_v2; common 50x50 grid; 300 K; empty bins unestimated",
    "exports": "SVG, PDF, 600 dpi LZW TIFF, 220 dpi PNG preview",
    "titles": "none; panel letters and axis labels only",
}
with (OUT / "FIGURE_TECHNICAL_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.writer(handle)
    writer.writerow(["field", "value"])
    writer.writerows(manifest.items())

print(OUT)
