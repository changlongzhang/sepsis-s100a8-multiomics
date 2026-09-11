from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PCA_FOLDER = sys.argv[1] if len(sys.argv) > 1 else "pooled_ca_pca_v3"
PCA = ROOT / "comparative_analysis_v3" / PCA_FOLDER
OUT = PCA / "fel_data"
OUT.mkdir(parents=True, exist_ok=True)
R_KJ_MOL_K = 0.008314462618
TEMPERATURE_K = 300.0
KBT = R_KJ_MOL_K * TEMPERATURE_K


def read_xvg_sets(path: Path) -> list[np.ndarray]:
    sets: list[list[list[float]]] = [[]]
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        if s == "&":
            if sets[-1]:
                sets.append([])
            continue
        try:
            sets[-1].append([float(x) for x in s.split()])
        except ValueError:
            continue
    return [np.asarray(x, dtype=float) for x in sets if x]


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


projections: dict[str, np.ndarray] = {}
projection_rows: list[dict] = []
for system in ("complex", "apo"):
    for rep in range(1, 4):
        tag = f"{system}_rep{rep}"
        sets = read_xvg_sets(PCA / f"{tag}_pc1_pc2_projection.xvg")
        if len(sets) != 2 or len(sets[0]) != len(sets[1]):
            raise ValueError(f"Expected two equal PC projection sets for {tag}")
        if not np.allclose(sets[0][:, 0], sets[1][:, 0]):
            raise ValueError(f"PC1/PC2 time mismatch for {tag}")
        data = np.column_stack((sets[0][:, 0] / 1000.0, sets[0][:, 1], sets[1][:, 1]))
        projections[tag] = data
        for time_ns, pc1, pc2 in data:
            projection_rows.append({
                "system": system,
                "replicate": rep,
                "time_ns": float(time_ns),
                "pc1_nm": float(pc1),
                "pc2_nm": float(pc2),
            })
write_csv(OUT / "pc1_pc2_projection_all_replicates.csv", projection_rows)

# The same bin edges are used for every replicate and condition so landscapes
# are directly comparable. Empty bins remain blank; they are not extrapolated.
all_xy = np.vstack([x[:, 1:3] for x in projections.values()])
bins = 50
x_edges = np.linspace(float(all_xy[:, 0].min()), float(all_xy[:, 0].max()), bins + 1)
y_edges = np.linspace(float(all_xy[:, 1].min()), float(all_xy[:, 1].max()), bins + 1)

groups: dict[str, np.ndarray] = dict(projections)
groups["complex_all_replicates"] = np.vstack([projections[f"complex_rep{i}"][:, 1:3] for i in range(1, 4)])
groups["apo_all_replicates"] = np.vstack([projections[f"apo_rep{i}"][:, 1:3] for i in range(1, 4)])
groups["all_six_replicates"] = all_xy

fel_rows: list[dict] = []
for label, raw in groups.items():
    xy = raw[:, 1:3] if raw.shape[1] == 3 else raw
    hist, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=(x_edges, y_edges))
    positive = hist > 0
    free = np.full(hist.shape, np.nan, dtype=float)
    free[positive] = -KBT * np.log(hist[positive] / np.max(hist))
    for ix in range(bins):
        for iy in range(bins):
            fel_rows.append({
                "group": label,
                "pc1_center_nm": float((x_edges[ix] + x_edges[ix + 1]) / 2),
                "pc2_center_nm": float((y_edges[iy] + y_edges[iy + 1]) / 2),
                "frame_count": int(hist[ix, iy]),
                "probability": float(hist[ix, iy] / len(xy)),
                "free_energy_kj_mol": ("" if not positive[ix, iy] else float(free[ix, iy])),
            })
write_csv(OUT / "free_energy_landscape_common_grid_300K.csv", fel_rows)

evals = read_xvg_sets(PCA / "pooled_ca_eigenvalues.xvg")[0]
values = evals[:, 1]
positive_sum = float(values[values > 0].sum())
variance_rows = []
cumulative = 0.0
for mode, value in evals:
    fraction = float(value / positive_sum) if value > 0 else 0.0
    cumulative += fraction
    variance_rows.append({
        "mode": int(mode),
        "eigenvalue_nm2": float(value),
        "variance_fraction": fraction,
        "variance_percent": fraction * 100.0,
        "cumulative_variance_percent": cumulative * 100.0,
    })
write_csv(OUT / "pooled_pca_explained_variance.csv", variance_rows)

(OUT / "PCA_FEL_COMPLETE.txt").write_text(
    f"completed; common pooled C-alpha basis; 186 atoms; six trajectories; 20-100 ns; 100 ps stride\n"
    f"FEL: 50x50 common grid; T={TEMPERATURE_K:.1f} K; kBT={KBT:.6f} kJ/mol; empty bins not extrapolated\n",
    encoding="utf-8",
)
print(f"PC1={variance_rows[0]['variance_percent']:.3f}% PC2={variance_rows[1]['variance_percent']:.3f}%")
