from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "comparative_analysis_v1" / "summaries"
OUT.mkdir(parents=True, exist_ok=True)


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            s = line.strip()
            if not s or s[0] in "#@":
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric data: {path}")
    return np.asarray(rows, dtype=float)


def sample_sd(x: np.ndarray) -> float:
    return float(np.std(x, ddof=1)) if len(x) > 1 else math.nan


def describe(x: np.ndarray) -> dict[str, float | int]:
    return {
        "n": int(len(x)),
        "mean": float(np.mean(x)),
        "sd": sample_sd(x),
        "median": float(np.median(x)),
        "q25": float(np.quantile(x, 0.25)),
        "q75": float(np.quantile(x, 0.75)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "last": float(x[-1]),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


common = {
    "whole_backbone_rmsd_nm": ("analysis", "protein_backbone_rmsd.xvg", "rmsd_backbone.xvg"),
    "chain_B_backbone_rmsd_nm": ("analysis", "chainB_backbone_rmsd.xvg", "rmsd_chainB_backbone.xvg"),
    "chain_D_backbone_rmsd_nm": ("analysis", "chainD_backbone_rmsd.xvg", "rmsd_chainD_backbone.xvg"),
    "chain_B_D_com_distance_nm": ("analysis", "chainB_chainD_com_distance.xvg", "chainB_chainD_COM_distance.xvg"),
    "protein_rg_nm": ("extended", "protein_rg_20_100ns.xvg", "protein_rg_20_100ns.xvg"),
    "protein_sasa_nm2": ("extended", "protein_sasa_20_100ns.xvg", "protein_sasa_20_100ns.xvg"),
    "interchain_hbonds_count": ("extended", "interchain_hbonds_20_100ns.xvg", "interchain_hbonds_20_100ns.xvg"),
}

complex_only = {
    "ligand_rmsd_after_proteinfit_nm": "ligand_rmsd_after_proteinfit.xvg",
    "ligand_original_pocket_com_distance_nm": "ligand_pocket_com_distance.xvg",
    "ligand_protein_min_distance_nm": "ligand_protein_mindist.xvg",
    "ligand_protein_contacts_count": "ligand_protein_contacts.xvg",
}

replicate_rows: list[dict] = []
block_rows: list[dict] = []

for system in ("complex", "apo"):
    for rep in range(1, 4):
        analysis = ROOT / "systems" / system / f"rep{rep}" / "analysis"
        extended = ROOT / "comparative_analysis_v1" / "per_trajectory" / f"{system}_rep{rep}"
        metric_paths: dict[str, Path] = {}
        for metric, (source, complex_name, apo_name) in common.items():
            metric_paths[metric] = (analysis if source == "analysis" else extended) / (complex_name if system == "complex" else apo_name)
        if system == "complex":
            metric_paths.update({m: analysis / f for m, f in complex_only.items()})
            metric_paths["ligand_protein_hbonds_count"] = extended / "ligand_protein_hbonds_20_100ns.xvg"

        for metric, path in metric_paths.items():
            arr = read_xvg(path)
            raw_time = arr[:, 0]
            # Earlier quick-gate files were exported in ns, whereas the later
            # extended metrics retain GROMACS' default ps time unit.
            time_ns = raw_time / 1000.0 if float(np.max(raw_time)) > 1000.0 else raw_time
            values = arr[:, 1]
            keep = (time_ns >= 20.0) & (time_ns <= 100.001)
            time_ns, values = time_ns[keep], values[keep]
            if len(values) == 0:
                raise ValueError(f"No 20-100 ns observations for {metric}: {path}")
            row = {"system": system, "replicate": rep, "metric": metric, **describe(values)}
            replicate_rows.append(row)
            for start in (20, 40, 60, 80):
                end = start + 20
                mask = (time_ns >= start) & ((time_ns < end) if end < 100 else (time_ns <= end + 0.001))
                b = values[mask]
                block_rows.append({
                    "system": system,
                    "replicate": rep,
                    "metric": metric,
                    "block_start_ns": float(start),
                    "block_end_ns": float(end),
                    **describe(b),
                })

write_csv(OUT / "per_replicate_metric_summary_20_100ns.csv", replicate_rows)
write_csv(OUT / "per_replicate_20ns_block_summary.csv", block_rows)

# Replicate-level consistency. Replicate means, not trajectory frames, are the
# independent units; no frame-level inferential P values are calculated.
consistency_rows: list[dict] = []
for system in ("complex", "apo"):
    metrics = sorted({r["metric"] for r in replicate_rows if r["system"] == system})
    for metric in metrics:
        means = np.array([r["mean"] for r in replicate_rows if r["system"] == system and r["metric"] == metric], dtype=float)
        med = float(np.median(means))
        consistency_rows.append({
            "system": system,
            "metric": metric,
            "n_replicates": len(means),
            "replicate_mean_median": med,
            "replicate_mean_min": float(np.min(means)),
            "replicate_mean_max": float(np.max(means)),
            "replicate_mean_range": float(np.ptp(means)),
            "replicate_mean_sd": sample_sd(means),
            "replicate_mean_cv_percent": (sample_sd(means) / abs(float(np.mean(means))) * 100.0) if float(np.mean(means)) != 0 else math.nan,
        })
write_csv(OUT / "replicate_consistency_summary.csv", consistency_rows)

# Within-trajectory late-block stability diagnostic.
late_rows: list[dict] = []
for system in ("complex", "apo"):
    for rep in range(1, 4):
        metrics = sorted({r["metric"] for r in block_rows if r["system"] == system and r["replicate"] == rep})
        for metric in metrics:
            selected = [r for r in block_rows if r["system"] == system and r["replicate"] == rep and r["metric"] == metric]
            b60 = next(r for r in selected if r["block_start_ns"] == 60.0)
            b80 = next(r for r in selected if r["block_start_ns"] == 80.0)
            delta = float(b80["mean"] - b60["mean"])
            late_rows.append({
                "system": system,
                "replicate": rep,
                "metric": metric,
                "mean_60_80ns": b60["mean"],
                "mean_80_100ns": b80["mean"],
                "late_block_delta": delta,
                "absolute_late_block_delta": abs(delta),
            })
write_csv(OUT / "late_block_stability_diagnostic.csv", late_rows)

# Residue map comes from ordered CA records. GROMACS renumbered topology chains
# A/B correspond to the source chain groups B/D used throughout the audit.
pdb = ROOT / "comparative_analysis_v1" / "per_trajectory" / "complex_rep1" / "ca_average_20_100ns.pdb"
residues = []
for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
    if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "CA":
        residues.append({
            "protein_position": len(residues) + 1,
            "topology_chain": line[21:22].strip(),
            "source_chain": "B" if len(residues) < 93 else "D",
            "residue_number_in_chain": int(line[22:26]),
            "residue_name": line[17:20].strip(),
        })
if len(residues) != 186:
    raise ValueError(f"Expected 186 CA residues, found {len(residues)}")

residue_rows: list[dict] = []
for system in ("complex", "apo"):
    for rep in range(1, 4):
        rmsf = read_xvg(ROOT / "comparative_analysis_v1" / "per_trajectory" / f"{system}_rep{rep}" / "ca_rmsf_20_100ns.xvg")
        if len(rmsf) != 186:
            raise ValueError(f"Unexpected RMSF length for {system} rep{rep}: {len(rmsf)}")
        contact = None
        if system == "complex":
            contact = read_xvg(ROOT / "comparative_analysis_v1" / "contact_occupancy_v2" / f"complex_rep{rep}_protein_residue_contact_occupancy_0p4nm_20_100ns.xvg")
            if len(contact) != 186:
                raise ValueError(f"Unexpected contact occupancy length for complex rep{rep}: {len(contact)}")
        for i, residue in enumerate(residues):
            residue_rows.append({
                "system": system,
                "replicate": rep,
                **residue,
                "ca_rmsf_nm": float(rmsf[i, 1]),
                "ligand_contact_occupancy_0p4nm": (float(contact[i, 1]) if contact is not None else ""),
            })
write_csv(OUT / "per_residue_rmsf_and_contact_occupancy.csv", residue_rows)

contact_rows: list[dict] = []
for i, residue in enumerate(residues):
    values = np.array([
        float(r["ligand_contact_occupancy_0p4nm"])
        for r in residue_rows
        if r["system"] == "complex" and r["protein_position"] == i + 1
    ])
    contact_rows.append({
        **residue,
        "occupancy_rep1": values[0],
        "occupancy_rep2": values[1],
        "occupancy_rep3": values[2],
        "occupancy_median": float(np.median(values)),
        "occupancy_mean": float(np.mean(values)),
        "occupancy_min": float(np.min(values)),
        "occupancy_max": float(np.max(values)),
        "replicates_ge_0p10": int(np.sum(values >= 0.10)),
        "replicates_ge_0p25": int(np.sum(values >= 0.25)),
    })
write_csv(OUT / "ligand_contact_occupancy_replicate_consensus.csv", contact_rows)

manifest = OUT / "SUMMARY_COMPLETE.txt"
manifest.write_text(
    "Replicate-aware summaries completed. Independent unit: trajectory replicate (n=3 per condition).\n"
    "Window: 20-100 ns; block diagnostics: 20 ns; contact cutoff: 0.4 nm.\n"
    "No frame-level inferential tests were performed.\n",
    encoding="utf-8",
)
print(str(manifest))
