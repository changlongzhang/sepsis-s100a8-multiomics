from __future__ import annotations

from io import BytesIO
from pathlib import Path
import csv
import json
import re

from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
OUT = ROOT / "results" / "comparative_analysis"
OUT.mkdir(parents=True, exist_ok=True)
OH_LOCAL = ROOT / "results" / "local_replicates"
TPL_LOCAL = TPL_ROOT / "results" / "local_replicates"
TPL_PREP = TPL_ROOT / "input" / "prepared"


def first_pose(path: Path) -> dict:
    coordinates, names, types, serials = [], [], [], []
    score = None
    input_to_serial = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            score = float(line.split()[3])
        if line.startswith("REMARK SMILES IDX"):
            values = [int(value) for value in line.split()[3:]]
            for index in range(0, len(values), 2):
                input_to_serial[values[index]] = values[index + 1]
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line.split()[-1]
            serial = int(line[6:11])
            if not atom_type.startswith("H"):
                serials.append(serial)
                names.append(line[12:16].strip())
                types.append(atom_type)
                coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {
        "coords": np.asarray(coordinates), "names": names, "types": types,
        "serials": serials, "score": score, "input_to_serial": input_to_serial,
        "serial_to_coord": {serial: coord for serial, coord in zip(serials, coordinates)},
    }


def symmetry_rmsd(pose: dict, reference: dict) -> float:
    total, count = 0.0, 0
    for atom_type in sorted(set(pose["types"])):
        left = [index for index, value in enumerate(pose["types"]) if value == atom_type]
        right = [index for index, value in enumerate(reference["types"]) if value == atom_type]
        if len(left) != len(right):
            continue
        costs = cdist(pose["coords"][left], reference["coords"][right]) ** 2
        rows, columns = linear_sum_assignment(costs)
        total += costs[rows, columns].sum()
        count += len(rows)
    return float(np.sqrt(total / count))


def contacts(pose: dict, receptor: str) -> list[dict]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(receptor, str(TPL_PREP / f"5HLO_S100A8_dimer_{receptor}_pH7p4.pdb"))
    rows = []
    for residue in structure.get_residues():
        if residue.id[0] != " ":
            continue
        distances = []
        for atom in residue.get_atoms():
            if atom.element in {"H", "D"}:
                continue
            distance = float(np.linalg.norm(pose["coords"] - atom.coord, axis=1).min())
            distances.append((distance, atom.name))
        distance, atom_name = min(distances)
        if distance <= 4.5:
            rows.append({
                "chain": residue.parent.id, "residue": residue.resname, "number": residue.id[1],
                "min_distance_A": round(distance, 3), "nearest_atom": atom_name,
            })
    return rows


def load_molecule(path: Path) -> Chem.Mol:
    molecule = next(iter(Chem.ForwardSDMolSupplier(BytesIO(path.read_bytes()), removeHs=False)))
    if molecule is None:
        raise RuntimeError(f"Could not parse {path}")
    return Chem.RemoveHs(molecule)


def direct_common_scaffold_rmsd(tpl_pose: dict, oh_pose: dict, tpl_mol: Chem.Mol, oh_mol: Chem.Mol) -> tuple[float, int]:
    mcs = rdFMCS.FindMCS(
        [tpl_mol, oh_mol], atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact, ringMatchesRingOnly=True,
        completeRingsOnly=True, matchChiralTag=False, timeout=60,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    tpl_matches = tpl_mol.GetSubstructMatches(query, uniquify=False)
    oh_matches = oh_mol.GetSubstructMatches(query, uniquify=False)
    best = float("inf")
    used = 0
    for tpl_match in tpl_matches:
        for oh_match in oh_matches:
            left, right = [], []
            for tpl_index, oh_index in zip(tpl_match, oh_match):
                tpl_serial = tpl_pose["input_to_serial"].get(tpl_index + 1)
                oh_serial = oh_pose["input_to_serial"].get(oh_index + 1)
                if tpl_serial in tpl_pose["serial_to_coord"] and oh_serial in oh_pose["serial_to_coord"]:
                    left.append(tpl_pose["serial_to_coord"][tpl_serial])
                    right.append(oh_pose["serial_to_coord"][oh_serial])
            if len(left) < 10:
                continue
            value = float(np.sqrt(np.mean(np.sum((np.asarray(left) - np.asarray(right)) ** 2, axis=1))))
            if value < best:
                best, used = value, len(left)
    if not np.isfinite(best):
        raise RuntimeError("Could not map common scaffold coordinates")
    return best, used


oh_summary = pd.read_csv(OH_LOCAL / "local_replicates_summary.csv")
tpl_summary = pd.read_csv(TPL_LOCAL / "local_replicates_summary.csv")
oh_poses = {row.pose_file: first_pose(OH_LOCAL / row.pose_file) for _, row in oh_summary.iterrows()}
tpl_poses = {row.pose_file: first_pose(TPL_LOCAL / row.pose_file) for _, row in tpl_summary.iterrows()}

for label, summary, poses in (("16-hydroxytriptolide", oh_summary, oh_poses), ("triptolide", tpl_summary, tpl_poses)):
    rmsd_values = []
    for _, row in summary.iterrows():
        reference = poses[f"{row.receptor}_{row.scoring}_seed107985.pdbqt"]
        rmsd_values.append(symmetry_rmsd(poses[row.pose_file], reference))
    summary["rmsd_to_condition_reference_A"] = rmsd_values
    summary["compound"] = label

combined = pd.concat([tpl_summary, oh_summary], ignore_index=True)
combined.to_csv(OUT / "matched_local_replicates_with_rmsd.csv", index=False, encoding="utf-8-sig")
statistics = combined.groupby(["compound", "receptor", "scoring"]).agg(
    n=("seed", "count"), mean_score=("best_score_kcal_mol", "mean"),
    sd_score=("best_score_kcal_mol", "std"), min_score=("best_score_kcal_mol", "min"),
    max_score=("best_score_kcal_mol", "max"),
    median_rmsd_A=("rmsd_to_condition_reference_A", "median"),
    max_rmsd_A=("rmsd_to_condition_reference_A", "max"),
).reset_index()
statistics.to_csv(OUT / "matched_replicate_statistics.csv", index=False, encoding="utf-8-sig")

tpl_molecule = load_molecule(TPL_ROOT / "input" / "prepared" / "triptolide_primary.sdf")
oh_molecule = load_molecule(ROOT / "input" / "prepared" / "16-hydroxytriptolide_primary.sdf")
comparison_rows = []
contact_rows = []
for receptor in ("AC", "BD"):
    for scoring in ("vina", "vinardo"):
        tpl_pose = tpl_poses[f"{receptor}_{scoring}_seed107985.pdbqt"]
        oh_pose = oh_poses[f"{receptor}_{scoring}_seed107985.pdbqt"]
        tpl_contacts = contacts(tpl_pose, receptor)
        oh_contacts = contacts(oh_pose, receptor)
        tpl_set = {f"{row['chain']}:{row['residue']}{row['number']}" for row in tpl_contacts}
        oh_set = {f"{row['chain']}:{row['residue']}{row['number']}" for row in oh_contacts}
        union = tpl_set | oh_set
        intersection = tpl_set & oh_set
        scaffold_rmsd, scaffold_atoms = direct_common_scaffold_rmsd(
            tpl_pose, oh_pose, tpl_molecule, oh_molecule)
        comparison_rows.append({
            "receptor": receptor, "scoring": scoring,
            "triptolide_score": tpl_pose["score"], "16oh_score": oh_pose["score"],
            "score_difference_16oh_minus_tpl": oh_pose["score"] - tpl_pose["score"],
            "pose_centroid_distance_A": float(np.linalg.norm(
                tpl_pose["coords"].mean(axis=0) - oh_pose["coords"].mean(axis=0))),
            "common_scaffold_direct_rmsd_A": scaffold_rmsd,
            "common_scaffold_atoms_mapped": scaffold_atoms,
            "triptolide_contact_count": len(tpl_set), "16oh_contact_count": len(oh_set),
            "shared_contact_count": len(intersection),
            "contact_jaccard": len(intersection) / len(union) if union else np.nan,
            "shared_contacts": ";".join(sorted(intersection)),
            "triptolide_only_contacts": ";".join(sorted(tpl_set - oh_set)),
            "16oh_only_contacts": ";".join(sorted(oh_set - tpl_set)),
        })
        for compound, rows in (("triptolide", tpl_contacts), ("16-hydroxytriptolide", oh_contacts)):
            for row in rows:
                contact_rows.append({"compound": compound, "receptor": receptor, "scoring": scoring, **row})

pd.DataFrame(comparison_rows).to_csv(OUT / "representative_pose_cross_compound_comparison.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(contact_rows).to_csv(OUT / "representative_pose_contacts_4p5A.csv", index=False, encoding="utf-8-sig")

tpl_blind = pd.read_csv(TPL_ROOT / "results" / "blind_tiling" / "blind_tiling_summary.csv")
oh_blind = pd.read_csv(ROOT / "results" / "blind_tiling" / "blind_tiling_summary.csv")
tpl_blind["compound"] = "triptolide"
oh_blind["compound"] = "16-hydroxytriptolide"
blind = pd.concat([tpl_blind, oh_blind], ignore_index=True)
blind["rank_within_compound_receptor"] = blind.groupby(["compound", "receptor"])["best_score_kcal_mol"].rank(method="min")
blind.to_csv(OUT / "matched_blind_tiling_comparison.csv", index=False, encoding="utf-8-sig")

summary = {
    "statistics": statistics.to_dict(orient="records"),
    "representative_pose_comparison": comparison_rows,
    "blind_global_best": blind.sort_values("best_score_kcal_mol").groupby(["compound", "receptor"]).first().reset_index().to_dict(orient="records"),
}
(OUT / "matched_docking_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
(OUT / "MATCHED_DOCKING_ANALYSIS_COMPLETE.txt").write_text(
    "completed; matched protocol; blind tiling plus 2 receptors x 2 scoring functions x 10 seeds; no 16OH MD or experiment\n",
    encoding="utf-8",
)
print(statistics.to_string(index=False))
print(pd.DataFrame(comparison_rows).to_string(index=False))
