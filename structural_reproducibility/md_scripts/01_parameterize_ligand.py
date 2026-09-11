from pathlib import Path
import hashlib
import json
import numpy as np
import tempfile

from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT.parent / "S100A8_triptolide_docking_topjournal" / "results" / "representative_pose_triptolide.sdf"
OUT = ROOT / "parameters" / "ligand"
OUT.mkdir(parents=True, exist_ok=True)

source_bytes = INPUT.read_bytes()
with tempfile.TemporaryDirectory(prefix="s100a8_md_") as tmpdir:
    ascii_input = Path(tmpdir) / "triptolide_poses.sdf"
    ascii_input.write_bytes(source_bytes)
    molecules = Molecule.from_file(ascii_input)
if not isinstance(molecules, list) or len(molecules) != 20:
    raise RuntimeError(f"Expected 20 docking-pose records, obtained {type(molecules)} / {len(molecules)}")

smiles = {m.to_smiles(isomeric=True, explicit_hydrogens=True) for m in molecules}
if len(smiles) != 1:
    raise RuntimeError("Docking records do not have one identical stereochemical identity")

mol = molecules[0]
mol.name = "TPL"
for atom in mol.atoms:
    atom.metadata["residue_name"] = "TPL"
    atom.metadata["residue_number"] = 1

charge_model = "openff-gnn-am1bcc-1.0.0.pt"
mol.assign_partial_charges(charge_model)
charges = np.asarray(mol.partial_charges.m, dtype=np.float64)
if abs(charges.sum()) > 1e-8:
    raise RuntimeError(f"Non-integral charge sum: {charges.sum()}")

force_field_file = "openff-2.2.1.offxml"
ff = ForceField(force_field_file)
topology = Topology.from_molecules([mol])
interchange = Interchange.from_smirnoff(force_field=ff, topology=topology, charge_from_molecules=[mol])
# A temporary box is required only for GROMACS 2020+ file export. The solvated
# protein system receives its physical box later during system construction.
interchange.box = np.eye(3) * 5.0 * unit.nanometer
interchange.to_gromacs(str(OUT / "TPL"))
mol.to_file(OUT / "TPL_charged.sdf", file_format="sdf")

meta = {
    "source": str(INPUT),
    "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
    "selected_record": 1,
    "records_in_source": len(molecules),
    "unique_stereochemical_identities": len(smiles),
    "isomeric_explicit_h_smiles": next(iter(smiles)),
    "n_atoms": mol.n_atoms,
    "formal_charge_e": float(mol.total_charge.m),
    "partial_charge_model": charge_model,
    "partial_charge_sum_e": float(charges.sum()),
    "partial_charge_sha256_float64": hashlib.sha256(charges.tobytes()).hexdigest(),
    "force_field": force_field_file,
}
(OUT / "parameterization_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
print(json.dumps(meta, indent=2))
