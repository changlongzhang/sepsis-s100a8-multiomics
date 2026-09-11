from __future__ import annotations

from io import BytesIO, StringIO
from pathlib import Path
import hashlib
import json
import subprocess

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, MACCSkeys, rdFMCS, rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
PYTHON = Path(r"C:\Users\ZCL\miniconda3\envs\s100a8_docking\python.exe")
MEEKO = Path(r"C:\Users\ZCL\miniconda3\envs\s100a8_docking\Scripts\mk_prepare_ligand.py")
SOURCE = ROOT / "input" / "downloaded" / "16-hydroxytriptolide_CID126556.sdf"
TPL_SOURCE = TPL_ROOT / "input" / "downloaded" / "triptolide_CID107985.sdf"
OUT = ROOT / "input" / "prepared"
ANALYSIS = ROOT / "results" / "comparative_analysis"
OUT.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)


def load_first(path: Path) -> Chem.Mol:
    molecule = next(iter(Chem.ForwardSDMolSupplier(BytesIO(path.read_bytes()), removeHs=False)))
    if molecule is None:
        raise RuntimeError(f"Could not parse {path}")
    molecule = Chem.RemoveHs(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule


def descriptors(molecule: Chem.Mol) -> dict:
    return {
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": Descriptors.MolWt(molecule),
        "heavy_atoms": molecule.GetNumHeavyAtoms(),
        "h_bond_donors": rdMolDescriptors.CalcNumHBD(molecule),
        "h_bond_acceptors": rdMolDescriptors.CalcNumHBA(molecule),
        "tpsa_A2": rdMolDescriptors.CalcTPSA(molecule),
        "clogp": Crippen.MolLogP(molecule),
        "formal_charge": Chem.GetFormalCharge(molecule),
    }


tpl = load_first(TPL_SOURCE)
oh16 = load_first(SOURCE)
tpl_smiles = Chem.MolToSmiles(tpl, isomericSmiles=True)
oh16_smiles = Chem.MolToSmiles(oh16, isomericSmiles=True)

generator = GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
morgan_similarity = DataStructs.TanimotoSimilarity(generator.GetFingerprint(tpl), generator.GetFingerprint(oh16))
generator_no_chirality = GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
morgan_similarity_no_chirality = DataStructs.TanimotoSimilarity(
    generator_no_chirality.GetFingerprint(tpl), generator_no_chirality.GetFingerprint(oh16)
)
rdkit_similarity = DataStructs.TanimotoSimilarity(Chem.RDKFingerprint(tpl), Chem.RDKFingerprint(oh16))
maccs_similarity = DataStructs.TanimotoSimilarity(
    MACCSkeys.GenMACCSKeys(tpl), MACCSkeys.GenMACCSKeys(oh16)
)
mcs = rdFMCS.FindMCS(
    [tpl, oh16],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareOrderExact,
    ringMatchesRingOnly=True,
    completeRingsOnly=True,
    matchChiralTag=True,
    timeout=60,
)
if mcs.canceled:
    raise RuntimeError("MCS search timed out")
mcs_atoms = int(mcs.numAtoms)
mcs_relaxed = rdFMCS.FindMCS(
    [tpl, oh16],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareOrderExact,
    ringMatchesRingOnly=True,
    completeRingsOnly=True,
    matchChiralTag=False,
    timeout=60,
)
if mcs_relaxed.canceled:
    raise RuntimeError("Relaxed MCS search timed out")
mcs_relaxed_atoms = int(mcs_relaxed.numAtoms)

# Rebuild the PubChem-defined stereoisomer using the same deterministic ETKDGv3/MMFF94s workflow as TPL.
base = Chem.AddHs(Chem.MolFromSmiles(oh16_smiles))
params = AllChem.ETKDGv3()
params.randomSeed = 126556
params.useSmallRingTorsions = True
conformer_ids = list(AllChem.EmbedMultipleConfs(base, numConfs=20, params=params))
if not conformer_ids:
    raise RuntimeError("No 16-hydroxytriptolide conformers generated")
properties = AllChem.MMFFGetMoleculeProperties(base, mmffVariant="MMFF94s")
energies = []
for conformer_id in conformer_ids:
    AllChem.MMFFOptimizeMolecule(base, mmffVariant="MMFF94s", confId=conformer_id, maxIters=2000)
    forcefield = AllChem.MMFFGetMoleculeForceField(base, properties, confId=conformer_id)
    energies.append((conformer_id, float(forcefield.CalcEnergy())))
energies.sort(key=lambda item: item[1])

all_sdf = StringIO()
writer = Chem.SDWriter(all_sdf)
for rank, (conformer_id, energy) in enumerate(energies, start=1):
    molecule = Chem.Mol(base)
    molecule.RemoveAllConformers()
    molecule.AddConformer(base.GetConformer(conformer_id), assignId=True)
    molecule.SetProp("_Name", f"16-hydroxytriptolide_conf_{rank:02d}")
    molecule.SetProp("MMFF94s_energy_kcal_mol", f"{energy:.6f}")
    writer.write(molecule)
writer.close()
(OUT / "16-hydroxytriptolide_conformers_MMFF94s.sdf").write_text(all_sdf.getvalue(), encoding="utf-8")

primary = Chem.Mol(base)
primary.RemoveAllConformers()
primary.AddConformer(base.GetConformer(energies[0][0]), assignId=True)
primary.SetProp("_Name", "16-hydroxytriptolide_CID126556_primary")
primary_sdf = StringIO()
primary_writer = Chem.SDWriter(primary_sdf)
primary_writer.write(primary)
primary_writer.close()
primary_path = OUT / "16-hydroxytriptolide_primary.sdf"
primary_path.write_text(primary_sdf.getvalue(), encoding="utf-8")

pdbqt_path = OUT / "16-hydroxytriptolide_primary.pdbqt"
command = [
    str(PYTHON), str(MEEKO), "-i", primary_path.name, "-o", pdbqt_path.name,
    "--charge_model", "gasteiger", "--add_index_map", "--rename_atoms",
]
process = subprocess.run(
    command, cwd=OUT, capture_output=True, text=True, encoding="utf-8", errors="replace"
)
(ROOT / "audit" / "meeko_ligand_preparation.log").write_text(
    "COMMAND\n" + subprocess.list2cmdline(command) + "\n\nSTDOUT\n" + process.stdout +
    "\nSTDERR\n" + process.stderr,
    encoding="utf-8",
)
if process.returncode != 0 or not pdbqt_path.exists():
    raise RuntimeError(f"Meeko failed: {process.stderr}")

source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
metadata = {
    "comparison_purpose": "Computational bridge only; no MD or experimental validation for 16-hydroxytriptolide",
    "triptolide": {
        "name": "triptolide", "pubchem_cid": 107985,
        "isomeric_smiles": tpl_smiles, "inchikey": Chem.MolToInchiKey(tpl), **descriptors(tpl),
    },
    "16-hydroxytriptolide": {
        "name": "16-hydroxytriptolide", "pubchem_cid": 126556,
        "isomeric_smiles": oh16_smiles, "inchikey": Chem.MolToInchiKey(oh16),
        "source_sha256": source_hash, **descriptors(oh16),
        "conformer_method": "RDKit ETKDGv3; fixed PubChem stereochemistry; 20 conformers; MMFF94s",
        "conformer_seed": 126556,
        "primary_MMFF94s_energy_kcal_mol": energies[0][1],
        "pdbqt": "Meeko 0.7.1; Gasteiger charges",
    },
    "similarity": {
        "morgan_radius2_2048_tanimoto_chiral": morgan_similarity,
        "morgan_radius2_2048_tanimoto_no_chirality": morgan_similarity_no_chirality,
        "rdkit_topological_tanimoto": rdkit_similarity,
        "maccs_tanimoto": maccs_similarity,
        "mcs_smarts": mcs.smartsString,
        "mcs_heavy_atoms": mcs_atoms,
        "mcs_fraction_of_triptolide_heavy_atoms": mcs_atoms / tpl.GetNumHeavyAtoms(),
        "mcs_fraction_of_16oh_heavy_atoms": mcs_atoms / oh16.GetNumHeavyAtoms(),
        "relaxed_mcs_smarts": mcs_relaxed.smartsString,
        "relaxed_mcs_heavy_atoms": mcs_relaxed_atoms,
        "relaxed_mcs_fraction_of_triptolide_heavy_atoms": mcs_relaxed_atoms / tpl.GetNumHeavyAtoms(),
        "relaxed_mcs_fraction_of_16oh_heavy_atoms": mcs_relaxed_atoms / oh16.GetNumHeavyAtoms(),
    },
}
(OUT / "compound_comparison_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

with (ANALYSIS / "compound_similarity.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    import csv
    writer = csv.writer(handle)
    writer.writerow(["metric", "value"])
    for key, value in metadata["similarity"].items():
        writer.writerow([key, value])

print(json.dumps(metadata, indent=2))
