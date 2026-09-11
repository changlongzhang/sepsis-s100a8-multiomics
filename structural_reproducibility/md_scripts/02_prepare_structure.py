from pathlib import Path
import hashlib
import json


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "S100A8_triptolide_docking_topjournal" / "results" / "S100A8_triptolide_representative_complex.pdb"
OUT = ROOT / "input" / "prepared"
OUT.mkdir(parents=True, exist_ok=True)

lines = SOURCE.read_text(encoding="ascii").splitlines()
protein, calcium, ligand = [], [], []
for line in lines:
    if not line.startswith(("ATOM  ", "HETATM")):
        continue
    rec, resn, chain = line[:6].strip(), line[17:20].strip(), line[21:22]
    if rec == "ATOM" and chain in {"B", "D"}:
        # pdb2gmx will rebuild hydrogens consistently with AMBER99SB-ILDN.
        if line[76:78].strip().upper() not in {"H", "D"}:
            protein.append(line)
    elif rec == "HETATM" and resn == "CA" and chain in {"B", "D"}:
        calcium.append(line)
    elif rec == "HETATM" and resn == "TPL":
        ligand.append(line)

if len(calcium) != 4 or len(ligand) != 50:
    raise RuntimeError(f"Unexpected retained heteroatoms: Ca={len(calcium)}, TPL={len(ligand)}")

def write_pdb(name, records):
    path = OUT / name
    path.write_text("\n".join(records + ["END"]) + "\n", encoding="ascii")
    return path

protein_path = write_pdb("S100A8_BD_protein.pdb", protein)
calcium_path = write_pdb("S100A8_BD_calcium.pdb", calcium)
ligand_path = write_pdb("triptolide_docked_pose1.pdb", ligand)
complex_path = write_pdb("S100A8_BD_TPL_Ca.pdb", protein + calcium + ligand)

meta = {
    "source": str(SOURCE),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "protein_heavy_atoms": len(protein),
    "retained_calcium_ions": len(calcium),
    "ligand_atoms": len(ligand),
    "retained_chains": ["B", "D"],
    "removed": ["all hydrogens for force-field rebuild", "Zn2+", "crystallization additives", "water"],
}
(OUT / "structure_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
print(json.dumps(meta, indent=2))

