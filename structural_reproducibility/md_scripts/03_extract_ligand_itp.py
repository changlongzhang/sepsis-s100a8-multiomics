from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
src = ROOT / "parameters" / "ligand" / "TPL.top"
dst = ROOT / "parameters" / "ligand" / "TPL.itp"
types_dst = ROOT / "parameters" / "ligand" / "TPL_atomtypes.itp"
lines = src.read_text(encoding="utf-8").splitlines()
start = next(i for i, x in enumerate(lines) if x.strip().lower() == "[ atomtypes ]")
end = next(i for i, x in enumerate(lines) if x.strip().lower() == "[ system ]")
mol_start = next(i for i, x in enumerate(lines) if x.strip().lower() == "[ moleculetype ]")
types_body = lines[start:mol_start]
mol_body = lines[mol_start:end]
types_dst.write_text("; OpenFF Sage 2.2.1 atom types for triptolide\n" + "\n".join(types_body) + "\n", encoding="utf-8")
dst.write_text("; OpenFF Sage 2.2.1 molecule parameters for triptolide\n" + "\n".join(mol_body) + "\n", encoding="utf-8")
if "[ moleculetype ]" not in dst.read_text(encoding="utf-8"):
    raise RuntimeError("Ligand moleculetype missing")
print(dst)
print(types_dst)
