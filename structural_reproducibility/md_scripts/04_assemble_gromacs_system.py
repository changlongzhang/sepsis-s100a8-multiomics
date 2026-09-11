from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "systems" / "complex" / "build"
LIG = ROOT / "parameters" / "ligand"


def gro_atoms(path):
    lines = path.read_text(encoding="ascii").splitlines()
    n = int(lines[1].strip())
    return lines[2 : 2 + n]


protein = gro_atoms(BUILD / "protein.gro")
calcium = gro_atoms(BUILD / "calcium.gro")
ligand = gro_atoms(LIG / "TPL.gro")
if len(calcium) != 4 or len(ligand) != 50:
    raise RuntimeError((len(calcium), len(ligand)))

# Renumber atom serials only; coordinates remain in the common docking frame.
combined = []
for i, line in enumerate(protein + calcium + ligand, 1):
    combined.append(line[:15] + f"{i:5d}" + line[20:])
(BUILD / "complex_unsolvated.gro").write_text(
    "S100A8 homodimer + 4 Ca2+ + triptolide\n"
    + f"{len(combined):5d}\n"
    + "\n".join(combined)
    + "\n   0.00000   0.00000   0.00000\n",
    encoding="ascii",
)

shutil.copy2(LIG / "TPL.itp", BUILD / "TPL.itp")
shutil.copy2(LIG / "TPL_atomtypes.itp", BUILD / "TPL_atomtypes.itp")
top = (BUILD / "topol.top").read_text(encoding="utf-8")
ffmarker = '#include "amber99sb-ildn.ff/forcefield.itp"'
top = top.replace(ffmarker, ffmarker + '\n\n; OpenFF ligand atom types must precede all molecule definitions\n#include "TPL_atomtypes.itp"')
marker = '; Include water topology\n#include "amber99sb-ildn.ff/tip3p.itp"'
top = top.replace(marker, '; Include triptolide molecule definition\n#include "TPL.itp"\n#ifdef POSRES_TPL\n#include "posre_TPL.itp"\n#endif\n\n' + marker)
top = top.replace("Protein\n\n[ molecules ]", "S100A8 homodimer with triptolide\n\n[ molecules ]")
top = top.rstrip() + "\nCA                  4\nTPL                 1\n"
(BUILD / "topol.top").write_text(top, encoding="utf-8")
print(BUILD / "complex_unsolvated.gro")
