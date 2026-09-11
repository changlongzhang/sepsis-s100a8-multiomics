from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "systems" / "complex" / "build"
BUILD = ROOT / "systems" / "apo" / "build"

if BUILD.exists():
    raise RuntimeError(f"Refusing to overwrite existing apo build: {BUILD}")
BUILD.mkdir(parents=True)


def gro_atoms(path: Path):
    lines = path.read_text(encoding="ascii").splitlines()
    count = int(lines[1].strip())
    atoms = lines[2 : 2 + count]
    if len(atoms) != count:
        raise RuntimeError(f"Truncated GRO file: {path}")
    return atoms


protein = gro_atoms(SOURCE / "protein.gro")
calcium = gro_atoms(SOURCE / "calcium.gro")
if len(protein) != 3061 or len(calcium) != 4:
    raise RuntimeError(f"Unexpected protein/calcium counts: {len(protein)}, {len(calcium)}")

combined = []
for atom_number, line in enumerate(protein + calcium, 1):
    combined.append(line[:15] + f"{atom_number:5d}" + line[20:])

(BUILD / "apo_unsolvated.gro").write_text(
    "S100A8 homodimer + 4 structural Ca2+ (apo)\n"
    + f"{len(combined):5d}\n"
    + "\n".join(combined)
    + "\n   0.00000   0.00000   0.00000\n",
    encoding="ascii",
)

for name in (
    "topol_Protein_chain_B.itp",
    "topol_Protein_chain_D.itp",
    "posre_Protein_chain_B.itp",
    "posre_Protein_chain_D.itp",
):
    shutil.copy2(SOURCE / name, BUILD / name)

topology = '''; Matched apo topology derived from the same S100A8 protein preparation
#include "amber99sb-ildn.ff/forcefield.itp"

#include "topol_Protein_chain_B.itp"
#include "topol_Protein_chain_D.itp"

#include "amber99sb-ildn.ff/tip3p.itp"

#ifdef POSRES_WATER
[ position_restraints ]
; i funct fcx fcy fcz
1 1 1000 1000 1000
#endif

#include "amber99sb-ildn.ff/ions.itp"

[ system ]
S100A8 homodimer apo with four structural Ca2+ in water

[ molecules ]
Protein_chain_B 1
Protein_chain_D 1
CA 4
'''
(BUILD / "topol.top").write_text(topology, encoding="ascii")
print(BUILD / "apo_unsolvated.gro")
