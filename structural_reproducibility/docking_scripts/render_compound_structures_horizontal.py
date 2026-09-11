from pathlib import Path
import shutil
import tempfile

from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
OUT = ROOT / "publication_figure_v2_complete" / "assets"
OUT.mkdir(parents=True, exist_ok=True)

tmp = Path(tempfile.mkdtemp(prefix="s100a8_compounds_horizontal_"))
tpl_path = tmp / "triptolide.sdf"
oh_path = tmp / "hydroxytriptolide.sdf"
shutil.copy2(TPL_ROOT / "input" / "prepared" / "triptolide_primary.sdf", tpl_path)
shutil.copy2(ROOT / "input" / "prepared" / "16-hydroxytriptolide_primary.sdf", oh_path)

tpl = Chem.SDMolSupplier(str(tpl_path), removeHs=True)[0]
oh = Chem.SDMolSupplier(str(oh_path), removeHs=True)[0]
if tpl is None or oh is None:
    raise RuntimeError("Could not read archived compound structures")

options = Draw.MolDrawOptions()
options.padding = 0.035
options.bondLineWidth = 3.0
options.legendFontSize = 60
options.minFontSize = 22
options.maxFontSize = 38
options.addStereoAnnotation = False

image = Draw.MolsToGridImage(
    [tpl, oh], molsPerRow=2, subImgSize=(900, 650),
    legends=["Triptolide", "16-hydroxytriptolide"],
    useSVG=False, drawOptions=options,
)
image.save(OUT / "compound_structures_horizontal.png", dpi=(600, 600))

# Vector master drawn directly from the archived SDF structures. This is the
# preferred manuscript source because bonds, atom labels and compound names
# remain sharp at any magnification.
svg = Draw.MolsToGridImage(
    [tpl, oh], molsPerRow=2, subImgSize=(900, 650),
    legends=["Triptolide", "16-hydroxytriptolide"],
    useSVG=True, drawOptions=options,
)
(OUT / "compound_structures_horizontal.svg").write_text(svg, encoding="utf-8")

vertical = Draw.MolsToGridImage(
    [tpl, oh], molsPerRow=1, subImgSize=(1000, 550),
    legends=["Triptolide", "16-hydroxytriptolide"],
    useSVG=False, drawOptions=options,
)
vertical.save(OUT / "compound_structures_vertical_v2.png", dpi=(600, 600))
shutil.rmtree(tmp, ignore_errors=True)
print(OUT / "compound_structures_horizontal.png")
