from __future__ import annotations

from pathlib import Path
import shutil
import tempfile

import pymol
from pymol import cmd


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
OUT = ROOT / "publication_figure_v1" / "assets"
OUT.mkdir(parents=True, exist_ok=True)

RECEPTOR = TPL_ROOT / "input" / "prepared" / "5HLO_S100A8_dimer_BD_pH7p4.pdb"
TPL = TPL_ROOT / "results" / "local_replicates" / "BD_vina_seed107985.pdbqt"
OH16 = ROOT / "results" / "local_replicates" / "BD_vina_seed107985.pdbqt"

TMP = Path(tempfile.mkdtemp(prefix="s100a8_compare_render_"))
TMP_RECEPTOR = TMP / "receptor.pdb"
TMP_TPL = TMP / "triptolide.pdbqt"
TMP_OH16 = TMP / "hydroxytriptolide.pdbqt"
shutil.copy2(RECEPTOR, TMP_RECEPTOR)
shutil.copy2(TPL, TMP_TPL)
shutil.copy2(OH16, TMP_OH16)


def load_scene() -> None:
    cmd.reinitialize()
    cmd.load(str(TMP_RECEPTOR), "receptor")
    cmd.load(str(TMP_TPL), "tpl")
    cmd.load(str(TMP_OH16), "oh16")
    cmd.remove("solvent")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.48)
    cmd.set("direct", 0.52)
    cmd.set("specular", 0.18)
    cmd.set("shininess", 16)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("two_sided_lighting", 1)
    cmd.hide("everything", "all")
    cmd.show("cartoon", "receptor and polymer.protein")
    cmd.color("0x6FA6C1", "receptor and chain B and polymer.protein")
    cmd.color("0x8DC7B8", "receptor and chain D and polymer.protein")
    cmd.show("sticks", "tpl")
    cmd.show("sticks", "oh16")
    cmd.set("stick_radius", 0.24, "tpl or oh16")
    cmd.color("0xE88A32", "tpl and elem C")
    cmd.color("0xB24A8A", "oh16 and elem C")
    cmd.color("0xD73027", "(tpl or oh16) and elem O")
    cmd.color("white", "(tpl or oh16) and elem H")


load_scene()
cmd.select("pocket", "byres (receptor and polymer.protein within 4.5 of (tpl or oh16))")
cmd.show("sticks", "pocket")
cmd.set("stick_radius", 0.12, "pocket")
cmd.util.cnc("pocket")
cmd.orient("receptor and polymer.protein")
cmd.turn("y", 25)
cmd.turn("x", -10)
cmd.zoom("receptor and polymer.protein", 5)
tmp_overview = TMP / "docking_overlay_overview_v3.png"
cmd.png(str(tmp_overview), width=2600, height=1800, dpi=600, ray=1)
shutil.copy2(tmp_overview, OUT / "docking_overlay_overview_v3.png")

load_scene()
cmd.set("cartoon_transparency", 0.34, "receptor and polymer.protein")
cmd.select("pocket", "receptor and ((chain B and resi 43+44+45+82+85+86+89) or (chain D and resi 2+6+10+13))")
cmd.show("sticks", "pocket")
cmd.set("stick_radius", 0.16, "pocket")
cmd.util.cnc("pocket")
cmd.distance("tpl_polar", "tpl and elem O", "pocket and elem N+O", 3.5, mode=2)
cmd.distance("oh16_polar", "oh16 and elem O", "pocket and elem N+O", 3.5, mode=2)
cmd.set("dash_color", "0xE88A32", "tpl_polar")
cmd.set("dash_color", "0xB24A8A", "oh16_polar")
cmd.set("dash_width", 2.0, "tpl_polar")
cmd.set("dash_width", 2.0, "oh16_polar")
cmd.set("dash_gap", 0.28, "tpl_polar")
cmd.set("dash_gap", 0.28, "oh16_polar")
cmd.hide("labels", "tpl_polar")
cmd.hide("labels", "oh16_polar")
cmd.orient("tpl or oh16 or pocket")
cmd.turn("y", -10)
cmd.turn("x", 10)
cmd.zoom("tpl or oh16 or pocket", 3.2)
tmp_closeup = TMP / "docking_overlay_closeup_v3.png"
cmd.png(str(tmp_closeup), width=2400, height=1800, dpi=600, ray=1)
shutil.copy2(tmp_closeup, OUT / "docking_overlay_closeup_v3.png")

shutil.rmtree(TMP, ignore_errors=True)

print(OUT)
