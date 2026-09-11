from __future__ import annotations

from pathlib import Path
import csv
import xml.etree.ElementTree as ET

from PIL import Image, ImageOps, ImageDraw
import fitz


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "comparative_analysis_v1" / "publication_figures_python_v1"
PREVIEW = OUT / "preview"
STEMS = [
    "Figure_1_MD_system_metrics",
    "Figure_2_ligand_behavior",
    "Figure_3_PCA_FEL",
    "RMSD_latest",
    "RMSF_latest",
    "Hydrogen_bonds_latest",
]

rows = []
previews = []
for stem in STEMS:
    paths = {ext: OUT / f"{stem}.{ext}" for ext in ("svg", "pdf", "tiff")}
    paths["png"] = PREVIEW / f"{stem}.png"
    missing = [str(p) for p in paths.values() if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"Missing or empty outputs for {stem}: {missing}")

    with Image.open(paths["tiff"]) as im:
        width, height = im.size
        dpi = im.info.get("dpi", (None, None))
        tiff_mode = im.mode
    if width < 1800 or height < 1200:
        raise ValueError(f"TIFF too small for {stem}: {width} x {height}")
    if dpi[0] is None or abs(float(dpi[0]) - 600.0) > 2:
        raise ValueError(f"TIFF DPI is not 600 for {stem}: {dpi}")

    pdf = fitz.open(paths["pdf"])
    if pdf.page_count != 1:
        raise ValueError(f"PDF is not one page for {stem}")
    rect = pdf[0].rect
    pdf_width_mm = float(rect.width) / 72.0 * 25.4
    pdf_height_mm = float(rect.height) / 72.0 * 25.4
    pdf.close()

    tree = ET.parse(paths["svg"])
    root = tree.getroot()
    text_nodes = [x for x in root.iter() if x.tag.endswith("text")]
    if len(text_nodes) < 5:
        raise ValueError(f"SVG has too few editable text nodes for {stem}: {len(text_nodes)}")

    rows.append({
        "figure": stem,
        "tiff_width_px": width,
        "tiff_height_px": height,
        "tiff_dpi_x": dpi[0],
        "tiff_dpi_y": dpi[1],
        "tiff_mode": tiff_mode,
        "pdf_pages": 1,
        "pdf_width_mm": round(pdf_width_mm, 2),
        "pdf_height_mm": round(pdf_height_mm, 2),
        "svg_editable_text_nodes": len(text_nodes),
        "status": "PASS",
    })
    previews.append(Image.open(paths["png"]).convert("RGB"))

with (OUT / "FIGURE_TECHNICAL_QC.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# Python/Pillow-only contact sheet for visual inspection.
thumbs = []
target_w = 1100
for im in previews:
    ratio = target_w / im.width
    thumb = im.resize((target_w, max(1, int(im.height * ratio))), Image.Resampling.LANCZOS)
    thumbs.append(ImageOps.expand(thumb, border=8, fill="white"))
gap = 24
sheet_w = target_w * 2 + gap * 3
row_heights = []
for i in range(0, len(thumbs), 2):
    row_heights.append(max(x.height for x in thumbs[i:i + 2]))
sheet_h = sum(row_heights) + gap * (len(row_heights) + 1)
sheet = Image.new("RGB", (sheet_w, sheet_h), "#ECECEC")
y = gap
for row, i in enumerate(range(0, len(thumbs), 2)):
    for col, im in enumerate(thumbs[i:i + 2]):
        sheet.paste(im, (gap + col * (target_w + gap), y))
    y += row_heights[row] + gap
sheet.save(PREVIEW / "all_figures_contact_sheet.png", dpi=(220, 220))

(OUT / "FIGURE_QC_COMPLETE.txt").write_text(
    "PASS: expected SVG/PDF/600 dpi TIFF/PNG files exist; PDFs are single-page; TIFF resolution passed; SVG editable text nodes detected.\n"
    "Visual contact sheet generated with Python/Pillow for manual inspection.\n",
    encoding="utf-8",
)
print(OUT / "FIGURE_QC_COMPLETE.txt")
