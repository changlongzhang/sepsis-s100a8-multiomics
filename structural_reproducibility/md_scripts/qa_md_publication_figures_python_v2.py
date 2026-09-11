from __future__ import annotations

from pathlib import Path
import csv
import xml.etree.ElementTree as ET

import fitz
from PIL import Image, ImageOps


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "comparative_analysis_v3" / "publication_figures_python_v2"
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
    paths = {extension: OUT / f"{stem}.{extension}" for extension in ("svg", "pdf", "tiff")}
    paths["png"] = PREVIEW / f"{stem}.png"
    missing = [str(path) for path in paths.values() if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"Missing or empty output: {missing}")

    with Image.open(paths["tiff"]) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (None, None))
        mode = image.mode
    if width < 1800 or height < 1200:
        raise ValueError(f"TIFF pixel dimensions too small: {stem} {width}x{height}")
    if dpi[0] is None or abs(float(dpi[0]) - 600) > 2:
        raise ValueError(f"TIFF is not 600 dpi: {stem} {dpi}")

    document = fitz.open(paths["pdf"])
    if document.page_count != 1:
        raise ValueError(f"PDF is not one page: {stem}")
    rect = document[0].rect
    pdf_width_mm = rect.width / 72 * 25.4
    pdf_height_mm = rect.height / 72 * 25.4
    document.close()
    if pdf_width_mm > 184 or pdf_height_mm > 171:
        raise ValueError(f"PDF exceeds target figure area: {stem} {pdf_width_mm:.1f}x{pdf_height_mm:.1f} mm")

    tree = ET.parse(paths["svg"])
    text_nodes = [node for node in tree.getroot().iter() if node.tag.endswith("text")]
    if len(text_nodes) < 5:
        raise ValueError(f"SVG text is not editable: {stem}")

    rows.append({
        "figure": stem,
        "tiff_width_px": width,
        "tiff_height_px": height,
        "tiff_dpi_x": dpi[0],
        "tiff_dpi_y": dpi[1],
        "tiff_mode": mode,
        "pdf_pages": 1,
        "pdf_width_mm": round(pdf_width_mm, 2),
        "pdf_height_mm": round(pdf_height_mm, 2),
        "svg_editable_text_nodes": len(text_nodes),
        "status": "PASS",
    })
    previews.append(Image.open(paths["png"]).convert("RGB"))

with (OUT / "FIGURE_TECHNICAL_QC.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

target_width = 1100
thumbs = []
for image in previews:
    ratio = target_width / image.width
    thumb = image.resize((target_width, int(image.height * ratio)), Image.Resampling.LANCZOS)
    thumbs.append(ImageOps.expand(thumb, border=8, fill="white"))
gap = 24
row_heights = [max(image.height for image in thumbs[index:index + 2])
               for index in range(0, len(thumbs), 2)]
sheet = Image.new("RGB", (target_width * 2 + gap * 3,
                          sum(row_heights) + gap * (len(row_heights) + 1)), "#ECECEC")
y = gap
for row, index in enumerate(range(0, len(thumbs), 2)):
    for column, image in enumerate(thumbs[index:index + 2]):
        sheet.paste(image, (gap + column * (target_width + gap), y))
    y += row_heights[row] + gap
sheet.save(PREVIEW / "all_figures_contact_sheet.png", dpi=(220, 220))

(OUT / "FIGURE_QC_COMPLETE.txt").write_text(
    "PASS: SVG/PDF/600 dpi TIFF/PNG bundles exist; PDFs are single-page and within target area; "
    "TIFF resolution passed; editable SVG text was detected. Visual contact sheet generated for manual inspection.\n",
    encoding="utf-8",
)
print(OUT / "FIGURE_QC_COMPLETE.txt")
