"""Compose corrected Fig. 1 from R-generated vector panels and unchanged B/C vectors."""
from pathlib import Path
import fitz

ROOT = Path(r"D:\桌面\sepsis")
P = ROOT / "review_revision" / "10_figure_sources" / "Fig1_S1_corrected_panels"
Q = ROOT / "review_revision" / "09_figures"
OUT = ROOT / "文章" / "返修投稿新主图" / "Fig1_corrected_candidate.pdf"

SOURCES = [
    P / "Fig1A_volcano_corrected.pdf",
    Q / "Fig1B_wgcna_readable.pdf",
    Q / "Fig1C_boxplot_readable.pdf",
    P / "Fig1D_overlap_200.pdf",
    P / "Fig1E_GO_corrected.pdf",
    P / "Fig1F_KEGG_corrected.pdf",
]
RECTS = [
    fitz.Rect(18, 14, 264, 176), fitz.Rect(280, 14, 526, 176),
    fitz.Rect(18, 192, 264, 350), fitz.Rect(280, 192, 526, 350),
    fitz.Rect(18, 366, 318, 540), fitz.Rect(330, 366, 526, 540),
]

doc = fitz.open()
page = doc.new_page(width=540, height=552)
page.draw_rect(page.rect, color=None, fill=(1, 1, 1))
for i, (source, rect) in enumerate(zip(SOURCES, RECTS)):
    src = fitz.open(source)
    page.show_pdf_page(rect, src, 0, keep_proportion=True, overlay=True)
    src.close()
    page.insert_text((rect.x0 - 11, rect.y0 + 10), chr(65 + i),
                     fontname="Helvetica-Bold", fontsize=12,
                     color=(0.125, 0.145, 0.165), overlay=True)

if OUT.exists():
    OUT.unlink()
doc.save(OUT, garbage=4, deflate=True, clean=True)
doc.close()
print(OUT)
