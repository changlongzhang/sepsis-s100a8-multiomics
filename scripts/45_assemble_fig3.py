"""Assemble six independently redrawn vector panels into final Fig. 3."""
from pathlib import Path
import fitz
from PIL import Image


ROOT = Path(r"D:\桌面\sepsis")
Q = ROOT / "review_revision" / "09_figures"
EXT = ROOT / "sepsis-s100a8-multiomics" / "data" / "08_External Validation" / "output"
OUT = ROOT / "文章" / "返修投稿新主图"
ARIAL_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")


def place(page, label, source, rect):
    doc = fitz.open(source)
    page.show_pdf_page(rect, doc, 0, keep_proportion=True, overlay=True)
    doc.close()
    page.insert_font(fontname="ArialBold", fontfile=str(ARIAL_BOLD))
    page.insert_text((rect.x0 - 13, rect.y0 + 10), label, fontname="ArialBold",
                     fontsize=10, color=(0.10, 0.14, 0.17), overlay=True)


def main():
    panels = [
        ("A", EXT / "Combined_ROC_Validation_Final.pdf", fitz.Rect(18, 16, 196, 190)),
        ("B", Q / "Fig_ClinicalComparator_EffectForest.pdf", fitz.Rect(210, 16, 530, 190)),
        ("C", Q / "Fig_ClinicalComparator_Expression.pdf", fitz.Rect(18, 210, 272, 405)),
        ("D", Q / "Fig_ClinicalComparator_ROC.pdf", fitz.Rect(286, 210, 530, 405)),
        ("E", Q / "Fig_DecisionCurve.pdf", fitz.Rect(18, 425, 272, 602)),
        ("F", Q / "Fig_S100A8_Prognosis_Extended.pdf", fitz.Rect(286, 425, 530, 602)),
    ]
    doc = fitz.open()
    page = doc.new_page(width=540, height=620)
    page.draw_rect(page.rect, color=None, fill=(1, 1, 1))
    for label, source, rect in panels:
        place(page, label, source, rect)
    pdf = OUT / "Fig3.pdf"
    doc.save(pdf, garbage=4, deflate=True, clean=True)
    scale = 2250 / page.rect.width
    pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
    image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    image.save(OUT / "Fig3.tif", dpi=(300, 300), compression="tiff_lzw")
    image.save(OUT / "Fig3_preview.png", optimize=True)
    doc.close()


if __name__ == "__main__":
    main()
