"""Assemble revised Fig. 5 from independent vector PDF panels."""
from pathlib import Path
import fitz
from PIL import Image

PROJECT=Path(r"D:\桌面\sepsis")
SRC=PROJECT/"review_revision"/"09_figures"
FINAL=PROJECT/"文章"/"返修投稿新主图"
STEMS={
 "A":"Fig_SCP548_DonorPseudobulk_S100A8","B":"Fig_SCP548_MS1_Abundance",
 "C":"Fig_SCP548_AntigenPresentation","D":"Fig_GSE205672_S100A8_Monocytes",
 "E":"Fig_GSE205672_AntigenPresentation","F":"Fig_Mechanistic_Triangulation"}
CELLS={
 "A":fitz.Rect(14,8,266,193), "B":fitz.Rect(278,8,530,193),
 "C":fitz.Rect(14,208,266,393), "D":fitz.Rect(278,208,530,393),
 "E":fitz.Rect(14,408,266,598), "F":fitz.Rect(278,408,530,598)}
LABELS={"A":(1.5,20),"B":(265,20),"C":(1.5,220),"D":(265,220),
        "E":(1.5,420),"F":(265,420)}

def main():
    doc=fitz.open(); page=doc.new_page(width=540,height=610)
    for k in "ABCDEF":
        src=fitz.open(SRC/f"{STEMS[k]}.pdf")
        page.show_pdf_page(CELLS[k],src,0,keep_proportion=True,overlay=True); src.close()
        x,y=LABELS[k]; page.insert_text((x,y),k,fontname="hebo",fontsize=12,color=(.09,.14,.18),overlay=True)
    master=FINAL/"Fig5_python_master.pdf"; doc.save(master,garbage=4,deflate=True,clean=True); doc.close()
    q=fitz.open(master); p=q[0]; z=2250/p.rect.width
    p.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False).save(FINAL/"Fig5_preview.png"); q.close()
    for suffix in ("pdf", "ai"):
        (FINAL/f"Fig5.{suffix}").write_bytes(master.read_bytes())
    q=fitz.open(master); p=q[0]; z=2250/p.rect.width
    pix=p.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(FINAL/"Fig5.tif",dpi=(300,300),compression="tiff_lzw"); q.close()
if __name__=="__main__": main()
