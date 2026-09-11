"""Assemble revised Fig. 7 and Fig. 8 from original images and vector plots."""
from pathlib import Path
import fitz
from PIL import Image

P=Path(r"D:\桌面\sepsis"); Q=P/"review_revision"/"09_figures"; F=P/"文章"/"返修投稿新主图"; E=P/"实验"

def place(page,label,path,rect,labelxy,clip=None):
    path=Path(path)
    if path.suffix.lower() in {".png",".tif",".tiff",".jpg",".jpeg"}:
        page.insert_image(rect,filename=str(path),keep_proportion=True,overlay=True)
    else:
        src=fitz.open(path); page.show_pdf_page(rect,src,0,keep_proportion=True,overlay=True,clip=clip); src.close()
    page.insert_text(labelxy,label,fontname="hebo",fontsize=12,color=(.09,.14,.18),overlay=True)

def fig7():
    d=fitz.open(); p=d.new_page(width=540,height=510)
    place(p,"A",Q/"Fig7_qPCR.pdf",fitz.Rect(14,8,530,150),(1.5,20))
    place(p,"B",E/"图六 WB.pdf",fitz.Rect(14,165,220,320),(1.5,177),
          fitz.Rect(699.54,339.99,1152.76,730.79))
    place(p,"C",Q/"Fig7_WB_quant.pdf",fitz.Rect(232,165,530,310),(219,177))
    place(p,"D",E/"S100A8_CETSA.pdf",fitz.Rect(14,338,290,500),(1.5,350))
    place(p,"E",E/"CETSA.png",fitz.Rect(310,360,530,480),(297,372))
    m=F/"Fig7.pdf"; d.save(m,garbage=4,deflate=True,clean=True)
    page=d[0]; z=2250/page.rect.width; pix=page.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(F/"Fig7.tif",dpi=(300,300),compression="tiff_lzw")
    im.thumbnail((1600,1600)); im.save(Q/"Fig7_title_audit.png")
    d.close()
def fig8():
    d=fitz.open(); p=d.new_page(width=540,height=490)
    # Retain additional right-side source margin so KD labels remain complete
    # after final-size font standardisation.
    clip1=fitz.Rect(78.98,62.09,530.00,402.64); clip2=fitz.Rect(95.27,43.54,490.00,382.31)
    place(p,"A",E/"补实验"/"补实验"/"PDF"/"补WB 1.pdf",fitz.Rect(14,8,190,168),(1.5,20),clip1)
    place(p,"B",Q/"Fig8_pNFkB.pdf",fitz.Rect(200,8,360,168),(187,20))
    place(p,"C",Q/"Fig8_S100A8_8group.pdf",fitz.Rect(370,8,530,168),(357,20))
    place(p,"D",E/"补实验"/"补实验"/"PDF"/"补WB 2.pdf",fitz.Rect(14,183,149,333),(1.5,195),clip2)
    place(p,"E",Q/"Fig8_CD74.pdf",fitz.Rect(159,183,273,333),(146,195))
    place(p,"F",Q/"Fig8_HLADRA.pdf",fitz.Rect(283,183,397,333),(270,195))
    place(p,"G",Q/"Fig8_S100A8_4group.pdf",fitz.Rect(407,183,521,333),(394,195))
    place(p,"H",Q/"Fig8_interaction_effect.pdf",fitz.Rect(60,348,522,480),(47,360))
    m=F/"Fig8.pdf"; d.save(m,garbage=4,deflate=True,clean=True)
    page=d[0]; z=2250/page.rect.width; pix=page.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(F/"Fig8.tif",dpi=(300,300),compression="tiff_lzw")
    im.thumbnail((1600,1600)); im.save(Q/"Fig8_final_QA.png")
    d.close()
def preview(master,out):
    d=fitz.open(master); p=d[0]; z=2250/p.rect.width; p.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False).save(out); d.close()
if __name__=="__main__": fig7(); fig8()
