"""Repair supplementary S7 layout and redraw S8 from deposited source data."""
from pathlib import Path
import fitz, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":9,
                             "xtick.labelsize":8,"ytick.labelsize":8})
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

P=Path(r"D:\桌面\sepsis"); FINAL=P/"文章"/"返修投稿新主图"; Q=P/"review_revision"/"09_figures"
S7=P/"S100A8_triptolide_MD_topjournal"/"comparative_analysis_v3"/"publication_figures_python_v2"/"Figure_3_PCA_FEL.pdf"

def uppercase_labels_inplace(path):
    doc=fitz.open(path); page=doc[0]; hits=[]
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines",[]):
            for span in line["spans"]:
                t=span["text"].strip()
                if t in list("abcdef"):
                    r=fitz.Rect(span["bbox"]); hits.append((r,t.upper(),span["size"])); page.add_redact_annot(fitz.Rect(r.x0-1,r.y0-1,r.x1+1,r.y1+1),fill=(1,1,1))
    if hits:
        page.apply_redactions()
        for r,t,size in hits: page.insert_text((r.x0,r.y1-1),t,fontname="Helvetica-Bold",fontsize=size,color=(.08,.08,.08),overlay=True)
        tmp=path.with_name(path.stem+"_tmp_upper.pdf"); doc.save(tmp,garbage=4,deflate=True,clean=True); doc.close(); tmp.replace(path)
    else: doc.close()

def s7():
    uppercase_labels_inplace(S7)
    d=fitz.open(); page=d.new_page(width=540,height=210); src=fitz.open(S7)
    page.show_pdf_page(fitz.Rect(8,8,532,202),src,0,keep_proportion=True); src.close()
    m=FINAL/"S7_fig.pdf"; d.save(m,garbage=4,deflate=True,clean=True)
    pix=page.get_pixmap(matrix=fitz.Matrix(2250/page.rect.width,2250/page.rect.width),alpha=False); from PIL import Image
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples); im.save(FINAL/"S7_fig.tif",dpi=(300,300),compression="tiff_lzw"); im.thumbnail((1300,1300)); im.save(Q/"S7_preview.png"); d.close()
def s8():
    d=pd.read_csv(P/"实验"/"PLOS_submission_ready_final"/"Fig7_source_data.csv"); d=d[d.panel=="A"].copy()
    doses=np.array([0,10,25,50,100]); d["dose"]=d.groupby("group").cumcount().mod(5).map(dict(enumerate(doses))); d=d.dropna(subset=["value"])
    groups=["LPS 100 ng/mL","LPS 200 ng/mL","LPS 500 ng/mL"]
    means=np.full((3,5),np.nan); ns=np.zeros((3,5),int)
    for i,g in enumerate(groups):
        z=d[d.group==g]
        for j,dose in enumerate(doses):
            v=z[z.dose==dose].value.to_numpy(float); ns[i,j]=len(v); means[i,j]=np.mean(v) if len(v) else np.nan
    fig,ax=plt.subplots(figsize=(7.5,4.2)); cmap=LinearSegmentedColormap.from_list("viability",["#D56A54","#F1E7DA","#F5F6F4","#8DB4BA","#315B78"])
    im=ax.imshow(means,aspect="auto",cmap=cmap,vmin=.65,vmax=1.15)
    ax.set_xticks(range(5),[str(x) for x in doses]); ax.set_yticks(range(3),groups,fontsize=8)
    ax.set_xlabel("Triptolide concentration (nmol/L)"); ax.set_ylabel("LPS condition")
    for i in range(3):
        for j in range(5): ax.text(j,i,f"{means[i,j]:.2f}\nn={ns[i,j]}",ha="center",va="center",fontsize=8,color="white" if abs(means[i,j]-1)>.20 else "#27333D")
    for s in ax.spines.values(): s.set_visible(False)
    cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.03); cb.set_label("Relative THP-1 viability",fontsize=8); cb.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.20,right=.92,bottom=.18,top=.97)
    fig.savefig(FINAL/"S8_fig_python_master.pdf",format="pdf")
    fig.savefig(Q/"S8_preview.png",dpi=300); plt.close(fig)
if __name__=="__main__": s7(); s8()
