"""Rebuild Fig. 6 from original vector/high-resolution assets and deposited MD data."""
from pathlib import Path
import shutil
import fitz
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(r"D:\桌面\sepsis")
OUT = ROOT / "文章" / "返修投稿新主图"
Q = ROOT / "review_revision" / "09_figures"
DOCK = ROOT / "S100A8_triptolide_16OH_comparative_docking_topjournal" / "publication_figure_v2_complete"
MD = ROOT / "S100A8_triptolide_MD_topjournal" / "comparative_analysis_v3" / "publication_figures_python_v2"
INK = "#27333D"; BLUE = "#4C78A8"; RED = "#D66A5E"; TEAL = "#4F938D"; ORANGE = "#D49A45"; TARGET = "#C97832"

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                     "pdf.fonttype":42,"svg.fonttype":"none","font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.75,
                     "text.color":INK,"axes.labelcolor":INK,"xtick.color":INK,"ytick.color":INK})

def redraw_a():
    drugs=["saracatinib","BI-2536","neratinib","NVP-AUY922","CGP-60474",
           "GSK-1059615","PD-0325901","BRD-K70748405","ZM 336372","Dorsomorphin",
           "16-OH triptolide","PKCβ inhibitor","trichostatin A","vorinostat",
           "L-sulforaphane","BNTX maleate","BRL 54443","EMETINE"]
    hits=[1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1,0]
    fig,ax=plt.subplots(figsize=(4.6,2.55)); ax.set_xlim(0,2); ax.set_ylim(9.25,-.45); ax.axis("off")
    for idx,(drug,hit) in enumerate(zip(drugs,hits)):
        col=idx//9; row=idx%9; x=col+.04; y=row
        is_target=drug=="16-OH triptolide"; face=TARGET if is_target else ("#F3B8B1" if hit else "#E7EBEE")
        ax.add_patch(plt.Rectangle((x,y-.34),.18,.68,facecolor=face,edgecolor="white",lw=.7))
        ax.text(x+.09,y,str(hit),ha="center",va="center",fontsize=9,fontweight="bold" if is_target else "normal",color="white" if is_target else INK)
        ax.text(x+.25,y,drug,ha="left",va="center",fontsize=9,fontweight="bold" if is_target else "normal",color=TARGET if is_target else INK)
    fig.subplots_adjust(left=.01,right=.99,bottom=.02,top=.98)
    for ext in ("pdf","svg"): fig.savefig(Q/f"Fig6A_drug_screen.{ext}")
    fig.savefig(Q/"Fig6A_drug_screen.tiff",dpi=600,pil_kwargs={"compression":"tiff_lzw"}); plt.close(fig)

def save_panel(fig, stem):
    for item in fig.findobj(match=matplotlib.text.Text):
        if item.get_text().strip() and item.get_fontsize() < 15.0:
            item.set_fontsize(15.0)
    try:
        fig.tight_layout(pad=.45)
    except Exception:
        pass
    for ext in ("pdf", "svg"):
        fig.savefig(Q / f"{stem}.{ext}")
    fig.savefig(Q / f"{stem}.tiff", dpi=600,
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)

def redraw_b():
    """Convert the RDKit/SDF-derived SVG master to a fully vector PDF panel."""
    source = DOCK / "assets" / "compound_structures_horizontal.svg"
    if not source.exists():
        raise FileNotFoundError(f"Missing SDF-derived vector structure master: {source}")
    target_svg = Q / "Fig6B_structures_labeled.svg"
    target_pdf = Q / "Fig6B_structures_labeled.pdf"
    shutil.copy2(source, target_svg)
    vector = fitz.open(source)
    target_pdf.write_bytes(vector.convert_to_pdf())
    vector.close()

def redraw_f():
    d = pd.read_csv(DOCK / "source_data" / "representative_pose_cross_compound_comparison.csv")
    d["condition"] = d.receptor + " / " + d.scoring.str.title()
    fig, ax = plt.subplots(figsize=(4.9, 3.15))
    colors = [BLUE, BLUE, TEAL, TEAL]
    bars = ax.bar(np.arange(len(d)), d.contact_jaccard, color=colors, width=.68, alpha=.90)
    for b, v in zip(bars, d.contact_jaccard):
        ax.text(b.get_x()+b.get_width()/2, v+.025, f"{v:.2f}", ha="center", va="bottom", fontsize=7.7)
    ax.set_ylim(0,1.12); ax.set_ylabel("Contact-set Jaccard")
    ax.set_xticks(range(len(d)), d.condition, rotation=28, ha="right", fontsize=7.2)
    ax.grid(axis="y",color="#E5EAED",lw=.6); ax.set_axisbelow(True)
    fig.subplots_adjust(left=.18,right=.98,bottom=.28,top=.95); save_panel(fig,"Fig6F_contact_similarity")

def redraw_e():
    d=pd.read_csv(DOCK/"source_data"/"matched_local_replicates_with_rmsd.csv")
    d["condition"]=d.receptor+" / "+d.scoring.str.title()
    conditions=["AC / Vina","BD / Vina","AC / Vinardo","BD / Vinardo"]
    wide=d.pivot_table(index=["condition","seed"],columns="compound",values="best_score_kcal_mol").reset_index()
    wide["delta"]=wide["16-hydroxytriptolide"]-wide["triptolide"]
    arrays=[wide.loc[wide.condition==c,"delta"].to_numpy() for c in conditions]
    fig,ax=plt.subplots(figsize=(5.35,3.25)); pos=np.arange(4)
    vp=ax.violinplot(arrays,positions=pos,widths=.78,showextrema=False,vert=False)
    palette=[BLUE,TEAL,ORANGE,TARGET]
    for body,c in zip(vp["bodies"],palette): body.set_facecolor(c); body.set_edgecolor(c); body.set_alpha(.25)
    bp=ax.boxplot(arrays,positions=pos,widths=.28,patch_artist=True,showfliers=False,vert=False,
                  medianprops={"color":INK,"linewidth":1.0},whiskerprops={"color":"#75818A","linewidth":.7},capprops={"color":"#75818A","linewidth":.7})
    for b,c in zip(bp["boxes"],palette): b.set_facecolor("white"); b.set_edgecolor(c); b.set_linewidth(1.0)
    rng=np.random.default_rng(12)
    for y,(v,c) in enumerate(zip(arrays,palette)):
        ax.scatter(v,y+rng.uniform(-.10,.10,len(v)),s=12,color=c,alpha=.70,edgecolor="white",linewidth=.25,zorder=3)
        ax.text(np.max(v)+.014,y,f"Δ={np.mean(v):+.2f}",va="center",fontsize=7,color=c)
    ax.axvline(0,color="#75818A",lw=.8)
    ax.set_yticks(pos,conditions,fontsize=7.2); ax.invert_yaxis()
    ax.set_xlabel("Δ docking score (16-OH − triptolide), kcal/mol")
    ax.grid(axis="x",color="#E5EAED",lw=.55); ax.set_axisbelow(True)
    ax.set_xlim(min(-.09,min(map(np.min,arrays))-.03),max(map(np.max,arrays))+.10)
    fig.subplots_adjust(left=.24,right=.98,bottom=.20,top=.96); save_panel(fig,"Fig6E_docking_effect")

def redraw_g():
    d = pd.read_csv(MD / "source_data" / "system_rmsd_replicate_median_range.csv")
    fig, ax = plt.subplots(figsize=(5.2, 3.15))
    for name, color, label in [("complex",TEAL,"Complex"),("apo",BLUE,"Apo")]:
        z=d[d.system==name].sort_values("time_ns").copy()
        z["smooth"]=z["median"].rolling(101,center=True,min_periods=1).mean()
        ax.plot(z.time_ns,z.smooth,color=color,lw=1.65,label=label)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Backbone RMSD (nm)")
    ax.legend(frameon=False,fontsize=7.5,ncol=2,loc="upper left")
    ax.grid(color="#E5EAED",lw=.55); ax.set_axisbelow(True)
    fig.subplots_adjust(left=.17,right=.98,bottom=.18,top=.95); save_panel(fig,"Fig6G_RMSD_clean")

def redraw_h():
    d=pd.read_csv(MD/"source_data"/"ligand_pocket.csv"); d["time_bin"]=(d.time_ns//5).astype(int)
    table=d.pivot_table(index="replicate",columns="time_bin",values="display_value",aggfunc="median")
    table=table.reindex(index=sorted(table.index),columns=range(int(d.time_bin.max())+1))
    fig,ax=plt.subplots(figsize=(5.2,3.15))
    cmap=LinearSegmentedColormap.from_list("pocket",["#F3F5F4","#A9C3C5",TEAL,ORANGE,RED])
    vmax=float(np.nanquantile(table.to_numpy(),.95)); im=ax.imshow(table,aspect="auto",cmap=cmap,vmin=0,vmax=vmax)
    ax.set_yticks(range(len(table.index)),[f"Rep {i}" for i in table.index],fontsize=7.5)
    ticks=np.arange(0,table.shape[1],4); ax.set_xticks(ticks,[str(int(t*5)) for t in ticks],fontsize=7.2)
    ax.set_xlabel("Time (ns; 5-ns bins)"); ax.set_ylabel("MD replicate")
    for s in ax.spines.values(): s.set_visible(False)
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.035); cb.ax.tick_params(labelsize=6.7)
    fig.subplots_adjust(left=.16,right=.91,bottom=.20,top=.95); save_panel(fig,"Fig6H_pocket_heatmap")

def place(page, label, path, rect, label_xy):
    path=Path(path)
    if path.suffix.lower() in {".png",".tif",".tiff",".jpg",".jpeg"}:
        page.insert_image(rect,filename=str(path),keep_proportion=True,overlay=True)
    else:
        src=fitz.open(path); page.show_pdf_page(rect,src,0,keep_proportion=True,overlay=True); src.close()
    page.insert_text(label_xy,label,fontname="Helvetica-Bold",fontsize=11,color=(.09,.14,.18),overlay=True)

def main():
    redraw_a(); redraw_b(); redraw_e(); redraw_f(); redraw_g(); redraw_h()
    a=Q/"Fig6A_drug_screen.pdf"
    b=Q/"Fig6B_structures_labeled.pdf"
    c=DOCK/"assets"/"docking_overlay_overview_v3.png"
    d=DOCK/"assets"/"docking_overlay_closeup_v3.png"
    e=Q/"Fig6E_docking_effect.pdf"
    f=Q/"Fig6F_contact_similarity.pdf"; g=Q/"Fig6G_RMSD_clean.pdf"; h=Q/"Fig6H_pocket_heatmap.pdf"
    doc=fitz.open(); page=doc.new_page(width=540,height=626); page.draw_rect(page.rect,color=None,fill=(1,1,1))
    cells={"A":fitz.Rect(14,8,254,168),"B":fitz.Rect(266,8,530,168),
           "C":fitz.Rect(14,180,266,320),"D":fitz.Rect(278,180,530,320),
           "E":fitz.Rect(14,332,266,458),"F":fitz.Rect(278,332,530,458),
           "G":fitz.Rect(14,472,266,615),"H":fitz.Rect(278,472,530,615)}
    paths=dict(zip("ABCDEFGH",[a,b,c,d,e,f,g,h]))
    for k in "ABCDEFGH":
        r=cells[k]; place(page,k,paths[k],r,(r.x0-12,r.y0+11))
    pdf=OUT/"Fig6.pdf"; doc.save(pdf,garbage=4,deflate=True,clean=True)
    scale=2250/page.rect.width; pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.save(OUT/"Fig6.tif",dpi=(300,300),compression="tiff_lzw")
    im.thumbnail((1600,1600)); im.save(Q/"Fig6_final_QA.png")
    doc.close()

if __name__=="__main__": main()
