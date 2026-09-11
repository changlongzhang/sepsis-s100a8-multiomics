"""Reassemble Fig. 1 from original vector PDFs and retain enrichment titles."""
from pathlib import Path
import fitz
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


ROOT = Path(r"D:\桌面\sepsis")
OUT = ROOT / "文章" / "返修投稿新主图"
Q = ROOT / "review_revision" / "09_figures"
SOURCES = [
    ROOT / "sepsis-s100a8-multiomics" / "data" / "02_DEG" / "output" / "Disease_vs_Normal_volcano_padj_S100A8_Fixed.pdf",
    ROOT / "sepsis-s100a8-multiomics" / "data" / "03_WGCNA_analysis" / "output" / "module_group_correlation_heatmap_all.pdf",
    ROOT / "sepsis-s100a8-multiomics" / "data" / "03_WGCNA_analysis" / "output" / "module_black_boxplot_half_border.pdf",
    ROOT / "sepsis-s100a8-multiomics" / "data" / "04_Extract overlapping genes" / "output" / "Venn_red_blue_horizontal.pdf",
    ROOT / "sepsis-s100a8-multiomics" / "data" / "05_Enrichment Analysis" / "output" / "Figure_A_GO.pdf",
    ROOT / "sepsis-s100a8-multiomics" / "data" / "05_Enrichment Analysis" / "output" / "Figure_B_KEGG.pdf",
]


def redraw_enrichment_panels():
    """Recreate the two original enrichment bar charts at their final slot size."""
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                         "pdf.fonttype":42,"font.size":8,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.linewidth":.7})
    go_terms=["NAD+ nucleosidase activity","NAD+ nucleosidase (cADPR generation)",
              "disulfide oxidoreductase activity","oxidoreductase activity, NAD(P)H to quinone",
              "oxidoreductase activity, sulfur-group donors","tertiary granule","specific granule",
              "secretory granule lumen","cytoplasmic vesicle lumen","vesicle lumen",
              "neutrophil activation","granulocyte activation","organophosphate catabolic process",
              "myeloid leukocyte activation","regulation of inflammatory response"]
    go_counts=[4,4,5,6,8,24,26,26,26,26,8,8,15,16,23]
    go_group=["MF"]*5+["CC"]*5+["BP"]*5
    colors={"BP":"#D66A5E","CC":"#4C78A8","MF":"#59A14F"}
    fig,ax=plt.subplots(figsize=(4.17,2.42)); y=np.arange(len(go_terms))
    ax.barh(y,go_counts,color=[colors[g] for g in go_group],height=.68)
    ax.set_yticks(y,go_terms,fontsize=8); ax.tick_params(axis="y",length=0,pad=2); ax.tick_params(axis="x",labelsize=8)
    ax.set_xlabel("Gene count",fontsize=8); ax.grid(axis="x",color="#E7EBEE",lw=.5); ax.set_axisbelow(True)
    for yy,v in zip(y,go_counts): ax.text(v+.35,yy,str(v),va="center",fontsize=8)
    fig.subplots_adjust(left=.61,right=.97,bottom=.18,top=.93); fig.savefig(Q/"Fig1E_GO_readable.pdf"); plt.close(fig)

    kegg_terms=["Fatty acid\nbiosynthesis","Ovarian\nsteroidogenesis","Fluid shear stress\nand atherosclerosis",
                "Fc epsilon RI\nsignaling pathway","Complement and\ncoagulation cascades"]
    kegg_counts=[3,5,7,5,6]
    kegg_group=["Metabolism","Related","Related","Immune","Immune"]
    kcols={"Related":"#4C78A8","Metabolism":"#59A14F","Immune":"#D66A5E"}
    fig,ax=plt.subplots(figsize=(2.72,2.42)); y=np.arange(len(kegg_terms))
    ax.barh(y,kegg_counts,color=[kcols[g] for g in kegg_group],height=.62)
    ax.set_yticks(y,kegg_terms,fontsize=8); ax.tick_params(axis="y",length=0,pad=2); ax.tick_params(axis="x",labelsize=8)
    ax.set_xlabel("Gene count",fontsize=8); ax.grid(axis="x",color="#E7EBEE",lw=.5); ax.set_axisbelow(True)
    for yy,v in zip(y,kegg_counts): ax.text(v+.12,yy,str(v),va="center",fontsize=8)
    fig.subplots_adjust(left=.64,right=.96,bottom=.18,top=.93); fig.savefig(Q/"Fig1F_KEGG_readable.pdf"); plt.close(fig)


def redraw_volcano_and_wgcna():
    deg=pd.read_csv(ROOT/"sepsis-s100a8-multiomics"/"data"/"02_DEG"/"output"/"Disease_vs_Normal_all_results_padj.csv")
    fig,ax=plt.subplots(figsize=(3.42,2.25)); reg=deg.Regulation.fillna("Not significant")
    palette={"Down":"#4C78A8","Not significant":"#C8CDD0","Up":"#D66A5E"}
    for key in ["Not significant","Down","Up"]:
        z=deg[reg.eq(key)]; ax.scatter(z.log2FoldChange,-np.log10(z.padj.clip(lower=1e-300)),s=5,color=palette[key],alpha=.65,edgecolors="none",label="NS" if key=="Not significant" else key)
    ax.axvline(-1,color="#7B858B",lw=.7,ls=(0,(3,3))); ax.axvline(1,color="#7B858B",lw=.7,ls=(0,(3,3))); ax.axhline(-np.log10(.05),color="#7B858B",lw=.7,ls=(0,(3,3)))
    hit=deg[deg.Gene.eq("S100A8")]
    if len(hit):
        r=hit.iloc[0]; ax.annotate("S100A8",(r.log2FoldChange,-np.log10(max(r.padj,1e-300))),xytext=(0,8),textcoords="offset points",ha="center",fontsize=8,fontweight="bold")
    ax.set_xlabel("Log2FC",fontsize=8); ax.set_ylabel("−log10(FDR)",fontsize=8); ax.tick_params(labelsize=8)
    ax.legend(frameon=False,fontsize=8,loc="center right",handletextpad=.3,labelspacing=.2)
    fig.subplots_adjust(left=.18,right=.97,bottom=.20,top=.95); fig.savefig(Q/"Fig1A_volcano_readable.pdf"); plt.close(fig)

    assoc=pd.read_csv(ROOT/"sepsis-s100a8-multiomics"/"data"/"03_WGCNA_analysis"/"output"/"module_group_association.csv")
    names=assoc.Module.str.replace("ME","",regex=False).tolist(); corr=assoc.Correlation.to_numpy(float)
    mat=np.column_stack([-corr,corr]); cmap=LinearSegmentedColormap.from_list("assoc",["#4C78A8","#F3F3F2","#D66A5E"])
    fig,ax=plt.subplots(figsize=(3.42,2.25)); im=ax.imshow(mat,aspect="auto",cmap=cmap,norm=TwoSlopeNorm(vmin=-.65,vcenter=0,vmax=.65))
    ax.set_yticks(range(len(names)),names,fontsize=8); ax.set_xticks([0,1],["Normal","Disease"],fontsize=8,rotation=35,ha="right")
    for i,r in assoc.iterrows():
        star=r.Significance if r.Significance!="NS" else "NS"
        for j,v in enumerate([-r.Correlation,r.Correlation]): ax.text(j,i,f"{v:+.2f} {star}",ha="center",va="center",fontsize=8,color="white" if abs(v)>.42 else "#26323A")
    for s in ax.spines.values():s.set_visible(False)
    cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03);cb.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.28,right=.91,bottom=.24,top=.96);fig.savefig(Q/"Fig1B_wgcna_readable.pdf");plt.close(fig)


def redact_redundant_group_label(source, output):
    doc=fitz.open(source); page=doc[0]
    for hit in page.search_for("Group"):
        page.add_redact_annot(fitz.Rect(hit.x0-2,hit.y0-2,hit.x1+2,hit.y1+2),fill=(1,1,1))
    page.apply_redactions(); doc.save(output,garbage=4,deflate=True,clean=True);doc.close()


def place(page, label, source, rect):
    doc = fitz.open(source)
    page.show_pdf_page(rect, doc, 0, keep_proportion=True, overlay=True)
    doc.close()
    page.insert_text((rect.x0 - 11, rect.y0 + 10), label, fontname="Helvetica-Bold",
                     fontsize=10, color=(0.10, 0.14, 0.17), overlay=True)


def main():
    redraw_enrichment_panels()
    redraw_volcano_and_wgcna()
    clean_c=Q/"Fig1C_boxplot_readable.pdf"; redact_redundant_group_label(SOURCES[2],clean_c)
    sources=[Q/"Fig1A_volcano_readable.pdf",Q/"Fig1B_wgcna_readable.pdf",clean_c,SOURCES[3],Q/"Fig1E_GO_readable.pdf",Q/"Fig1F_KEGG_readable.pdf"]
    rects = [
        fitz.Rect(18, 14, 264, 176), fitz.Rect(280, 14, 526, 176),
        fitz.Rect(18, 192, 264, 350), fitz.Rect(280, 192, 526, 350),
        fitz.Rect(18, 366, 318, 540), fitz.Rect(330, 366, 526, 540),
    ]
    doc = fitz.open(); page = doc.new_page(width=540, height=552)
    page.draw_rect(page.rect, color=None, fill=(1, 1, 1))
    for i, (source, rect) in enumerate(zip(sources, rects)):
        place(page, chr(65 + i), source, rect)
    # These two panel-level analytical headings are intentionally retained.
    page.insert_textbox(fitz.Rect(54, 367, 282, 382), "GO Enrichment Analysis",
                        fontname="Helvetica-Bold", fontsize=8.5, align=1,
                        color=(0.10, 0.14, 0.17), overlay=True)
    page.insert_textbox(fitz.Rect(348, 367, 512, 382), "KEGG Enrichment Analysis",
                        fontname="Helvetica-Bold", fontsize=8.5, align=1,
                        color=(0.10, 0.14, 0.17), overlay=True)
    pdf = OUT / "Fig1.pdf"
    doc.save(pdf, garbage=4, deflate=True, clean=True)
    scale = 2250 / page.rect.width
    pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
    image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    image.save(OUT / "Fig1.tif", dpi=(300, 300), compression="tiff_lzw")
    image.thumbnail((1600, 1600)); image.save(Q / "Fig1_title_audit.png")
    doc.close()


if __name__ == "__main__":
    main()
