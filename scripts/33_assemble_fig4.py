"""Assemble the revised four-panel Fig. 4 as a vector PDF master."""
from pathlib import Path

import fitz
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, ListedColormap


PROJECT = Path(r"D:\桌面\sepsis")
SC = (PROJECT / "sepsis-s100a8-multiomics" / "data" /
      "09_scRNA-seq & Bulk Immune Infiltration" / "output")
REV = PROJECT / "review_revision" / "09_figures"
FINAL = PROJECT / "文章" / "返修投稿新主图"

SOURCES = {
    "A": REV / "Fig4A_Cleaned_UMAP_readable.pdf",
    "B": SC / "Plot_05_S100A8_FeaturePlot.pdf",
    "C": REV / "Fig4C_Annotation_DotPlot_readable.pdf",
    "D": REV / "Fig_Monocyte_UnsupervisedClusters.pdf",
    "E": REV / "Fig_Monocyte_StatePrograms.pdf",
}
CELLS = {
    "A": fitz.Rect(14, 8, 266, 166),
    "B": fitz.Rect(278, 8, 530, 166),
    "C": fitz.Rect(14, 176, 530, 363),
    "D": fitz.Rect(14, 378, 266, 598),
    "E": fitz.Rect(278, 378, 530, 598),
}
LABELS = {"A": (1.5, 20), "B": (265, 20), "C": (1.5, 188),
          "D": (1.5, 390), "E": (265, 390)}


def redraw_panel_d() -> None:
    """Show the original donor-level activities with a compact effect summary."""
    stats = pd.read_csv(REV / "Fig_DonorLevel_TFActivity_source_data.csv")
    raw = pd.read_csv(PROJECT / "review_revision" / "05_single_cell_pseudobulk" / "donor_level_TF_activity.csv")
    donors = raw.sort_values(["group", "donor_id"], ascending=[True, True]).donor_id.drop_duplicates().tolist()
    tfs = stats.sort_values("effect", ascending=False).tf.tolist()
    matrix = raw.pivot(index="tf", columns="donor_id", values="activity").reindex(index=tfs, columns=donors)
    # Within-TF standardisation preserves donor heterogeneity while making TF rows comparable.
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
    groups = raw.drop_duplicates("donor_id").set_index("donor_id").reindex(donors).group
    ink = "#27333D"
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9,
                         "text.color": ink, "axes.labelcolor": ink, "xtick.color": ink, "ytick.color": ink})
    fig = plt.figure(figsize=(7.2, 3.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[.13, 1], width_ratios=[3.1, 1.15], hspace=.05, wspace=.28)
    strip = fig.add_subplot(gs[0,0]); ax = fig.add_subplot(gs[1,0]); eff = fig.add_subplot(gs[1,1])
    cmap = LinearSegmentedColormap.from_list("tf", ["#4C78A8", "#F2F3F2", "#D66A5E"])
    im = ax.imshow(z, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2))
    ax.set_xticks(range(len(donors)), donors, rotation=35, ha="right", fontsize=7.2)
    ax.set_yticks(range(len(tfs)), tfs, fontsize=7.6)
    for y in range(len(tfs)):
        for x in range(len(donors)):
            v=z.iloc[y,x]
            ax.text(x,y,f"{v:+.1f}",ha="center",va="center",fontsize=6.5,
                    color="white" if abs(v)>1.05 else ink)
    ax.set_xlabel("Donor-level TF activity (within-TF z-score)", fontsize=7.8)
    ax.set_xticks(np.arange(-.5, len(donors), 1), minor=True); ax.set_yticks(np.arange(-.5, len(tfs), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2); ax.tick_params(which="minor", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.025); cb.set_label("Within-TF z-score", fontsize=7); cb.ax.tick_params(labelsize=6.5)

    group_codes=np.array([[0 if g=="Control" else 1 for g in groups]])
    strip.imshow(group_codes,aspect="auto",cmap=ListedColormap(["#4C78A8","#D66A5E"]),vmin=0,vmax=1)
    strip.set_xticks(range(len(donors)),[]); strip.set_yticks([])
    strip.text(-.12,.5,"Group",transform=strip.transAxes,ha="right",va="center",fontsize=7.2,fontweight="bold")
    for s in strip.spines.values(): s.set_visible(False)

    ordered=stats.set_index("tf").reindex(tfs)
    y=np.arange(len(tfs)); colors=np.where(ordered.effect>=0,"#D66A5E","#4C78A8")
    eff.barh(y,ordered.effect,color=colors,alpha=.82,height=.62)
    eff.axvline(0,color="#7F898F",lw=.75)
    eff.set_yticks(y,[]); eff.invert_yaxis(); eff.set_xlabel("Sepsis − control\neffect",fontsize=7.5)
    eff.grid(axis="x",color="#E7EBEE",lw=.55); eff.set_axisbelow(True)
    for yy,r in enumerate(ordered.itertuples()):
        eff.text(r.effect+(.04 if r.effect>=0 else -.04),yy,f"{r.effect:+.2f}",va="center",
                 ha="left" if r.effect>=0 else "right",fontsize=6.8)
    eff.text(.5,1.03,"FDR range 0.286–0.905",transform=eff.transAxes,ha="center",fontsize=6.7,color="#65727E")
    fig.subplots_adjust(left=.13, right=.98, bottom=.24, top=.94)
    base = REV / "Fig_DonorLevel_TFActivity"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.03)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", pad_inches=.03)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=.03,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=.03)
    plt.close(fig)


def remove_source_title(path: Path, title: str) -> None:
    doc=fitz.open(path); page=doc[0]; hits=page.search_for(title)
    for r in hits:
        rr=fitz.Rect(r.x0-2,r.y0-2,r.x1+2,r.y1+2); page.add_redact_annot(rr,fill=(1,1,1))
    if hits:
        page.apply_redactions(); tmp=path.with_name(path.stem+"_tmp_notitle.pdf")
        doc.save(tmp,garbage=4,deflate=True,clean=True); doc.close(); tmp.replace(path)
    else: doc.close()


def prepare_panel_c() -> None:
    """Retain the original dot plot while replacing its colliding group headers."""
    source = SC / "Plot_02_Annotation_DotPlot.pdf"
    output = SOURCES["C"]
    doc = fitz.open(source)
    page = doc[0]
    # Only the six top group headers are redacted; the matrix, genes and legends remain intact.
    header_boxes = [
        fitz.Rect(58, 0, 180, 12), fitz.Rect(193, 0, 290, 12),
        fitz.Rect(300, 0, 381, 12), fitz.Rect(438, 0, 471, 12),
        fitz.Rect(554, 0, 584, 12), fitz.Rect(601, 0, 635, 12),
    ]
    for box in header_boxes:
        page.add_redact_annot(box, fill=(1, 1, 1))
    page.apply_redactions()
    headers = [
        ("Classical myeloid", 119, (0.30, 0.47, 0.66)),
        ("Non-classical", 241, (0.31, 0.58, 0.55)),
        ("Antigen-presenting", 340, (0.83, 0.60, 0.27)),
        ("T / NK", 454, (0.77, 0.38, 0.34)),
        ("B cell", 569, (0.49, 0.42, 0.65)),
        ("Platelet", 618, (0.45, 0.51, 0.55)),
    ]
    for label, center, color in headers:
        width = fitz.get_text_length(label, fontname="Helvetica-Bold", fontsize=6.2)
        page.draw_line((center - width / 2, 1.0), (center + width / 2, 1.0), color=color, width=1.2)
        page.insert_text((center - width / 2, 8.2), label, fontname="Helvetica-Bold",
                         fontsize=6.2, color=(0.12, 0.17, 0.20), overlay=True)
    doc.save(output, garbage=4, deflate=True, clean=True)
    doc.close()


def prepare_panel_a() -> None:
    """Shift the original UMAP legend left so final-size type remains in-panel."""
    source = SC / "Plot_03_Cleaned_UMAP.pdf"
    output = SOURCES["A"]
    doc = fitz.open(source)
    page = doc[0]
    page.add_redact_annot(fitz.Rect(480, 163, 576, 248), fill=(1, 1, 1))
    page.apply_redactions()
    entries = [
        ("T Cells", 187, (0.17, 0.63, 0.17)),
        ("B Cells", 202, (0.12, 0.47, 0.71)),
        ("NK Cells", 217, (0.63, 0.82, 0.39)),
        ("Monocytes", 232, (0.89, 0.10, 0.11)),
    ]
    for label, y, color in entries:
        page.draw_circle((454, y - 4), 3.3, color=None, fill=color)
        page.insert_text((465, y), label, fontname="Helvetica", fontsize=12,
                         color=(0.10, 0.14, 0.17), overlay=True)
    doc.save(output, garbage=4, deflate=True, clean=True)
    doc.close()


def main() -> None:
    prepare_panel_a()
    prepare_panel_c()
    remove_source_title(SOURCES["B"], "S100A8 Expression")
    document = fitz.open()
    page = document.new_page(width=540, height=606)
    for label in ("A", "B", "C", "D", "E"):
        source = fitz.open(SOURCES[label])
        page.show_pdf_page(CELLS[label], source, 0, keep_proportion=True, overlay=True)
        source.close()
        x, y = LABELS[label]
        page.insert_text((x, y), label, fontname="hebo", fontsize=12,
                         color=(0.09, 0.14, 0.18), overlay=True)
    master = FINAL / "Fig4_python_master.pdf"
    document.save(master, garbage=4, deflate=True, clean=True)
    document.close()

    check = fitz.open(master)
    page = check[0]
    zoom = 2250 / page.rect.width
    page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False).save(
        FINAL / "Fig4_preview.png"
    )
    check.close()

    # Submission deliverables: editable PDF/AI plus RGB 300-dpi LZW TIFF.
    for suffix in ("pdf", "ai"):
        target = FINAL / f"Fig4.{suffix}"
        target.write_bytes(master.read_bytes())
    check = fitz.open(master); page = check[0]
    zoom = 2250 / page.rect.width
    pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    image.save(FINAL / "Fig4.tif", dpi=(300, 300), compression="tiff_lzw")
    check.close()


if __name__ == "__main__":
    main()
