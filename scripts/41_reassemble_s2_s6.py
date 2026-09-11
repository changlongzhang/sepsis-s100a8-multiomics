"""Redraw S2/S3B as independent panels and assemble S2-S6 without nested mini-figures."""
from pathlib import Path
import fitz
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Patch

ROOT=Path(r"D:\桌面\sepsis"); FIG=ROOT/"review_revision"/"09_figures"; OUT=ROOT/"文章"/"返修投稿新主图"
MD=ROOT/"S100A8_triptolide_MD_topjournal"/"comparative_analysis_v3"/"publication_figures_python_v2"
INK="#27333D"; BLUE="#3E718C"; RED="#C65F4A"; ORANGE="#C98B42"; TEAL="#3F8587"; MUTED="#65727E"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"pdf.fonttype":42,"svg.fonttype":"none","font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.75,"text.color":INK,"axes.labelcolor":INK,"xtick.color":INK,"ytick.color":INK})

S2_STEMS=["Fig_ClinicalComparator_PR","Fig_S100A8_APACHEII","Fig_S100A8_Mortality","Fig_S100A8_Longitudinal_Trajectory","Fig_S100A8_TreatmentResponse","Fig_Prognosis_Forest"]
S3_STEMS=["Fig_Donor_UMAP","Fig_Donor_Composition","Fig_Donor_CellNumber","Fig_Donor_StateAbundance","Fig_Monocyte_Clustree","Fig_Monocyte_UMAP_ByDonor"]
S4_STEMS=["Fig_Pseudobulk_S100A8_Donor","Fig_IndependentSignature_DonorScores","Fig_Bulk_IndependentSignature","Fig_SCP548_S100A8_HLAII_Correlation","Fig_GSE205672_S100A8_HLAII_Correlation","Fig_GSE205672_NFkBScore"]

def read(stem): return pd.read_csv(FIG/f"{stem}_source_data.csv")
def save(fig,stem):
    for ext in ("pdf","svg"): fig.savefig(FIG/f"{stem}.{ext}",bbox_inches="tight",pad_inches=.03)
    fig.savefig(FIG/f"{stem}.tiff",dpi=600,bbox_inches="tight",pad_inches=.03,pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(FIG/f"{stem}.png",dpi=300,bbox_inches="tight",pad_inches=.03); plt.close(fig)

def s2a():
    d=read("Fig_ClinicalComparator_PR"); rows=[]
    for key,z in d.groupby("comparison_id",sort=False):
        z=z.sort_values("recall"); au=float(np.trapezoid(z.precision,z.recall)); base=float(z.prevalence.iloc[0]); rows.append((key,au,base,au/base))
    z=pd.DataFrame(rows,columns=["comparison","AUPRC","prevalence","lift"]).sort_values("lift",ascending=False)
    mat=z[["AUPRC","prevalence","lift"]].copy(); mat["lift"]=np.clip((mat["lift"]-1)/1.0,0,1)
    fig,ax=plt.subplots(figsize=(5.4,4.0)); cmap=LinearSegmentedColormap.from_list("pr",["#F2F4F4","#A8C2CC",BLUE])
    im=ax.imshow(mat,aspect="auto",cmap=cmap,vmin=0,vmax=1)
    labels=[x.replace("_"," | ") for x in z.comparison]
    ax.set_yticks(range(len(z)),labels,fontsize=6.8); ax.set_xticks(range(3),["AUPRC","Baseline","Lift vs baseline"],fontsize=7.2)
    for i,r in z.reset_index(drop=True).iterrows():
        for j,t in enumerate([f"{r.AUPRC:.2f}",f"{r.prevalence:.2f}",f"{r.lift:.2f}×"]): ax.text(j,i,t,ha="center",va="center",fontsize=6.8,color="white" if mat.iloc[i,j]>.58 else INK)
    for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.40,right=.98,bottom=.14,top=.98); save(fig,"Fig_ClinicalComparator_PR")

def s2b():
    d=read("Fig_S100A8_APACHEII"); d["severity"]=pd.cut(d.apache_ii,[-np.inf,15,24,np.inf],labels=["≤15","16–24","≥25"]); d["outcome"]=np.where(d.death==1,"Non-survivor","Survivor")
    rows=["Survivor","Non-survivor"]; cols=["≤15","16–24","≥25"]
    med=np.full((2,3),np.nan); n=np.zeros((2,3),int)
    for i,r in enumerate(rows):
        for j,c in enumerate(cols):
            v=d[(d.outcome==r)&(d.severity.astype(str)==c)].expression; n[i,j]=len(v); med[i,j]=v.median() if len(v) else np.nan
    fig,ax=plt.subplots(figsize=(5.0,3.8)); cmap=LinearSegmentedColormap.from_list("sev",["#EAF0F2",BLUE,RED]); im=ax.imshow(med,aspect="auto",cmap=cmap)
    ax.set_xticks(range(3),cols); ax.set_yticks(range(2),rows); ax.set_xlabel("APACHE II category")
    for i in range(2):
        for j in range(3):
            ax.text(j,i,"—" if np.isnan(med[i,j]) else f"{med[i,j]:.2f}\nn={n[i,j]}",ha="center",va="center",fontsize=8,color="white" if not np.isnan(med[i,j]) and med[i,j]>np.nanmedian(med) else INK)
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.05); cb.set_label("Median S100A8 expression",fontsize=7.3); cb.ax.tick_params(labelsize=6.8)
    for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.26,right=.88,bottom=.20,top=.97); save(fig,"Fig_S100A8_APACHEII")

def s2c():
    d=read("Fig_S100A8_Mortality"); datasets=list(d.dataset.drop_duplicates()); groups=["Survivor","Non-survivor"]
    mat=np.zeros((len(datasets),2)); ns=np.zeros_like(mat,int)
    for i,ds in enumerate(datasets):
        for j,g in enumerate(groups):
            v=d[(d.dataset==ds)&(d.mortality_label==g)].S100A8_z; mat[i,j]=v.median(); ns[i,j]=len(v)
    fig,ax=plt.subplots(figsize=(5.0,3.8)); im=ax.imshow(mat,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1))
    ax.set_xticks(range(2),groups); ax.set_yticks(range(len(datasets)),datasets)
    for i in range(len(datasets)):
        for j in range(2): ax.text(j,i,f"{mat[i,j]:+.2f}\nn={ns[i,j]}",ha="center",va="center",fontsize=7.5,color="white" if abs(mat[i,j])>.55 else INK)
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.05); cb.set_label("Median within-cohort z-score",fontsize=7.2); cb.ax.tick_params(labelsize=6.7)
    for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.25,right=.88,bottom=.16,top=.97); save(fig,"Fig_S100A8_Mortality")

def s2d():
    d=read("Fig_S100A8_Longitudinal_Trajectory"); p=d.pivot_table(index=["outcome","patient_id"],columns="day",values="expression",aggfunc="mean")
    p=p.reindex(sorted(p.index,key=lambda x:(x[0]!="Non-survivor",str(x[1])))); arr=p.to_numpy()
    fig,ax=plt.subplots(figsize=(5.0,4.0)); cmap=plt.get_cmap("viridis").copy(); cmap.set_bad("#E6E9EA"); im=ax.imshow(np.ma.masked_invalid(arr),aspect="auto",cmap=cmap)
    ax.set_xticks(range(len(p.columns)),[str(x) for x in p.columns]); ax.set_xlabel("Day after enrollment")
    ticks=np.linspace(0,len(p)-1,min(8,len(p)),dtype=int); ax.set_yticks(ticks,[f"{p.index[i][0][:2]}-{p.index[i][1]}" for i in ticks],fontsize=6.5); ax.set_ylabel("Patient (outcome prefix)")
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.035); cb.set_label("S100A8 expression",fontsize=7.2); cb.ax.tick_params(labelsize=6.6)
    for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.23,right=.91,bottom=.17,top=.98); save(fig,"Fig_S100A8_Longitudinal_Trajectory")

def s2e():
    d=read("Fig_S100A8_TreatmentResponse"); groups=["NR","R"]; fig,ax=plt.subplots(figsize=(4.8,3.8)); vals=[d[d.response==g].delta_T2_minus_T1.to_numpy() for g in groups]
    vp=ax.violinplot(vals,[0,1],widths=.75,showextrema=False)
    for body,c in zip(vp["bodies"],[RED,BLUE]): body.set_facecolor(c); body.set_edgecolor(c); body.set_alpha(.20)
    bp=ax.boxplot(vals,positions=[0,1],widths=.28,patch_artist=True,showfliers=False,medianprops={"color":INK})
    for b,c in zip(bp["boxes"],[RED,BLUE]): b.set_facecolor("white"); b.set_edgecolor(c)
    rng=np.random.default_rng(9)
    for i,(v,c) in enumerate(zip(vals,[RED,BLUE])): ax.scatter(i+rng.uniform(-.12,.12,len(v)),v,s=18,color=c,alpha=.72,edgecolor="white",linewidth=.3)
    ax.axhline(0,color="#89949C",lw=.8); ax.set_xticks([0,1],["Non-responder","Responder"]); ax.set_ylabel("ΔS100A8 logCPM (T2 − T1)")
    ax.grid(axis="y",color="#E5EAED",lw=.6); ax.set_axisbelow(True); fig.subplots_adjust(left=.20,right=.98,bottom=.16,top=.97); save(fig,"Fig_S100A8_TreatmentResponse")

def s2f():
    d=read("Fig_Prognosis_Forest").copy(); d["short"]=d.dataset+" | "+d.model.str.replace("_"," ",regex=False); d=d.iloc[::-1]
    mat=d[["log_effect"]].to_numpy(); lim=max(.6,float(np.max(np.abs(mat))))
    fig,(ax,txt)=plt.subplots(1,2,figsize=(6.0,3.8),gridspec_kw={"width_ratios":[1,2.3]}); im=ax.imshow(mat,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    ax.set_xticks([0],["log(OR/HR)"]); ax.set_yticks(range(len(d)),d.short,fontsize=6.5)
    for i,r in enumerate(d.itertuples()): ax.text(0,i,f"{r.log_effect:+.2f}",ha="center",va="center",fontsize=7,color="white" if abs(r.log_effect)>.35 else INK)
    for s in ax.spines.values(): s.set_visible(False)
    txt.set_xlim(0,1); txt.set_ylim(len(d)-.5,-.5); txt.axis("off"); txt.text(.02,-.65,"Estimate [95% CI]",fontsize=7.5,fontweight="bold"); txt.text(.98,-.65,"p",fontsize=7.5,fontweight="bold",ha="right")
    for i,r in enumerate(d.itertuples()):
        txt.add_patch(plt.Rectangle((0,i-.48),1,.96,facecolor="#F5F7F7" if i%2==0 else "white",edgecolor="none")); txt.text(.02,i,f"{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]",va="center",fontsize=6.7); txt.text(.98,i,f"{r.p_value:.3f}",ha="right",va="center",fontsize=6.7)
    fig.subplots_adjust(left=.34,right=.98,bottom=.16,top=.92,wspace=.38); save(fig,"Fig_Prognosis_Forest")

def s3b():
    d=read("Fig_Donor_Composition"); donors=list(d.donor_id.drop_duplicates()); cells=["T Cells","B Cells","NK Cells","Monocytes"]
    mat=d.pivot(index="donor_id",columns="cell_type",values="proportion").reindex(index=donors,columns=cells)
    fig,ax=plt.subplots(figsize=(5.2,3.8)); cmap=LinearSegmentedColormap.from_list("comp",["#F3F5F5","#9FBEC9",BLUE]); im=ax.imshow(mat,aspect="auto",cmap=cmap,vmin=0,vmax=.65)
    ax.set_xticks(range(len(cells)),cells,rotation=25,ha="right",fontsize=7.4); ax.set_yticks(range(len(donors)),donors,fontsize=7.4)
    for i in range(len(donors)):
        for j in range(len(cells)):
            v=mat.iloc[i,j]; ax.text(j,i,f"{v*100:.1f}%",ha="center",va="center",fontsize=6.8,color="white" if v>.38 else INK)
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.035); cb.set_label("Cell-type proportion",fontsize=7.2); cb.ax.tick_params(labelsize=6.6)
    for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.17,right=.91,bottom=.23,top=.98); save(fig,"Fig_Donor_Composition")

def s3c():
    d=read("Fig_Donor_CellNumber").sort_values("n_cells"); colors=np.where(d.group.eq("Control"),BLUE,RED)
    fig,ax=plt.subplots(figsize=(5.1,3.8)); bars=ax.barh(np.arange(len(d)),d.n_cells,color=colors,alpha=.82,height=.66)
    ax.set_yticks(range(len(d)),d.donor_id,fontsize=7.7); ax.set_xlabel("Cells per donor")
    for b,v in zip(bars,d.n_cells): ax.text(v+max(d.n_cells)*.018,b.get_y()+b.get_height()/2,f"{v:,}",va="center",fontsize=7.1)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=BLUE,label="Control"),Patch(facecolor=RED,label="Sepsis")],frameon=False,ncol=2,fontsize=7.2,loc="lower right")
    ax.grid(axis="x",color="#E5EAED",lw=.6); ax.set_axisbelow(True); ax.set_xlim(0,max(d.n_cells)*1.18)
    fig.subplots_adjust(left=.18,right=.97,bottom=.16,top=.96); save(fig,"Fig_Donor_CellNumber")

def s3d():
    d=read("Fig_Donor_StateAbundance").sort_values("effect_sepsis_minus_control",ascending=False)
    mat=d[["effect_sepsis_minus_control"]].to_numpy(); fig,(ax,txt)=plt.subplots(1,2,figsize=(5.8,3.5),gridspec_kw={"width_ratios":[1,2.4]})
    im=ax.imshow(mat,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-.8,vcenter=0,vmax=.8))
    ax.set_xticks([0],["Sepsis − control"]); ax.set_yticks(range(len(d)),[x.replace(" / "," /\n") for x in d.state_label],fontsize=6.8)
    for i,r in enumerate(d.itertuples()): ax.text(0,i,f"{r.effect_sepsis_minus_control*100:+.1f}%",ha="center",va="center",fontsize=7.4,color="white" if abs(r.effect_sepsis_minus_control)>.45 else INK)
    for s in ax.spines.values(): s.set_visible(False)
    txt.set_xlim(0,1); txt.set_ylim(len(d)-.5,-.5); txt.axis("off"); txt.text(.02,-.65,"Bootstrap 95% CI",fontsize=7.4,fontweight="bold"); txt.text(.98,-.65,"FDR",fontsize=7.4,fontweight="bold",ha="right")
    for i,r in enumerate(d.itertuples()): txt.add_patch(plt.Rectangle((0,i-.46),1,.92,facecolor="#F4F6F6" if i%2==0 else "white",edgecolor="none")); txt.text(.02,i,f"[{r.ci_low*100:+.1f}, {r.ci_high*100:+.1f}]%",va="center",fontsize=6.9); txt.text(.98,i,f"{r.fdr:.3f}",ha="right",va="center",fontsize=6.9)
    fig.subplots_adjust(left=.38,right=.98,bottom=.18,top=.88,wspace=.36); save(fig,"Fig_Donor_StateAbundance")

def svg_pdf(path):
    if path.suffix.lower()==".pdf": return fitz.open(path)
    svg=fitz.open("svg",path.read_bytes()); return fitz.open("pdf",svg.convert_to_pdf())

def make_supp(stem,stems):
    out=fitz.open(); page=out.new_page(width=540,height=720); page.draw_rect(page.rect,color=None,fill=(1,1,1))
    xs=[(22,264),(278,520)]; ys=[(18,224),(252,458),(486,692)]
    for i,name in enumerate(stems):
        x0,x1=xs[i%2]; y0,y1=ys[i//2]; rect=fitz.Rect(x0,y0,x1,y1); src=svg_pdf(FIG/f"{name}.svg"); page.show_pdf_page(rect,src,0,keep_proportion=True); src.close()
        page.insert_text((x0-10,y0+8),chr(65+i),fontname="Helvetica-Bold",fontsize=10,color=(.08,.08,.08))
    pdf=OUT/f"{stem}_fig.pdf"; out.save(pdf,garbage=4,deflate=True,clean=True)
    scale=2250/page.rect.width; pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False); im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples); im.save(OUT/f"{stem}_fig.tif",dpi=(300,300),compression="tiff_lzw"); im.thumbnail((1300,1300)); im.save(FIG/f"{stem}_preview.png"); out.close()

def uppercase_labels_inplace(path):
    doc=fitz.open(path); page=doc[0]; hits=[]
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines",[]):
            for span in line["spans"]:
                t=span["text"].strip()
                if t in list("abcdef"):
                    r=fitz.Rect(span["bbox"]); hits.append((r,t.upper(),span["size"]))
                    rr=fitz.Rect(r.x0-1,r.y0-1,r.x1+1,r.y1+1); page.add_redact_annot(rr,fill=(1,1,1))
    if hits:
        page.apply_redactions()
        for r,t,size in hits: page.insert_text((r.x0,r.y1-1),t,fontname="Helvetica-Bold",fontsize=size,color=(.08,.08,.08),overlay=True)
        temp=path.with_name(path.stem+"_tmp_upper.pdf"); doc.save(temp,garbage=4,deflate=True,clean=True); doc.close(); temp.replace(path)
    else: doc.close()

def preserve_md(stem,filename):
    path=MD/filename; uppercase_labels_inplace(path); src=fitz.open(path); out=fitz.open(); page=out.new_page(width=540,height=540*src[0].rect.height/src[0].rect.width); page.show_pdf_page(page.rect,src,0,keep_proportion=True); src.close()
    pdf=OUT/f"{stem}_fig.pdf"; out.save(pdf,garbage=4,deflate=True,clean=True); scale=2250/page.rect.width; pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False); im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples); im.save(OUT/f"{stem}_fig.tif",dpi=(300,300),compression="tiff_lzw"); im.thumbnail((1300,1300)); im.save(FIG/f"{stem}_preview.png"); out.close()

def normalize_submission_tiffs():
    """Keep raster submissions within the PLOS 2625-pixel height limit."""
    for path in OUT.glob("*.tif"):
        with Image.open(path) as src:
            im=src.convert("RGB")
        if im.height>2625:
            width=round(im.width*2625/im.height)
            im=im.resize((width,2625),Image.Resampling.LANCZOS)
            im.save(path,dpi=(300,300),compression="tiff_lzw")

def _panel(ax, label):
    ax.text(-.16,1.08,label,transform=ax.transAxes,fontsize=12,fontweight="bold",
            ha="left",va="top",color=INK,clip_on=False)

def make_s2_direct():
    """Direct, source-data redraw: no pre-composed panels are nested in S2."""
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":9,
                         "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8})
    fig=plt.figure(figsize=(7.5,8.75),facecolor="white")
    gs=fig.add_gridspec(3,2,left=.15,right=.975,bottom=.07,top=.975,
                        hspace=.48,wspace=.38,height_ratios=[1.02,1.05,.93])

    # A: comparator performance matrix. Values are more legible than nine overlapping PR curves.
    ax=fig.add_subplot(gs[0,0]); d=read("Fig_ClinicalComparator_PR"); rows=[]
    for key,z in d.groupby("comparison_id",sort=False):
        z=z.sort_values("recall"); au=float(np.trapezoid(z.precision,z.recall)); base=float(z.prevalence.iloc[0])
        rows.append((key,au,base,au/base))
    z=pd.DataFrame(rows,columns=["comparison","AUPRC","baseline","lift"]).sort_values("lift",ascending=False)
    vals=z[["AUPRC","baseline","lift"]].copy(); score=vals.copy()
    score["lift"]=np.clip((score["lift"]-.85)/1.0,0,1)
    cmap=LinearSegmentedColormap.from_list("clinical",["#F3F5F5","#A8C3C5",TEAL])
    ax.imshow(score,aspect="auto",cmap=cmap,vmin=0,vmax=1)
    labels=[]
    for s in z.comparison:
        ds,grp=s.split("_",1); labels.append(ds+" | "+grp.replace("CardiogenicShock","Cardiogenic shock").replace("NonInfectiousICU","Non-infectious ICU").replace("ChronicDisease","Chronic disease"))
    labels=[s.replace("Cardiogenic shock","C. shock").replace("Non-infectious ICU","ICU ctrl").replace("NoninfectiousICU","ICU ctrl").replace("Chronic disease","Chronic").replace("Postoperative","Post-op") for s in labels]
    ax.set_yticks(range(len(z)),labels,fontsize=8); ax.set_xticks(range(3),["AUPRC","Baseline","Lift"])
    for i,r in z.reset_index(drop=True).iterrows():
        for j,t in enumerate([f"{r.AUPRC:.2f}",f"{r.baseline:.2f}",f"{r.lift:.2f}×"]):
            ax.text(j,i,t,ha="center",va="center",fontsize=8,color="white" if score.iloc[i,j]>.62 else INK)
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"A")

    # B: show every observation; medians and IQR remain visible without implying density in tiny groups.
    ax=fig.add_subplot(gs[0,1]); d=read("Fig_S100A8_APACHEII").copy()
    d["severity"]=pd.cut(d.apache_ii,[-np.inf,15,24,np.inf],labels=["≤15","16–24","≥25"])
    d["outcome"]=np.where(d.death==1,"Non-survivor","Survivor")
    rng=np.random.default_rng(21); positions=[]; data=[]; colors=[]; labels=[]
    k=0
    for sev in ["≤15","16–24","≥25"]:
        for out,c in [("Survivor",BLUE),("Non-survivor",RED)]:
            v=d[(d.severity.astype(str)==sev)&(d.outcome==out)].expression.to_numpy()
            if len(v):
                positions.append(k); data.append(v); colors.append(c); labels.append((sev,out,len(v)))
            k+=1
    bp=ax.boxplot(data,positions=positions,widths=.58,patch_artist=True,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.1},whiskerprops={"color":MUTED},capprops={"color":MUTED})
    for b,c in zip(bp["boxes"],colors): b.set(facecolor="white",edgecolor=c,linewidth=1.0)
    for pos,v,c in zip(positions,data,colors): ax.scatter(pos+rng.uniform(-.18,.18,len(v)),v,s=20,color=c,alpha=.72,edgecolor="white",linewidth=.35,zorder=3)
    ax.set_xticks([.5,2.5,4.5],["≤15","16–24","≥25"]); ax.set_xlabel("APACHE II category"); ax.set_ylabel("S100A8 expression")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=BLUE,label="Survivor"),Patch(facecolor=RED,label="Non-survivor")],frameon=False,ncol=2,loc="upper left")
    ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"B")

    # C: within-cohort z-scores make cross-study mortality contrasts comparable.
    ax=fig.add_subplot(gs[1,0]); d=read("Fig_S100A8_Mortality").copy(); datasets=list(d.dataset.drop_duplicates())
    positions=[]; arrays=[]; colors=[]
    for i,ds in enumerate(datasets):
        for j,(grp,c) in enumerate([("Survivor",BLUE),("Non-survivor",RED)]):
            positions.append(i*3+j); arrays.append(d[(d.dataset==ds)&(d.mortality_label==grp)].S100A8_z.to_numpy()); colors.append(c)
    vp=ax.violinplot(arrays,positions=positions,widths=.72,showextrema=False)
    for body,c in zip(vp["bodies"],colors): body.set_facecolor(c); body.set_edgecolor(c); body.set_alpha(.18)
    bp=ax.boxplot(arrays,positions=positions,widths=.30,patch_artist=True,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.1},whiskerprops={"color":MUTED},capprops={"color":MUTED})
    for b,c in zip(bp["boxes"],colors): b.set(facecolor="white",edgecolor=c,linewidth=1.0)
    rng=np.random.default_rng(8)
    for p,v,c in zip(positions,arrays,colors):
        if len(v)<=120: ax.scatter(p+rng.uniform(-.14,.14,len(v)),v,s=9,color=c,alpha=.38,edgecolors="none")
    ax.axhline(0,color="#9AA4AA",lw=.7); ax.set_xticks([.5,3.5,6.5],datasets); ax.set_ylabel("Within-cohort S100A8 z-score")
    ax.legend(handles=[Patch(facecolor=BLUE,label="Survivor"),Patch(facecolor=RED,label="Non-survivor")],frameon=False,ncol=2,loc="upper right")
    ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"C")

    # D: patient-level longitudinal heat map with a dedicated outcome side bar.
    ax=fig.add_subplot(gs[1,1]); d=read("Fig_S100A8_Longitudinal_Trajectory")
    p=d.pivot_table(index=["outcome","patient_id"],columns="day",values="expression",aggfunc="mean")
    p=p.reindex(sorted(p.index,key=lambda x:(x[0]!="Non-survivor",str(x[1])))); arr=p.to_numpy()
    cmap=LinearSegmentedColormap.from_list("trajectory",["#315B78","#F0EBD8","#D4A13A"]); cmap.set_bad("#ECEFF0")
    im=ax.imshow(np.ma.masked_invalid(arr),aspect="auto",cmap=cmap)
    ax.set_xticks(range(len(p.columns)),[str(x) for x in p.columns]); ax.set_xlabel("Day after enrollment")
    ticks=np.linspace(0,len(p)-1,min(7,len(p)),dtype=int); ax.set_yticks(ticks,[str(p.index[i][1]) for i in ticks]); ax.set_ylabel("Patient")
    split=sum(1 for x in p.index if x[0]=="Non-survivor")
    if 0<split<len(p): ax.axhline(split-.5,color="white",lw=1.3)
    ax.text(.985,.78,"Non-survivor",transform=ax.transAxes,rotation=90,va="center",ha="right",fontsize=8,color=RED,
            bbox=dict(facecolor="white",edgecolor="none",alpha=.72,pad=.5))
    ax.text(.985,.20,"Survivor",transform=ax.transAxes,rotation=90,va="center",ha="right",fontsize=8,color=BLUE,
            bbox=dict(facecolor="white",edgecolor="none",alpha=.72,pad=.5))
    ax.text(.5,1.015,"Low  ←  S100A8 expression  →  High",transform=ax.transAxes,ha="center",va="bottom",fontsize=8,color=MUTED)
    [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"D")

    # E: raw changes, median and IQR; omit a smoothed violin for the small subgroups.
    ax=fig.add_subplot(gs[2,0]); d=read("Fig_S100A8_TreatmentResponse"); groups=[("NR","Non-responder",RED),("R","Responder",BLUE)]
    arrays=[d[d.response==g].delta_T2_minus_T1.to_numpy() for g,_,_ in groups]
    bp=ax.boxplot(arrays,positions=[0,1],widths=.42,patch_artist=True,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.2},whiskerprops={"color":MUTED},capprops={"color":MUTED})
    for b,(_,_,c) in zip(bp["boxes"],groups): b.set(facecolor="white",edgecolor=c,linewidth=1.1)
    rng=np.random.default_rng(9)
    for i,(v,(_,_,c)) in enumerate(zip(arrays,groups)): ax.scatter(i+rng.uniform(-.17,.17,len(v)),v,s=28,color=c,alpha=.78,edgecolor="white",linewidth=.4,zorder=3)
    ax.axhline(0,color="#8E989E",lw=.8); ax.set_xticks([0,1],[x[1] for x in groups]); ax.set_ylabel("ΔS100A8 logCPM (T2 − T1)")
    ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"E")

    # F: conventional forest plot retains estimate, uncertainty and p values.
    fg=gs[2,1].subgridspec(1,2,width_ratios=[1.55,1.0],wspace=.05)
    ax=fig.add_subplot(fg[0,0]); tab=fig.add_subplot(fg[0,1],sharey=ax)
    d=read("Fig_Prognosis_Forest").copy().iloc[::-1].reset_index(drop=True)
    y=np.arange(len(d)); lo=d.estimate-d.ci_low; hi=d.ci_high-d.estimate
    sig=d.p_value<.05; cols=np.where(sig,RED,BLUE)
    for i,r in d.iterrows():
        ax.errorbar(r.estimate,i,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt="o",ms=5,
                    color=cols[i],ecolor=cols[i],elinewidth=1.1,capsize=2.5,mec="white",mew=.4)
    ax.axvline(1,color="#8D979D",lw=.8,ls=(0,(3,2))); ax.set_xscale("log")
    model_short={"unadjusted":"Unadj.","age_sex_adjusted":"Age/sex","age_adjusted":"Age"}
    labels=[f"{r.dataset} | {model_short.get(r.model,r.model)}" for r in d.itertuples()]
    ax.set_yticks(y,labels,fontsize=8); ax.set_xlabel("OR/HR per 1 SD")
    ax.set_xlim(.08,4.0); ax.set_xticks([.1,.5,1,2],["0.1","0.5","1","2"])
    tab.set_xlim(0,1); tab.set_ylim(-.5,len(d)-.5); tab.axis("off")
    tab.text(.02,1.02,"Estimate [95% CI]",transform=tab.transAxes,fontsize=8,fontweight="bold",ha="left")
    for i,r in d.iterrows():
        tab.text(.02,i+.13,f"{r.estimate:.2f} [{r.ci_low:.2f}–{r.ci_high:.2f}]",ha="left",va="center",fontsize=8)
        tab.text(.02,i-.17,f"p = {r.p_value:.3f}",ha="left",va="center",fontsize=8,color=MUTED)
    ax.grid(axis="x",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"F")

    master=OUT/"S2_fig_python_master.pdf"; fig.savefig(master,format="pdf")
    fig.savefig(FIG/"S2_preview.png",dpi=300); plt.close(fig)

def make_s3_direct():
    """Donor robustness and monocyte-state annotation from independent source tables."""
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":9,
                         "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8})
    fig=plt.figure(figsize=(7.5,8.75),facecolor="white")
    gs=fig.add_gridspec(3,2,left=.105,right=.975,bottom=.12,top=.975,
                        hspace=.48,wspace=.40,height_ratios=[1.05,.92,1.03])
    donors=["HC1","HC2","NSES","P25","P50","S2","S3"]
    pal=dict(zip(donors,["#4D7891","#77A6B6","#D4A13A","#C66B55","#8B77A8","#6AA889","#C487A9"]))

    ax=fig.add_subplot(gs[0,0]); d=read("Fig_Donor_UMAP")
    for donor in donors:
        z=d[d.donor_id==donor]
        if len(z)>4000: z=z.sample(4000,random_state=17)
        ax.scatter(z.umap_1,z.umap_2,s=1.5,color=pal[donor],alpha=.42,edgecolors="none",label=donor,rasterized=True)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(frameon=False,ncol=4,markerscale=4,loc="upper center",bbox_to_anchor=(.52,1.02),handletextpad=.15,columnspacing=.65)
    [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"A")

    ax=fig.add_subplot(gs[0,1]); d=read("Fig_Donor_Composition")
    cells=["T Cells","B Cells","NK Cells","Monocytes"]
    mat=d.pivot(index="donor_id",columns="cell_type",values="proportion").reindex(index=donors,columns=cells)
    cmap=LinearSegmentedColormap.from_list("composition",["#F3F5F5","#A8C3C5",TEAL])
    ax.imshow(mat,aspect="auto",cmap=cmap,vmin=0,vmax=.65)
    ax.set_xticks(range(4),["T cells","B cells","NK cells","Monocytes"],rotation=24,ha="right")
    ax.set_yticks(range(len(donors)),donors)
    for i in range(len(donors)):
        for j in range(4):
            v=mat.iloc[i,j]; ax.text(j,i,f"{v*100:.1f}%",ha="center",va="center",fontsize=8,color="white" if v>.38 else INK)
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"B")

    ax=fig.add_subplot(gs[1,0]); d=read("Fig_Donor_CellNumber").sort_values("n_cells")
    cols=[BLUE if g=="Control" else RED for g in d.group]; bars=ax.barh(np.arange(len(d)),d.n_cells,color=cols,height=.64,alpha=.88)
    ax.set_yticks(range(len(d)),d.donor_id); ax.set_xlabel("Cells per donor"); ax.set_xlim(0,d.n_cells.max()*1.18)
    for b,v in zip(bars,d.n_cells): ax.text(v+d.n_cells.max()*.018,b.get_y()+b.get_height()/2,f"{v:,}",va="center",fontsize=8)
    ax.legend(handles=[Patch(facecolor=BLUE,label="Control"),Patch(facecolor=RED,label="Sepsis")],frameon=False,ncol=2,loc="lower right")
    ax.grid(axis="x",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"C")

    d=read("Fig_Donor_StateAbundance").sort_values("effect_sepsis_minus_control").reset_index(drop=True)
    dg=gs[1,1].subgridspec(1,2,width_ratios=[1.15,2.35],wspace=.03)
    lab=fig.add_subplot(dg[0,0]); ax=fig.add_subplot(dg[0,1],sharey=lab)
    y=np.arange(len(d)); colors=np.where(d.effect_sepsis_minus_control>0,RED,BLUE)
    for i,r in d.iterrows():
        ax.errorbar(r.effect_sepsis_minus_control,i,xerr=[[r.effect_sepsis_minus_control-r.ci_low],[r.ci_high-r.effect_sepsis_minus_control]],
                    fmt="o",ms=7,color=colors[i],ecolor=colors[i],elinewidth=1.4,capsize=3,mec="white",mew=.5)
        ax.text(.98,i-.16 if i else i+.16,f"FDR = {r.fdr:.3f}",transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=8,color=MUTED)
    state_labels=["VCAN+ classical\nmonocyte", "Emergency myeloid /\nimmature-like"]
    lab.set_xlim(0,1); lab.set_ylim(-.5,len(d)-.5); lab.axis("off")
    for yy,label in zip(y,state_labels):
        lab.text(.98,yy,label,ha="right",va="center",fontsize=8)
    ax.set_yticks(y,[])
    ax.set_xlim(min(-.70,float(d.ci_low.min())-.06),max(.86,float(d.ci_high.max())+.06))
    ax.axvline(0,color="#929DA3",lw=.8,ls=(0,(3,2))); ax.set_xlabel("Sepsis − control proportion")
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1,decimals=0)); ax.grid(axis="x",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(lab,"D")

    ax=fig.add_subplot(gs[2,0]); d=read("Fig_Monocyte_UMAP_ByDonor")
    for donor in donors:
        z=d[d.donor_id==donor]; ax.scatter(z.UMAP_1,z.UMAP_2,s=5,color=pal[donor],alpha=.65,edgecolors="none",label=donor,rasterized=True)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
    # Donor colors are shared with panel A; a second legend would obscure the sparse monocyte map.
    [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"E")

    # F: exploratory donor-level TF activity. None of these contrasts passes FDR.
    ax=fig.add_subplot(gs[2,1])
    stats=pd.read_csv(FIG/"Fig_DonorLevel_TFActivity_source_data.csv").sort_values("effect")
    y=np.arange(len(stats)); colors=np.where(stats.effect>=0,RED,BLUE)
    ax.barh(y,stats.effect,color=colors,alpha=.82,height=.62)
    ax.axvline(0,color="#8D979D",lw=.8)
    ax.set_xlim(stats.effect.min()-0.22, stats.effect.max()+0.32)
    ax.set_yticks(y,stats.tf); ax.set_xlabel("Sepsis − control effect")
    ax.grid(axis="x",color="#E7EBED",lw=.6); ax.set_axisbelow(True)
    for yy,r in enumerate(stats.itertuples()):
        if r.effect < -0.45:
            ax.text(r.effect+.055,yy,f"{r.effect:+.2f}",va="center",ha="left",
                    fontsize=8,color="white")
        else:
            ax.text(r.effect+(.035 if r.effect>=0 else -.035),yy,f"{r.effect:+.2f}",va="center",
                    ha="left" if r.effect>=0 else "right",fontsize=8)
    ax.text(.5,1.02,"All FDR > 0.05 (range 0.286–0.905)",transform=ax.transAxes,
            ha="center",va="bottom",fontsize=8,color=MUTED)
    _panel(ax,"F")
    master=OUT/"S3_fig_python_master.pdf"
    fig.savefig(master,format="pdf"); fig.savefig(FIG/"S3_preview.png",dpi=300)
    plt.close(fig)
    for suffix in ("pdf","ai"):
        (OUT/f"S3_fig.{suffix}").write_bytes(master.read_bytes())
    doc=fitz.open(master); page=doc[0]; scale=2250/page.rect.width
    pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    if im.height>2625:
        im=im.resize((round(im.width*2625/im.height),2625),Image.Resampling.LANCZOS)
    im.save(OUT/"S3_fig.tif",dpi=(300,300),compression="tiff_lzw"); doc.close()

def make_s4_direct():
    """Independent-cohort validation without nested facets or duplicated mini-figures."""
    from scipy.stats import spearmanr
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":9,
                         "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8})
    fig=plt.figure(figsize=(7.5,8.75),facecolor="white")
    gs=fig.add_gridspec(3,2,left=.11,right=.975,bottom=.075,top=.975,hspace=.48,wspace=.42)

    ax=fig.add_subplot(gs[0,0]); d=read("Fig_Pseudobulk_S100A8_Donor")
    groups=[("Control",BLUE),("Sepsis",RED)]; rng=np.random.default_rng(5)
    arrays=[d[d.group==g].logCPM.to_numpy() for g,_ in groups]
    bp=ax.boxplot(arrays,positions=[0,1],widths=.42,patch_artist=True,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.2},whiskerprops={"color":MUTED},capprops={"color":MUTED})
    for b,(_,c) in zip(bp["boxes"],groups): b.set(facecolor="white",edgecolor=c,linewidth=1.1)
    for i,((g,c),v) in enumerate(zip(groups,arrays)):
        z=d[d.group==g]
        ax.scatter(i+rng.uniform(-.12,.12,len(v)),v,s=34,color=c,edgecolor="white",linewidth=.5,zorder=3)
        for xx,yy,label in zip(i+rng.uniform(-.12,.12,len(v)),v,z.donor_id): ax.text(xx,yy+.08,label,ha="center",va="bottom",fontsize=8,color=c)
    ax.set_xticks([0,1],["Control","Sepsis"]); ax.set_ylabel("S100A8 donor pseudo-bulk logCPM"); ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"A")

    ax=fig.add_subplot(gs[0,1]); d=read("Fig_IndependentSignature_DonorScores")
    sigs=list(d.signature.drop_duplicates()); x=np.arange(len(sigs)); rng=np.random.default_rng(12)
    for off,(grp,c) in zip([-.18,.18],groups):
        for i,sig in enumerate(sigs):
            v=d[(d.group==grp)&(d.signature==sig)].score.to_numpy()
            ax.scatter(np.full(len(v),i+off)+rng.uniform(-.055,.055,len(v)),v,s=20,color=c,alpha=.78,edgecolor="white",linewidth=.35)
            med=np.median(v); ax.plot([i+off-.10,i+off+.10],[med,med],color=c,lw=2)
    short=["Inflam. −A8","Inflam. −A8/A9","Antigen pres.","Combined −A8","Combined −A8/A9"]
    ax.set_xticks(x,short,rotation=28,ha="right"); ax.set_ylabel("Donor-level score")
    ax.legend(handles=[Patch(facecolor=BLUE,label="Control"),Patch(facecolor=RED,label="Sepsis")],frameon=False,ncol=2,loc="upper left")
    ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"B")

    ax=fig.add_subplot(gs[1,0]); d=read("Fig_Bulk_IndependentSignature")
    datasets=["GSE131411","GSE134347","GSE28750","GSE65682"]
    cats=["Healthy","Cardiogenic shock","Non-infectious ICU","Postoperative sterile inflammation","Sepsis"]
    catlab=["Healthy","C. shock","ICU control","Post-op","Sepsis"]
    mat=np.full((4,5),np.nan); nn=np.zeros((4,5),int)
    for i,ds in enumerate(datasets):
        for j,g in enumerate(cats):
            v=d[(d.dataset==ds)&(d.group==g)].independent_score
            if len(v): mat[i,j]=v.median(); nn[i,j]=len(v)
    lim=max(.5,float(np.nanmax(np.abs(mat)))); im=ax.imshow(mat,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    ax.set_yticks(range(4),datasets); ax.set_xticks(range(5),catlab,rotation=28,ha="right")
    for i in range(4):
        for j in range(5):
            if not np.isnan(mat[i,j]): ax.text(j,i,f"{mat[i,j]:+.2f}\nn={nn[i,j]}",ha="center",va="center",fontsize=8,color="white" if abs(mat[i,j])>.62*lim else INK)
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]; _panel(ax,"C")

    ax=fig.add_subplot(gs[1,1]); d=read("Fig_SCP548_S100A8_HLAII_Correlation")
    cohorts=list(d.cohort.drop_duplicates()); palette=dict(zip(cohorts,["#4D7891","#77A6B6","#D4A13A","#B67845","#8B77A8","#C66B55","#9B4F43"]))
    for g in cohorts:
        z=d[d.cohort==g]; ax.scatter(z.S100A8,z.antigen_presentation_score,s=28,color=palette[g],alpha=.82,edgecolor="white",linewidth=.35,label=g)
    rho,p=spearmanr(d.S100A8,d.antigen_presentation_score,nan_policy="omit"); ax.text(.03,.96,f"Spearman ρ = {rho:.2f}\np = {p:.3g}",transform=ax.transAxes,ha="left",va="top",fontsize=8)
    ax.set_xlabel("S100A8 donor logCPM"); ax.set_ylabel("Antigen-presentation score"); ax.legend(frameon=False,ncol=2,loc="lower left",fontsize=8)
    ax.grid(color="#E7EBED",lw=.5); ax.set_axisbelow(True); _panel(ax,"D")

    ax=fig.add_subplot(gs[2,0]); d=read("Fig_GSE205672_S100A8_HLAII_Correlation"); cols={"Control":BLUE,"Sepsis":RED}
    for g in ["Control","Sepsis"]:
        z=d[d.group==g]; ax.scatter(z.S100A8,z.antigen_presentation_score,s=30,color=cols[g],alpha=.82,edgecolor="white",linewidth=.4,label=g)
    rho,p=spearmanr(d.S100A8,d.antigen_presentation_score,nan_policy="omit"); ax.text(.03,.96,f"Spearman ρ = {rho:.2f}\np = {p:.3g}",transform=ax.transAxes,ha="left",va="top",fontsize=8)
    ax.set_xlabel("log2(S100A8 FPKM + 1)"); ax.set_ylabel("Antigen-presentation score"); ax.legend(frameon=False,loc="lower left")
    ax.grid(color="#E7EBED",lw=.5); ax.set_axisbelow(True); _panel(ax,"E")

    ax=fig.add_subplot(gs[2,1]); d=read("Fig_GSE205672_NFkBScore"); arrays=[d[d.group==g].NFkB_score.to_numpy() for g,_ in groups]
    bp=ax.boxplot(arrays,positions=[0,1],widths=.42,patch_artist=True,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.2},whiskerprops={"color":MUTED},capprops={"color":MUTED})
    for b,(_,c) in zip(bp["boxes"],groups): b.set(facecolor="white",edgecolor=c,linewidth=1.1)
    rng=np.random.default_rng(18)
    for i,(v,(_,c)) in enumerate(zip(arrays,groups)): ax.scatter(i+rng.uniform(-.16,.16,len(v)),v,s=24,color=c,alpha=.78,edgecolor="white",linewidth=.35)
    ax.set_xticks([0,1],["Control","Sepsis"]); ax.set_ylabel("NF-κB transcriptional score"); ax.grid(axis="y",color="#E7EBED",lw=.6); ax.set_axisbelow(True); _panel(ax,"F")
    fig.savefig(OUT/"S4_fig_python_master.pdf",format="pdf"); fig.savefig(FIG/"S4_preview.png",dpi=300); plt.close(fig)

def main():
    make_s2_direct(); make_s3_direct(); make_s4_direct(); s3b(); s3c(); s3d()
    preserve_md("S5","Figure_1_MD_system_metrics.pdf"); preserve_md("S6","Figure_2_ligand_behavior.pdf")
    normalize_submission_tiffs()

if __name__=="__main__": main()
