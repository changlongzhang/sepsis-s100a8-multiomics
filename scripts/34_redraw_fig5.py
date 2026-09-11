"""Redraw Fig. 5 source panels with distinct evidence-appropriate graphics."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from scipy.stats import gaussian_kde
import numpy as np
import pandas as pd

ROOT = Path(r"D:\桌面\sepsis\review_revision\09_figures")
INK, MUTED = "#27333D", "#65727E"
BLUE, ORANGE, RED, PURPLE, TEAL = "#4C78A8", "#D49A45", "#D66A5E", "#7A7180", "#4F938D"
COHORTS = ["Control", "Leuk-UTI", "URO", "Int-URO", "ICU-NoSEP", "ICU-SEP", "Bac-SEP"]
COHORT_COLORS = dict(zip(COHORTS, [BLUE, "#7F9CB7", ORANGE, "#D9AE70", "#8793A2", RED, "#B85650"]))

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
})


def data(stem):
    return pd.read_csv(ROOT / f"{stem}_source_data.csv")


def save(fig, stem):
    target = {
        "Fig_Monocyte_UnsupervisedClusters": 17.0,
        "Fig_Monocyte_StatePrograms": 15.0,
        "Fig_SCP548_DonorPseudobulk_S100A8": 16.5,
        "Fig_SCP548_MS1_Abundance": 16.5,
        "Fig_SCP548_AntigenPresentation": 15.0,
        "Fig_GSE205672_S100A8_Monocytes": 16.0,
        "Fig_GSE205672_AntigenPresentation": 16.0,
        "Fig_Mechanistic_Triangulation": 16.5,
    }[stem]
    for item in fig.findobj(match=matplotlib.text.Text):
        if item.get_text().strip() and item.get_fontsize() < target:
            item.set_fontsize(target)
    try:
        fig.tight_layout(pad=0.45)
    except Exception:
        pass
    base = ROOT / stem
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base.with_suffix(".png"), dpi=300)
    plt.close(fig)


def boot_ci(values, stat=np.median, seed=20260824):
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    boot = np.array([stat(rng.choice(values, len(values), replace=True)) for _ in range(2500)])
    return stat(values), *np.quantile(boot, [0.025, 0.975])


def panel_a():
    d = data("Fig_Monocyte_UnsupervisedClusters")
    colors = {0: BLUE, 1: ORANGE, 2: RED}
    labels = {0: "C0", 1: "C1", 2: "C2"}
    fig, ax = plt.subplots(figsize=(5.4, 4.1))
    for cluster in [0, 1, 2]:
        z = d[d.cluster == cluster]
        x, y = z.UMAP_1.to_numpy(), z.UMAP_2.to_numpy()
        if len(z) > 25:
            mx, my = np.median(x), np.median(y)
            sx = max(np.median(np.abs(x - mx)) * 1.4826, 0.25)
            sy = max(np.median(np.abs(y - my)) * 1.4826, 0.25)
            core = np.sqrt(((x - mx) / sx) ** 2 + ((y - my) / sy) ** 2) < 3.2
            xc, yc = x[core], y[core]
            xx, yy = np.mgrid[xc.min()-0.35:xc.max()+0.35:80j, yc.min()-0.35:yc.max()+0.35:80j]
            dens = gaussian_kde(np.vstack([xc, yc]))(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
            level = np.quantile(dens[dens > 0], 0.48)
            ax.contourf(xx, yy, dens, levels=[level, dens.max()], colors=[colors[cluster]], alpha=0.10)
            ax.contour(xx, yy, dens, levels=[level], colors=[colors[cluster]], linewidths=1.1, alpha=0.75)
        ax.scatter(x, y, s=9, color=colors[cluster], alpha=0.55, edgecolor="white", linewidth=0.2,
                   label=f"{labels[cluster]}  (n={len(z)})")
        mx, my = np.median(x), np.median(y)
        ax.text(mx, my, labels[cluster], ha="center", va="center", color=colors[cluster],
                fontsize=9.2, fontweight="bold", zorder=6,
                bbox={"boxstyle":"round,pad=0.24", "facecolor":"white",
                      "edgecolor":colors[cluster], "linewidth":0.9, "alpha":0.92})
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    # Extra horizontal breathing room prevents the terminal 10.0 tick from
    # touching the half-width panel boundary.
    ax.margins(x=.08, y=.08)
    fig.subplots_adjust(left=.16, right=.97, bottom=.17, top=.97)
    save(fig, "Fig_Monocyte_UnsupervisedClusters")


def panel_b():
    d = data("Fig_Monocyte_StatePrograms")
    programs = list(dict.fromkeys(d.program))
    clusters = list(dict.fromkeys(d.cluster_label))
    # Vertical layout: programs read top-to-bottom and the three states form
    # compact columns, which is clearer in a portrait half-width panel.
    mat = d.pivot(index="program", columns="cluster_label", values="mean_z_score").reindex(index=programs, columns=clusters)
    fig, ax = plt.subplots(figsize=(4.4, 5.2))
    cmap = LinearSegmentedColormap.from_list("program", [BLUE, "#EEF1F2", RED])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-1.7, vcenter=0, vmax=1.7))
    short_program = {
        "inflammatory": "Inflammatory", "antigen_presentation": "Antigen pres.",
        "immunosuppression": "Immunosuppress.", "immature_emergency": "Emergency",
        "monocyte_core": "Monocyte core", "classical": "Classical",
        "nonclassical": "Non-classical", "lymphoid_contamination": "Lymphoid",
        "combined_inflammatory_HLADRlow": "Inflam./HLA-DRlow",
    }
    col_labels = ["C0\nEmergency", "C1\nClassical", "C2\nLymphoid"]
    ax.set_xticks(range(len(clusters)), col_labels[:len(clusters)], fontsize=7.7)
    ax.set_yticks(range(len(programs)), [short_program[p] for p in programs], fontsize=7.5)
    for j in range(mat.shape[1]):
        key = np.nanargmax(np.abs(mat.iloc[:, j].to_numpy()))
        ax.scatter(j, key, s=115, facecolor="none", edgecolor=INK, linewidth=1.0)
    ax.set_xticks(np.arange(-.5, len(clusters), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(programs), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2); ax.tick_params(which="minor", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.025, pad=.025); cb.set_label("z", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    fig.subplots_adjust(left=.40, right=.90, bottom=.14, top=.97)
    save(fig, "Fig_Monocyte_StatePrograms")


def panel_c():
    d=data("Fig_SCP548_DonorPseudobulk_S100A8"); fig,ax=plt.subplots(figsize=(7.1,3.4))
    vals=[d.loc[d.cohort==c,"S100A8"].dropna().to_numpy(float) for c in COHORTS]
    bp=ax.boxplot(vals,patch_artist=True,widths=.62,showfliers=False,
                  medianprops={"color":INK,"linewidth":1.2},
                  whiskerprops={"color":"#7D8991","linewidth":.8},
                  capprops={"color":"#7D8991","linewidth":.8})
    for box,c in zip(bp["boxes"],COHORTS): box.set(facecolor=COHORT_COLORS[c],alpha=.28,edgecolor=COHORT_COLORS[c],linewidth=1.1)
    rng=np.random.default_rng(21)
    for i,(c,v) in enumerate(zip(COHORTS,vals),1):
        ax.scatter(i+rng.uniform(-.13,.13,len(v)),v,s=13,color=COHORT_COLORS[c],alpha=.55,edgecolor="white",linewidth=.25)
    ax.set_xticks(range(1,len(COHORTS)+1),COHORTS,rotation=28,ha="right",fontsize=7.7)
    ax.set_ylabel("S100A8 donor pseudobulk logCPM"); ax.grid(axis="y",color="#E7EBEE",lw=.6); ax.set_axisbelow(True)
    fig.subplots_adjust(left=.12,right=.98,bottom=.25,top=.97); save(fig,"Fig_SCP548_DonorPseudobulk_S100A8")


def panel_d():
    d=data("Fig_SCP548_MS1_Abundance"); bins=np.linspace(0,1,11); mat=[]
    for cohort in COHORTS:
        v=d.loc[d.cohort==cohort,"MS1_fraction"].dropna().to_numpy(float)
        h,_=np.histogram(v,bins=bins); mat.append(h/max(h.sum(),1))
    mat=np.asarray(mat); fig,ax=plt.subplots(figsize=(7.1,3.35))
    cmap=LinearSegmentedColormap.from_list("density",["#F5F7F7","#A9C1CB",BLUE])
    im=ax.imshow(mat,aspect="auto",cmap=cmap,vmin=0,vmax=max(.35,mat.max()))
    ax.set_yticks(range(len(COHORTS)),COHORTS,fontsize=7.8)
    ax.set_xticks(np.arange(10),[f"{(bins[i]+bins[i+1])/2:.2f}" for i in range(10)],rotation=35,ha="right",fontsize=7)
    ax.set_xlabel("MS1-like monocyte fraction per donor (bin midpoint)")
    for y in range(mat.shape[0]):
        for x in range(mat.shape[1]):
            if mat[y,x]>0: ax.text(x,y,f"{mat[y,x]*100:.0f}%",ha="center",va="center",fontsize=6.4,color="white" if mat[y,x]>.24 else INK)
    for s in ax.spines.values(): s.set_visible(False)
    # Cell percentages carry the scale directly; a second colorbar is redundant
    # at final half-width publication size.
    fig.subplots_adjust(left=.18,right=.94,bottom=.25,top=.97); save(fig,"Fig_SCP548_MS1_Abundance")


def panel_e():
    d=data("Fig_SCP548_AntigenPresentation"); ref=d.loc[d.cohort=="Control","antigen_presentation_score"].dropna().to_numpy(float)
    fig,ax=plt.subplots(figsize=(6.4,3.5)); rng=np.random.default_rng(17)
    rows=[]
    for cohort in COHORTS[1:]:
        v=d.loc[d.cohort==cohort,"antigen_presentation_score"].dropna().to_numpy(float)
        boot=np.array([np.median(rng.choice(v,len(v),True))-np.median(rng.choice(ref,len(ref),True)) for _ in range(3000)])
        rows.append((cohort,np.median(v)-np.median(ref),*np.quantile(boot,[.025,.975])))
    rr=rows[::-1]; y=np.arange(len(rr)); est=np.array([r[1] for r in rr]); colors=[COHORT_COLORS[r[0]] for r in rr]
    ax.barh(y,est,color=colors,alpha=.78,height=.64,edgecolor="white")
    for yy,r in enumerate(rr):
        xpos = r[1] + .04 if r[1] >= 0 else -.04
        ax.text(xpos,yy,f"{r[1]:+.2f}",va="center",ha="left" if r[1]>=0 else "right",fontsize=7.2,color=INK)
    ax.axvline(0,color="#79858E",lw=.8); ax.set_yticks(y,[r[0] for r in rr],fontsize=7.8)
    ax.set_xlabel("Median difference from Control in antigen-presentation score")
    ax.grid(axis="x",color="#E7EBEE",lw=.65); ax.set_axisbelow(True); fig.subplots_adjust(left=.18,right=.98,bottom=.20,top=.95)
    save(fig,"Fig_SCP548_AntigenPresentation")


def panel_f():
    d=data("Fig_GSE205672_S100A8_Monocytes"); groups=["Control","Sepsis"]; colors=[BLUE,RED]
    fig,ax=plt.subplots(figsize=(4.6,3.8)); rng=np.random.default_rng(6)
    for i,(g,c) in enumerate(zip(groups,colors)):
        v=d.loc[d.group==g,"S100A8"].dropna().to_numpy(float)
        parts=ax.violinplot(v,[i],widths=.72,showextrema=False); [b.set_facecolor(c) or b.set_alpha(.18) for b in parts["bodies"]]
        ax.scatter(i+rng.uniform(-.13,.13,len(v)),v,s=18,color=c,alpha=.68,edgecolor="white",linewidth=.35)
        med,lo,hi=boot_ci(v,seed=50+i); ax.plot([i,i],[lo,hi],color=INK,lw=1.3); ax.scatter(i,med,s=45,facecolor="white",edgecolor=c,linewidth=1.6,zorder=5)
    ax.set_xticks([0,1],groups); ax.set_ylabel("log2(S100A8 FPKM + 1)"); ax.grid(axis="y",color="#E7EBEE",lw=.65); ax.set_axisbelow(True)
    fig.subplots_adjust(left=.18,right=.97,bottom=.14,top=.96); save(fig,"Fig_GSE205672_S100A8_Monocytes")


def panel_g():
    d=data("Fig_GSE205672_AntigenPresentation"); fig,ax=plt.subplots(figsize=(4.8,3.8))
    allv=d.antigen_presentation_score.dropna().to_numpy(float); xx=np.linspace(allv.min()-.2,allv.max()+.2,240)
    for base,(g,c) in enumerate([("Control",BLUE),("Sepsis",RED)]):
        v=d.loc[d.group==g,"antigen_presentation_score"].dropna().to_numpy(float); yy=gaussian_kde(v)(xx); yy=yy/yy.max()*.72
        ax.fill_between(xx,base,base+yy,color=c,alpha=.28); ax.plot(xx,base+yy,color=c,lw=1.25)
        ax.text(xx.min(),base+.08,f"{g}  n={len(v)}",fontsize=7.7,color=c,va="bottom")
    ax.set_yticks([]); ax.set_xlabel("Antigen-presentation score"); ax.set_ylim(-.05,1.9)
    ax.grid(axis="x",color="#E7EBEE",lw=.6); ax.set_axisbelow(True); fig.subplots_adjust(left=.10,right=.97,bottom=.16,top=.96)
    save(fig,"Fig_GSE205672_AntigenPresentation")


def panel_h():
    d=data("Fig_Mechanistic_Triangulation").sort_values("rho"); fig,ax=plt.subplots(figsize=(7.1,3.0))
    corr_cmap=LinearSegmentedColormap.from_list("corr",[BLUE,"#F2F3F2",RED])
    mat=d[["rho"]].to_numpy(); im=ax.imshow(mat,aspect="auto",cmap=corr_cmap,norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1))
    short_labels = d.label.str.replace("Antigen-presentation program", "Ag-presentation", regex=False).str.replace("Emergency-myeloid program", "Emergency-myeloid", regex=False)
    ax.set_xticks([0],["Spearman ρ"],fontsize=7.7); ax.set_yticks(range(len(d)),short_labels,fontsize=7.4)
    for y,r in enumerate(d.itertuples()):
        ax.text(0,y,f"{r.rho:+.2f}\n[{r.ci_low:+.2f}, {r.ci_high:+.2f}]",ha="center",va="center",fontsize=6.7,color="white" if abs(r.rho)>.55 else INK)
    for s in ax.spines.values(): s.set_visible(False)
    # Every cell is numerically annotated, so omit the redundant colorbar.
    fig.subplots_adjust(left=.48,right=.88,bottom=.20,top=.97); save(fig,"Fig_Mechanistic_Triangulation")


FUNCS={"A":panel_a,"B":panel_b,"C":panel_c,"D":panel_d,"E":panel_e,"F":panel_f,"G":panel_g,"H":panel_h}
if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--panel",choices=list(FUNCS)+["all"],default="A"); a=p.parse_args()
    for k in (FUNCS if a.panel=="all" else [a.panel]): FUNCS[k](); print("rendered",k)
