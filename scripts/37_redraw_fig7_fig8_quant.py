"""Regenerate Fig. 7/8 quantitative panels with audited significance labels."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(r"D:\桌面\sepsis")
SRC=ROOT/"实验"/"PLOS_submission_ready_final"/"Fig7_source_data.csv"
OUT=ROOT/"review_revision"/"09_figures"; D=pd.read_csv(SRC)
CON=pd.read_csv(OUT/"Fig7_Fig8_planned_contrasts.csv")
INK="#27333D"; COL={"Control":"#4C78A8","LPS":"#D66A5E","LPS+TPL":"#4F938D","TPL":"#7FB5AE","siNC":"#4C78A8","siS100A8":"#C97832"}
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"pdf.fonttype":42,"svg.fonttype":"none","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.75,"text.color":INK,"axes.labelcolor":INK,"xtick.color":INK,"ytick.color":INK})

def ptxt(p):
    if p < .0001: return "****"
    if p < .001: return "***"
    if p < .01: return "**"
    if p < .05: return "*"
    return "ns"
def pval(panel, text): return float(CON[(CON.panel==panel)&(CON.contrast==text)].p_adjusted.iloc[0])
def bracket(ax,x1,x2,y,h,text):
    ax.plot([x1,x1,x2,x2],[y,y+h,y+h,y],color=INK,lw=.75,clip_on=False); ax.text((x1+x2)/2,y+h,text,ha="center",va="bottom",fontsize=6.8)
def save(fig,stem):
    target={"Fig7_qPCR":10.0,"Fig7_WB_quant":13.0,"Fig8_pNFkB":16.0,
            "Fig8_S100A8_8group":16.0,"Fig8_CD74":17.0,"Fig8_HLADRA":17.0,
            "Fig8_S100A8_4group":17.0,"Fig8_interaction_effect":13.0}[stem]
    for item in fig.findobj(match=matplotlib.text.Text):
        if item.get_text().strip() and item.get_fontsize()<target: item.set_fontsize(target)
    try: fig.tight_layout(pad=.4)
    except Exception: pass
    for ext in ["pdf","svg"]: fig.savefig(OUT/f"{stem}.{ext}")
    fig.savefig(OUT/f"{stem}.tiff",dpi=600,pil_kwargs={"compression":"tiff_lzw"}); fig.savefig(OUT/f"{stem}.png",dpi=300); plt.close(fig)
def draw_groups(ax,z,groups,colors,ylabel,brackets):
    rng=np.random.default_rng(42)
    means=[]; sems=[]
    for i,g in enumerate(groups):
        v=z[z.group==g].value.dropna().to_numpy(float); means.append(v.mean()); sems.append(v.std(ddof=1)/np.sqrt(len(v)))
        ax.bar(i,means[-1],width=.62,color=colors[i],alpha=.82,edgecolor="none"); ax.errorbar(i,means[-1],yerr=sems[-1],fmt="none",ecolor=INK,capsize=3,lw=.9)
        ax.scatter(i+rng.uniform(-.09,.09,len(v)),v,s=16,color=colors[i],alpha=.86,edgecolor="white",linewidth=.3,zorder=3)
    ymax=max(np.array(means)+np.array(sems)); yr=max(ymax,1e-6); step=yr*.16
    for j,(x1,x2,label) in enumerate(brackets): bracket(ax,x1,x2,ymax+step*(j+.18),step*.12,label)
    ax.set_ylim(bottom=min(0,ax.get_ylim()[0]),top=ymax+step*(len(brackets)+.75)); ax.set_xticks(range(len(groups)),groups,rotation=30,ha="right",fontsize=7.3); ax.set_ylabel(ylabel,fontsize=8)

def fig7_qpcr():
    panels=[("B","IL-6"),("C","TNF-α"),("D","IL-1β"),("E","S100A8")]; fig,axs=plt.subplots(1,4,figsize=(8.0,2.75))
    for ax,(p,title) in zip(axs,panels):
        z=D[D.panel==p]; bs=[(0,1,ptxt(pval(p,"Control vs LPS"))),(1,2,ptxt(pval(p,"LPS vs LPS+TPL")))]
        draw_groups(ax,z,["Control","LPS","LPS+TPL"],[COL["Control"],COL["LPS"],COL["LPS+TPL"]],"Relative mRNA expression",bs)
        ax.set_xlabel(title, fontweight="bold", fontsize=8)
        if ax is not axs[0]: ax.set_ylabel("")
    fig.subplots_adjust(left=.08,right=.99,bottom=.25,top=.80,wspace=.48); save(fig,"Fig7_qPCR")
def fig7_wb():
    panels=[("G","S100A8 / β-actin"),("H","NF-κB / β-actin"),("I","p-NF-κB / NF-κB")]; fig,axs=plt.subplots(1,3,figsize=(6.7,2.9))
    for ax,(p,title) in zip(axs,panels):
        z=D[D.panel==p]; bs=[(0,1,ptxt(pval(p,"Control vs LPS"))),(1,2,ptxt(pval(p,"LPS vs LPS+TPL")))]
        draw_groups(ax,z,["Control","LPS","LPS+TPL"],[COL["Control"],COL["LPS"],COL["LPS+TPL"]],"Normalized band intensity",bs)
        ax.set_xlabel(title, fontweight="bold", fontsize=7.5)
        if ax is not axs[0]: ax.set_ylabel("")
    fig.subplots_adjust(left=.09,right=.99,bottom=.25,top=.78,wspace=.42); save(fig,"Fig7_WB_quant")
def kd_plot(panel,stem,title,eight):
    z=D[D.panel==panel].copy(); order=[("siNC","Control"),("siNC","TPL"),("siNC","LPS"),("siNC","LPS+TPL"),("siS100A8","Control"),("siS100A8","TPL"),("siS100A8","LPS"),("siS100A8","LPS+TPL")] if eight else [("siNC","Control"),("siNC","LPS"),("siS100A8","Control"),("siS100A8","LPS")]
    fig,ax=plt.subplots(figsize=(5.1 if eight else 3.8,3.5)); rng=np.random.default_rng(7); means=[];sems=[]
    for i,(si,tr) in enumerate(order):
        v=z[(z.siRNA==si)&(z.treatment==tr)].value.to_numpy(float); means.append(v.mean()); sems.append(v.std(ddof=1)/np.sqrt(len(v))); c=COL[si]
        ax.bar(i,means[-1],.62,color=c,alpha=.72); ax.errorbar(i,means[-1],yerr=sems[-1],fmt="none",ecolor=INK,capsize=3,lw=.9); ax.scatter(i+rng.uniform(-.07,.07,len(v)),v,s=19,color=c,edgecolor="white",linewidth=.3,zorder=3)
    ymax=max(np.array(means)+np.array(sems)); step=ymax*.14
    if eight:
        specs=[(0,2,"siNC: LPS vs Control"),(2,3,"siNC: LPS+TPL vs LPS"),(2,6,"LPS: siS100A8 vs siNC")]
    else: specs=[(0,1,"siNC: LPS vs Control"),(1,3,"LPS: siS100A8 vs siNC")]
    for j,(a,b,name) in enumerate(specs): bracket(ax,a,b,ymax+step*(j+.15),step*.11,ptxt(pval(panel,name)))
    treatment_short={"Control":"C","TPL":"T","LPS":"L","LPS+TPL":"L+T"}
    ax.set_ylim(0,ymax+step*(len(specs)+.75)); ax.set_xticks(range(len(order)),[treatment_short[b] for a,b in order],rotation=0,ha="center",fontsize=6.8); ax.set_ylabel(title); ax.grid(axis="y",color="#E7EBEE",lw=.55); ax.set_axisbelow(True)
    if eight:
        ax.text(1.5,-.22,"siNC",transform=ax.get_xaxis_transform(),ha="center",va="top",clip_on=False)
        ax.text(5.5,-.22,"siS100A8",transform=ax.get_xaxis_transform(),ha="center",va="top",clip_on=False)
        ax.axvline(3.5,color="#D9DEE1",lw=.7)
    else:
        ax.text(.5,-.22,"siNC",transform=ax.get_xaxis_transform(),ha="center",va="top",clip_on=False)
        ax.text(2.5,-.22,"siS100A8",transform=ax.get_xaxis_transform(),ha="center",va="top",clip_on=False)
        ax.axvline(1.5,color="#D9DEE1",lw=.7)
    fig.subplots_adjust(left=.18,right=.98,bottom=.31,top=.97); save(fig,stem)
def interaction_plot():
    z=pd.read_csv(ROOT/"文章"/"返修投稿新主图"/"interaction_effects.csv")
    z["endpoint_short"]=np.where(z.endpoint.str.contains("p-NF"),"p-NF-κB/NF-κB","S100A8/β-actin")
    rows=["p-NF-κB/NF-κB","S100A8/β-actin"]
    cols=["siNC","siS100A8","Dependency contrast"]
    est=np.full((2,3),np.nan); lo=np.full_like(est,np.nan); hi=np.full_like(est,np.nan)
    for _,r in z.iterrows():
        i=rows.index(r.endpoint_short); j=cols.index(r.siRNA)
        est[i,j]=r.estimate; lo[i,j]=r.lower; hi[i,j]=r.upper
    from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list("interaction",["#4C78A8","#F2F3F2","#C97832"])
    lim=max(.45,float(np.nanmax(np.abs(est))))
    fig,ax=plt.subplots(figsize=(7.0,2.8)); im=ax.imshow(est,aspect="auto",cmap=cmap,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    ax.set_xticks(range(3),["siNC","siS100A8","Dependency\ncontrast"],fontsize=7.5)
    ax.set_yticks(range(2),rows,fontsize=7.5)
    for i in range(2):
        for j in range(3):
            ax.text(j,i,f"{est[i,j]:+.2f}\n[{lo[i,j]:+.2f}, {hi[i,j]:+.2f}]",ha="center",va="center",fontsize=6.8,
                    color="white" if abs(est[i,j])>lim*.52 else INK,fontweight="bold" if j==2 else "normal")
    ax.set_xlabel("TPL effect under LPS: estimate [95% CI]",fontsize=8)
    ax.set_xticks(np.arange(-.5,3,1),minor=True); ax.set_yticks(np.arange(-.5,2,1),minor=True)
    ax.grid(which="minor",color="white",linewidth=2); ax.tick_params(which="minor",length=0)
    for s in ax.spines.values(): s.set_visible(False)
    # The effect and interval are printed in every cell; the colorbar would be
    # redundant and too narrow after final multi-panel scaling.
    fig.subplots_adjust(left=.24,right=.93,bottom=.24,top=.96); save(fig,"Fig8_interaction_effect")
def main():
    fig7_qpcr(); fig7_wb(); kd_plot("N","Fig8_pNFkB","p-NF-κB / total NF-κB",True); kd_plot("O","Fig8_S100A8_8group","S100A8 / β-actin",True); kd_plot("Q","Fig8_CD74","CD74 / β-actin",False); kd_plot("R","Fig8_HLADRA","HLA-DRα / β-actin",False); kd_plot("S","Fig8_S100A8_4group","S100A8 / β-actin",False); interaction_plot()
if __name__=="__main__": main()
