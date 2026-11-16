#!/usr/bin/env python3
# fig_cfd_accepted.py ORG   (rice | human)
# CFD panel for the discordance figure's bottom block: assembly-CFD expected cuts for
# ACCEPTED guides, per tool, for the two headline classes. One organism per call.
# Mirrors the head-to-head block: 2 rows (classes) x 3 tools (CRISPOR/GuideScan2/kRISP-meR),
# box + strip. Figure-2 fonts, soothing palette. Vector PDF + PNG.
#
# Data: data/asm_cfd_dist_{rice_n22,human_chm13}.tsv  (cols: cfd_cuts class group guide tool truth verdict)
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORG = sys.argv[1] if len(sys.argv)>1 else "rice"
ROWTAG = (len(sys.argv)<3) or (sys.argv[2] not in ("0","no","false"))
DATA="/mnt/windows-data/Thesis/data"
OUT=os.path.dirname(os.path.abspath(__file__))
TSV = {"rice":f"{DATA}/asm_cfd_dist_rice_n22.tsv","human":f"{DATA}/asm_cfd_dist_human_chm13.tsv"}[ORG]

plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#4a4a4a",
    "axes.labelcolor":"#2b2b2b","text.color":"#2b2b2b","xtick.color":"#4a4a4a","ytick.color":"#4a4a4a",
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":18,"axes.linewidth":1.3})
# soothing palette, per tool (consistent with head-to-head accent green for kRISP-meR)
TOOLCOL={"CRISPOR":"#E4A66B","GuideScan2":"#C9A7C2","kRISP-meR":"#6BA292"}
tools=["CRISPOR","GuideScan2","kRISP-meR"]
CASES=["multiple_v_single","present_v_absent"]
ROW_TAG={"multiple_v_single":"Duplicated in individual","present_v_absent":"Individual-specific"}
CAP=8.0   # y-axis cap; guides above are shown as an overflow marker

df=pd.read_csv(TSV,sep="\t")
df=df[df["verdict"]=="accept"]

fig,axes=plt.subplots(2,1,figsize=(6.6,8.4),sharex=True)
for ri,case in enumerate(CASES):
    ax=axes[ri]; sub=df[df["class"]==case]
    for xi,t in enumerate(tools):
        vals=sub[sub["tool"]==t]["cfd_cuts"].values
        if len(vals)==0:
            ax.text(xi,0.5,"designs\nnone",ha="center",va="center",fontsize=12,color="#999",style="italic")
            continue
        capped=np.clip(vals,None,CAP); over=int((vals>CAP).sum())
        # box
        bp=ax.boxplot([capped],positions=[xi],widths=0.55,patch_artist=True,showfliers=False,
                      medianprops=dict(color="#2b2b2b",lw=2),
                      whiskerprops=dict(color="#4a4a4a",lw=1.3),capprops=dict(color="#4a4a4a",lw=1.3),
                      boxprops=dict(facecolor=TOOLCOL[t],edgecolor="#4a4a4a",lw=1.2,alpha=0.75))
        # strip
        jit=(np.random.RandomState(xi+ri*10).rand(len(capped))-0.5)*0.28
        ax.scatter(np.full(len(capped),xi)+jit,capped,s=10,color=TOOLCOL[t],alpha=0.45,edgecolors="none",zorder=3)
        # mean diamond
        ax.scatter([xi],[np.mean(capped)],marker="D",s=60,facecolor="white",edgecolor=TOOLCOL[t],lw=1.8,zorder=5)
        if over: ax.text(xi,CAP+0.15,f"+{over}",ha="center",va="bottom",fontsize=12,color=TOOLCOL[t])
    ax.axhline(1.0,ls=":",lw=1.4,color="#3a8f7a",zorder=1)   # single-copy safe
    ax.axhline(2.0,ls="--",lw=1.4,color="#888",zorder=1)     # duplicated
    ax.set_ylim(0,CAP+0.6); ax.set_yticks([0,2,4,6,8])
    ax.set_ylabel("Assembly-CFD\nexpected cuts",fontsize=15)
    ax.tick_params(labelsize=14,width=1.3,length=5)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    (ROWTAG) and ax.annotate(ROW_TAG[case],xy=(-0.30,0.5),xycoords="axes fraction",rotation=90,
                ha="center",va="center",fontsize=16,color="#3a3a3a",weight="bold")
axes[1].set_xticks(range(len(tools))); axes[1].set_xticklabels(tools,fontsize=15)
axes[0].set_xlim(-0.6,len(tools)-0.4)
fig.subplots_adjust(left=0.20,right=0.96,top=0.96,bottom=0.08,hspace=0.16)
for ext in ("pdf","png"):
    p=f"{OUT}/fig_cfd_{ORG}.{ext}"; fig.savefig(p,dpi=300,bbox_inches="tight"); print(f"wrote {p}")
plt.close()
