#!/usr/bin/env python3
# panel_cfd.py ORG CLASS
#   ORG   = rice | human
#   CLASS = multiple_v_single | present_v_absent
# ONE atomic assembly-CFD panel (ACCEPTED guides): 3 tools, box + strip + mean diamond,
# reference lines at 1 (single-copy) and 2 (duplicated). No title/rowtag/legend.
# Figure-2 fonts, soothing palette. Output: cfd_<org>_<class>.pdf/.png
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORG, CLS = sys.argv[1], sys.argv[2]
DATA="/mnt/windows-data/Thesis/data"
OUT=os.path.dirname(os.path.abspath(__file__))
TSV={"rice":f"{DATA}/asm_cfd_dist_rice_n22.tsv","human":f"{DATA}/asm_cfd_dist_human_chm13.tsv"}[ORG]

plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#4a4a4a",
    "axes.labelcolor":"#2b2b2b","text.color":"#2b2b2b","xtick.color":"#4a4a4a","ytick.color":"#4a4a4a",
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":18,"axes.linewidth":1.3})
TOOLCOL={"CRISPOR":"#E4A66B","GuideScan2":"#C9A7C2","kRISP-meR":"#6BA292"}
tools=["CRISPOR","GuideScan2","kRISP-meR"]
CAP=8.0

df=pd.read_csv(TSV,sep="\t")
df=df[(df["verdict"]=="accept") & (df["class"]==CLS)]

fig,ax=plt.subplots(figsize=(5.8,5.6))
for xi,t in enumerate(tools):
    vals=df[df["tool"]==t]["cfd_cuts"].values
    if len(vals)==0:
        ax.text(xi,0.6,"designs\nnone",ha="center",va="center",fontsize=13,color="#999",style="italic")
        continue
    capped=np.clip(vals,None,CAP); over=int((vals>CAP).sum())
    ax.boxplot([capped],positions=[xi],widths=0.55,patch_artist=True,showfliers=False,
               medianprops=dict(color="#2b2b2b",lw=2),
               whiskerprops=dict(color="#4a4a4a",lw=1.3),capprops=dict(color="#4a4a4a",lw=1.3),
               boxprops=dict(facecolor=TOOLCOL[t],edgecolor="#4a4a4a",lw=1.2,alpha=0.75))
    jit=(np.random.RandomState(xi).rand(len(capped))-0.5)*0.28
    ax.scatter(np.full(len(capped),xi)+jit,capped,s=11,color=TOOLCOL[t],alpha=0.45,edgecolors="none",zorder=3)
    ax.scatter([xi],[np.mean(capped)],marker="D",s=65,facecolor="white",edgecolor=TOOLCOL[t],lw=1.8,zorder=5)
    if over: ax.text(xi,CAP+0.15,f"+{over}",ha="center",va="bottom",fontsize=13,color=TOOLCOL[t])
ax.axhline(1.0,ls=":",lw=1.5,color="#3a8f7a",zorder=1)   # single-copy safe
ax.axhline(2.0,ls="--",lw=1.5,color="#888",zorder=1)     # duplicated
ax.set_xlim(-0.6,len(tools)-0.4)
ax.set_ylim(0,CAP+0.6); ax.set_yticks([0,2,4,6,8])
ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools,fontsize=16)
ax.set_ylabel("CFD expected cuts",fontsize=16)
ax.tick_params(labelsize=15,width=1.3,length=5)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.16,right=0.98,top=0.98,bottom=0.11)
for ext in ("pdf","png"):
    p=f"{OUT}/cfd_{ORG}_{CLS}.{ext}"; fig.savefig(p,dpi=300); print(f"wrote {p}")
plt.close()
