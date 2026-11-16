#!/usr/bin/env python3
# 09_plot_headtohead_n22.py — rice N22 head-to-head figures from the fresh TSV.
# Reads fresh/headtohead_common_n22.tsv (no hardcoded counts) -> two PNGs:
#   fig4b_3category_n22.png  (accept / reject / not-generated, stacked)
#   fig5b_2category_n22.png  (recommended vs not-recommended)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

FRESH = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
TSV   = f"{FRESH}/headtohead_common_n22.tsv"
OUT   = f"{FRESH}/plots"; os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#86B3C4"; REJ="#E3A15B"; NOG="#CDD2DA"; DONT="#AAB1BD"
tools=["kRISP-meR","CRISPOR","GuideScan2"]

# ---- load counts from TSV: {class:{N, tool:(accept,reject,not_designed)}} ----
C={}
with open(TSV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        cls=row["class"]; C.setdefault(cls,{"N":int(row["N"])})
        C[cls][row["tool"]]=(int(row["accept"]),int(row["reject"]),int(row["not_designed"]))

TITLE={"multiple_v_single":"Multiple in individual / Single in reference\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in reference\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in reference\n(should accept)",
"present_v_absent":"Present in individual / Absent in reference\n(should accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]

def pct(cls,t):
    N=C[cls]["N"]; a,r,g=C[cls][t]; return (a/N*100,r/N*100,g/N*100) if N else (0,0,0)

# ---------- FIG 4b: 3 categories ----------
fig,axes=plt.subplots(2,2,figsize=(9.2,7.2))
for ax,cls in zip(axes.flat,ORDER):
    x=np.arange(3)
    A=[pct(cls,t)[0] for t in tools]; R=[pct(cls,t)[1] for t in tools]; G=[pct(cls,t)[2] for t in tools]
    ax.bar(x,A,0.62,color=ACC,edgecolor="white",lw=0.6,label="generated & accepted")
    ax.bar(x,R,0.62,bottom=A,color=REJ,edgecolor="white",lw=0.6,label="generated & rejected")
    ax.bar(x,G,0.62,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=0.6,label="not generated")
    for i in range(3):
        for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
            if val>=7: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=7.5,color="#333")
    ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=8)
    ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=8)
    ax.set_title(f"{TITLE[cls]}  (N={C[cls]['N']})",fontsize=8.5,color="#2b2f3a")
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")
h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="upper center",ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
fig.suptitle("Rice N22 vs japonica reference — 3-tool head-to-head",fontsize=11,y=1.0,x=0.5)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig(f"{OUT}/fig4b_3category_n22.png",dpi=220); plt.close()

# ---------- FIG 5b: 2 categories ----------
fig,axes=plt.subplots(2,2,figsize=(9.2,7.2))
for ax,cls in zip(axes.flat,ORDER):
    x=np.arange(3)
    REC=[pct(cls,t)[0] for t in tools]; DON=[pct(cls,t)[1]+pct(cls,t)[2] for t in tools]
    ax.bar(x,REC,0.62,color=ACC,edgecolor="white",lw=0.6,label="recommended (accept)")
    ax.bar(x,DON,0.62,bottom=REC,color=DONT,edgecolor="white",lw=0.6,label="not recommended (reject or not generated)")
    for i in range(3):
        if REC[i]>=6: ax.text(i,REC[i]/2,f"{int(round(REC[i]))}",ha="center",va="center",fontsize=8,color="white" if REC[i]>18 else "#333")
    ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=8)
    ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=8)
    ax.set_title(f"{TITLE[cls]}  (N={C[cls]['N']})",fontsize=8.5,color="#2b2f3a")
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")
h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="upper center",ncol=2,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
fig.suptitle("Rice N22 vs japonica reference — recommended vs not",fontsize=11,y=1.0,x=0.5)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig(f"{OUT}/fig5b_2category_n22.png",dpi=220); plt.close()
print(f"wrote {OUT}/fig4b_3category_n22.png and fig5b_2category_n22.png")
