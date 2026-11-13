#!/usr/bin/env python3
# plot_anchor_comparison_rice.py — RICE: individual-anchored (N22) vs
# reference-anchored (japonica) target windows, side by side. Same tools, same
# classes, both-strand kRISP-meR; only the target SOURCE differs.
# 3-category stacked bars: generated&accepted / generated&rejected / not-generated.
# Reads the two local TSVs; writes plots/reference-based-targets/.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

OUT="/mnt/windows-data/Thesis/plots/reference-based-targets"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#86B3C4"; REJ="#E3A15B"; NOG="#CDD2DA"
tools=["kRISP-meR","CRISPOR","GuideScan2"]

def load(tsv):
    C={}
    with open(tsv) as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            cls=row["class"]; C.setdefault(cls,{"N":int(row["N"])})
            C[cls][row["tool"]]=(int(row["accept"]),int(row["reject"]),int(row["not_designed"]))
    return C
IND =load("/mnt/windows-data/Thesis/data/headtohead_classN_n22_bothstrand.tsv")  # N22-anchored
REF =load("/mnt/windows-data/Thesis/data/headtohead_classN_n22_REF.tsv")         # japonica-anchored

TITLE={"multiple_v_single":"Multiple in individual / Single in reference\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in reference\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in reference\n(should accept)",
"present_v_absent":"Present in individual / Absent in reference\n(should accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
# which anchor gives EXACT windows for each class (for the annotation)
EXACT={"multiple_v_single":"both","present_v_absent":"IND",
       "ref_present_v_absent":"REF","ref_multiple_v_single":"REF"}

def pct(D,cls,t):
    N=D[cls]["N"]; a,r,g=D[cls][t]; return (a/N*100,r/N*100,g/N*100) if N else (0,0,0)

# 2 rows (IND top, REF bottom) x 4 classes
fig,axes=plt.subplots(2,4,figsize=(15.5,7.6))
POS={("IND",c):(0,i) for i,c in enumerate(ORDER)}
POS.update({("REF",c):(1,i) for i,c in enumerate(ORDER)})
for tag,D in [("IND",IND),("REF",REF)]:
    for cls in ORDER:
        r,c=POS[(tag,cls)]; ax=axes[r,c]
        x=np.arange(3)
        A=[pct(D,cls,t)[0] for t in tools]; R=[pct(D,cls,t)[1] for t in tools]; G=[pct(D,cls,t)[2] for t in tools]
        ax.bar(x,A,0.62,color=ACC,edgecolor="white",lw=0.6,label="generated & accepted")
        ax.bar(x,R,0.62,bottom=A,color=REJ,edgecolor="white",lw=0.6,label="generated & rejected")
        ax.bar(x,G,0.62,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=0.6,label="not generated")
        for i in range(3):
            for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
                if val>=8: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=7,color="#333")
        ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=7.5)
        ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=7.5)
        if r==0: ax.set_title(TITLE[cls],fontsize=8,color="#2b2f3a")
        for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")

fig.text(0.008,0.72,"Targets from\nN22 (individual)",ha="left",va="center",fontsize=10.5,color="#2b2f3a",weight="bold")
fig.text(0.008,0.30,"Targets from\njaponica (reference)",ha="left",va="center",fontsize=10.5,color="#2b2f3a",weight="bold")
fig.suptitle("Rice: target windows anchored to the individual (N22) vs the reference (japonica)",
             fontsize=12.5,color="#2b2f3a",weight="bold",y=0.995)

h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="lower center",ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,-0.005))
fig.tight_layout(rect=[0.11,0.03,1,0.94])
fig.savefig(f"{OUT}/fig_anchor_comparison_rice.png",dpi=200,bbox_inches="tight"); plt.close()
print(f"wrote {OUT}/fig_anchor_comparison_rice.png")
