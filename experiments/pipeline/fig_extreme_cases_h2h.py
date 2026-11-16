#!/usr/bin/env python3
# fig_extreme_cases_h2h.py
# Manuscript Results subsection figure: the TWO EXTREME discordance cases where the
# reference misleads, kRISP-meR vs CRISPOR vs GuideScan2, both organisms.
#   multiple_v_single  (duplicated in individual, single in reference; should REJECT)
#   present_v_absent   (present in individual, absent from reference; should ACCEPT)
# 3-category stacked: generated&accepted / generated&rejected / not-generated.
# Layout: 2 rows (cases) x 2 cols (human | rice).
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os

OUT="/mnt/windows-data/Thesis/plots"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#86B3C4"; REJ="#E3A15B"; NOG="#CDD2DA"
tools=["kRISP-meR","CRISPOR","GuideScan2"]

# (accept, reject, not_designed), N  -- from headtohead_classN_*_bothstrand.tsv
DATA={
 ("human","multiple_v_single"): {"N":116,"kRISP-meR":(3,86,27),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
 ("human","present_v_absent"):  {"N":713,"kRISP-meR":(447,259,7),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
 ("rice","multiple_v_single"):  {"N":133,"kRISP-meR":(7,121,5),"CRISPOR":(64,69,0),"GuideScan2":(107,26,0)},
 ("rice","present_v_absent"):   {"N":237,"kRISP-meR":(161,67,9),"CRISPOR":(140,97,0),"GuideScan2":(0,0,237)},
}
CASE_TITLE={
 "multiple_v_single":"Duplicated in individual, single in reference\n(guide cuts multiple times — should REJECT)",
 "present_v_absent":"Present in individual, absent from reference\n(individual-only target — should ACCEPT)"}
ORG_TITLE={"human":"Human  (CHM13 vs GRCh38)","rice":"Rice  (N22 vs japonica)"}
CASES=["multiple_v_single","present_v_absent"]; ORGS=["human","rice"]

def pct(D,cls_key,t):
    N=D["N"]; a,r,g=D[t]; return (100*a/N,100*r/N,100*g/N) if N else (0,0,0)

fig,axes=plt.subplots(2,2,figsize=(10.4,8.6))
for ri,case in enumerate(CASES):
    for ci,org in enumerate(ORGS):
        ax=axes[ri,ci]; D=DATA[(org,case)]
        x=np.arange(3)
        A=[pct(D,case,t)[0] for t in tools]; R=[pct(D,case,t)[1] for t in tools]; G=[pct(D,case,t)[2] for t in tools]
        ax.bar(x,A,0.62,color=ACC,edgecolor="white",lw=0.6,label="generated & accepted")
        ax.bar(x,R,0.62,bottom=A,color=REJ,edgecolor="white",lw=0.6,label="generated & rejected")
        ax.bar(x,G,0.62,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=0.6,label="not generated")
        for i in range(3):
            for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
                if val>=7: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=8,color="#333")
        ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=8.5)
        ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=8)
        for s in ("top","right"): ax.spines[s].set_visible(False)
        if ri==0: ax.set_title(ORG_TITLE[org],fontsize=11,color="#2b2f3a",weight="bold",pad=8)
    axes[ri,0].set_ylabel("% of class guides",fontsize=10)
    # case label spanning the row (left of the two panels)
    axes[ri,0].annotate(CASE_TITLE[case],xy=(-0.30,0.5),xycoords="axes fraction",
                        rotation=90,ha="center",va="center",fontsize=9.2,color="#2b2f3a",weight="bold")

h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="lower center",ncol=3,frameon=False,fontsize=9.5,bbox_to_anchor=(0.5,-0.005))
fig.subplots_adjust(left=0.14,right=0.97,top=0.92,bottom=0.09,hspace=0.28,wspace=0.20)
png=f"{OUT}/fig_extreme_cases_h2h.png"
fig.savefig(png,dpi=200,bbox_inches="tight"); plt.close()
print(f"wrote {png}")
