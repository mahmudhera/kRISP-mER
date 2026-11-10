#!/usr/bin/env python3
# fig_hh_panel.py ORG [rowtag]
#   ORG    = rice | human
#   rowtag = 1 (show left-side class row tags) | 0 (omit; for the right column)
# Per-organism head-to-head panel: 2 headline classes x 3 tools (CRISPOR/GuideScan2/kRISP-meR),
# 3-category stacked (accept/reject/not-generated). Figure-2 fonts, soothing palette.
# No organism title (added in LaTeX). Row tags only on the left (rice) panel.
# Vector PDF + PNG.
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORG    = sys.argv[1] if len(sys.argv)>1 else "rice"
ROWTAG = (len(sys.argv)<3) or (sys.argv[2] not in ("0","no","false"))
OUT=os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#4a4a4a",
    "axes.labelcolor":"#2b2b2b","text.color":"#2b2b2b","xtick.color":"#4a4a4a","ytick.color":"#4a4a4a",
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":18,"axes.linewidth":1.3})
ACC="#6BA292"; REJ="#E4A66B"; NOG="#D4D2CC"; LABEL="#3a3a3a"
tools=["CRISPOR","GuideScan2","kRISP-meR"]
CASES=["multiple_v_single","present_v_absent"]
ROW_TAG={"multiple_v_single":"Duplicated in individual","present_v_absent":"Individual-specific"}

# (accept, reject, not_generated), N  -- from headtohead_classN_*_bothstrand.tsv
DATA={
 ("human","multiple_v_single"): {"N":116,"kRISP-meR":(3,86,27),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
 ("human","present_v_absent"):  {"N":713,"kRISP-meR":(447,259,7),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
 ("rice","multiple_v_single"):  {"N":133,"kRISP-meR":(7,121,5),"CRISPOR":(64,69,0),"GuideScan2":(107,26,0)},
 ("rice","present_v_absent"):   {"N":237,"kRISP-meR":(161,67,9),"CRISPOR":(140,97,0),"GuideScan2":(0,0,237)},
}
def pct(D,t):
    N=D["N"]; a,r,g=D[t]; return (100*a/N,100*r/N,100*g/N) if N else (0,0,0)

fig,axes=plt.subplots(2,1,figsize=(6.6,8.4),sharex=True)
for ri,case in enumerate(CASES):
    ax=axes[ri]; D=DATA[(ORG,case)]; x=np.arange(3)
    A=[pct(D,t)[0] for t in tools]; R=[pct(D,t)[1] for t in tools]; G=[pct(D,t)[2] for t in tools]
    ax.bar(x,A,0.64,color=ACC,edgecolor="white",lw=1.0,label="generated & accepted")
    ax.bar(x,R,0.64,bottom=A,color=REJ,edgecolor="white",lw=1.0,label="generated & rejected")
    ax.bar(x,G,0.64,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=1.0,label="not generated")
    for i in range(3):
        for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
            if val>=9: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=15,color="#2b2b2b")
    ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=14,width=1.3,length=5)
    ax.set_ylabel("% of class guides",fontsize=15)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    if ROWTAG:
        ax.annotate(ROW_TAG[case],xy=(-0.30,0.5),xycoords="axes fraction",rotation=90,
                    ha="center",va="center",fontsize=16,color=LABEL,weight="bold")
axes[1].set_xticks(range(len(tools))); axes[1].set_xticklabels(tools,fontsize=15)
fig.subplots_adjust(left=0.20,right=0.96,top=0.96,bottom=0.08,hspace=0.16)
for ext in ("pdf","png"):
    p=f"{OUT}/fig_hh_{ORG}.{ext}"; fig.savefig(p,dpi=300,bbox_inches="tight"); print(f"wrote {p}")
plt.close()
