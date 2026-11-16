#!/usr/bin/env python3
# panel_hh.py ORG CLASS
#   ORG   = rice | human
#   CLASS = multiple_v_single | present_v_absent
# ONE atomic head-to-head panel: 3 tools (CRISPOR/GuideScan2/kRISP-meR), 3-category
# stacked (accept/reject/not-generated). No title, no row tag, no legend (added in LaTeX).
# Figure-2 fonts, soothing palette. Output: hh_<org>_<class>.pdf/.png
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORG, CLS = sys.argv[1], sys.argv[2]
OUT=os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#4a4a4a",
    "axes.labelcolor":"#2b2b2b","text.color":"#2b2b2b","xtick.color":"#4a4a4a","ytick.color":"#4a4a4a",
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":18,"axes.linewidth":1.3})
ACC="#6BA292"; REJ="#E4A66B"; NOG="#D4D2CC"
tools=["CRISPOR","GuideScan2","kRISP-meR"]

DATA={
 ("human","multiple_v_single"): {"N":116,"kRISP-meR":(3,86,27),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
 ("human","present_v_absent"):  {"N":713,"kRISP-meR":(447,259,7),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
 ("rice","multiple_v_single"):  {"N":133,"kRISP-meR":(7,121,5),"CRISPOR":(64,69,0),"GuideScan2":(107,26,0)},
 ("rice","present_v_absent"):   {"N":237,"kRISP-meR":(161,67,9),"CRISPOR":(140,97,0),"GuideScan2":(0,0,237)},
}
D=DATA[(ORG,CLS)]
def pct(t):
    N=D["N"]; a,r,g=D[t]; return (100*a/N,100*r/N,100*g/N) if N else (0,0,0)

fig,ax=plt.subplots(figsize=(5.8,5.6))
x=np.arange(3)
A=[pct(t)[0] for t in tools]; R=[pct(t)[1] for t in tools]; G=[pct(t)[2] for t in tools]
ax.bar(x,A,0.66,color=ACC,edgecolor="white",lw=1.0,label="generated & accepted")
ax.bar(x,R,0.66,bottom=A,color=REJ,edgecolor="white",lw=1.0,label="generated & rejected")
ax.bar(x,G,0.66,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=1.0,label="not generated")
for i in range(3):
    for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
        if val>=9: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=16,color="#2b2b2b")
ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=16)
ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=15,width=1.3,length=5)
ax.set_ylabel("% of guides",fontsize=17)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.16,right=0.98,top=0.98,bottom=0.11)
for ext in ("pdf","png"):
    p=f"{OUT}/hh_{ORG}_{CLS}.{ext}"; fig.savefig(p,dpi=300); print(f"wrote {p}")
plt.close()
