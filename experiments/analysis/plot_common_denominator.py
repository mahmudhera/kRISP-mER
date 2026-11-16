#!/usr/bin/env python3
# plot_common_denominator.py — HUMAN head-to-head figures (fig4/fig5).
# CLASS-GUIDE denominator (N = all class-defining guides in the 50 windows).
#
# SOURCE OF TRUTH = the VM TSV headtohead_classN_human.tsv (from
# 08_headtohead_classN_human.py). Reads that TSV if found; otherwise falls back
# to the baked-in counts below, which MUST match the TSV. Do NOT hand-edit the
# fallback without re-deriving it from the TSV.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

OUT="/mnt/windows-data/Thesis/plots"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#86B3C4"; REJ="#E3A15B"; NOG="#CDD2DA"; DONT="#AAB1BD"
tools=["kRISP-meR","CRISPOR","GuideScan2"]

TSV_CANDIDATES=[
    os.path.expanduser("~/human_assemblies_kmer_count/july-10/headtohead_classN_human.tsv"),
    os.path.join(os.path.dirname(__file__),"headtohead_classN_human.tsv"),
    "/mnt/windows-data/Thesis/data/headtohead_classN_human.tsv",
]
# fallback (MUST equal the TSV): class -> {N, tool:(accept,reject,not_generated)}
FALLBACK={
"multiple_v_single":{"N":116,"kRISP-meR":(2,48,66),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
"ref_present_v_absent":{"N":191,"kRISP-meR":(0,10,181),"CRISPOR":(40,151,0),"GuideScan2":(120,71,0)},
"ref_multiple_v_single":{"N":160,"kRISP-meR":(18,51,91),"CRISPOR":(1,159,0),"GuideScan2":(118,42,0)},
"present_v_absent":{"N":713,"kRISP-meR":(199,117,397),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
}
def load_counts():
    for tsv in TSV_CANDIDATES:
        if os.path.exists(tsv):
            C={}
            with open(tsv) as fh:
                for row in csv.DictReader(fh,delimiter="\t"):
                    cls=row["class"]; C.setdefault(cls,{"N":int(row["N"])})
                    C[cls][row["tool"]]=(int(row["accept"]),int(row["reject"]),int(row["not_designed"]))
            print(f"[source] read counts from {tsv}"); return C
    print("[source] TSV not found; using baked-in fallback (verify it matches the VM TSV)")
    return FALLBACK
C=load_counts()

TITLE={"multiple_v_single":"Multiple in individual / Single in reference\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in reference\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in reference\n(should accept)",
"present_v_absent":"Present in individual / Absent in reference\n(should accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]

def pct(cls,t):
    N=C[cls]["N"]; a,r,g=C[cls][t]; return a/N*100,r/N*100,g/N*100

# ---------- FIG 4: 3 categories ----------
fig,axes=plt.subplots(2,2,figsize=(9.2,7.2))
for ax,cls in zip(axes.flat,ORDER):
    x=np.arange(3)
    A=[pct(cls,t)[0] for t in tools]; R=[pct(cls,t)[1] for t in tools]; G=[pct(cls,t)[2] for t in tools]
    ax.bar(x,A,0.62,color=ACC,edgecolor="white",lw=0.6,label="generated & accepted")
    ax.bar(x,R,0.62,bottom=A,color=REJ,edgecolor="white",lw=0.6,label="generated & rejected")
    ax.bar(x,G,0.62,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=0.6,label="not generated")
    for i in range(3):
        for val,bot,dark in [(A[i],0,False),(R[i],A[i],False),(G[i],A[i]+R[i],True)]:
            if val>=7: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=7.5,color="#333")
    ax.set_xticks(x); ax.set_xticklabels(tools,fontsize=8)
    ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=8)
    ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")
h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="upper center",ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=[0,0,1,0.955]); fig.savefig(f"{OUT}/fig4_3category.png",dpi=220); plt.close()

# ---------- FIG 5: 2 categories (reject + not-generated merged) ----------
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
    ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")
h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="upper center",ncol=2,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=[0,0,1,0.955]); fig.savefig(f"{OUT}/fig5_2category.png",dpi=220); plt.close()
print("done")
