#!/usr/bin/env python3
# 10_plot_headtohead_combined.py — two-panel human|rice head-to-head figure,
# THREE-category view (generated&accepted / generated&rejected / not-generated).
# Left block = human (CHM13 vs GRCh38), right block = rice (N22 vs japonica).
# Counts baked in from krispmer_class_tables.txt (Tables 5 human / 5b rice).
# Writes /mnt/windows-data/Thesis/plots/fig_h2h_combined.png
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

# class -> {N, tool:(gen&accept, gen&reject, not-generated)}   [3-category]
# CLASS-GUIDE denominator (N = all class guides present in the 50 windows).
#
# SOURCE OF TRUTH = the two VM TSVs (headtohead_classN_human.tsv /
# headtohead_classN_n22.tsv). Read if found; else the baked-in fallbacks below,
# which MUST match the TSVs. Do NOT hand-edit the fallbacks without re-deriving.
import csv
def _load(cands, fb):
    for tsv in cands:
        if os.path.exists(tsv):
            C={}
            with open(tsv) as fh:
                for row in csv.DictReader(fh,delimiter="\t"):
                    cls=row["class"]; C.setdefault(cls,{"N":int(row["N"])})
                    C[cls][row["tool"]]=(int(row["accept"]),int(row["reject"]),int(row["not_designed"]))
            print(f"[source] {tsv}"); return C
    print("[source] fallback (verify it matches the VM TSV)"); return fb

_HUMAN_FB={
"multiple_v_single":    {"N":116,"kRISP-meR":(2,48,66),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
"ref_present_v_absent": {"N":191,"kRISP-meR":(0,10,181),"CRISPOR":(40,151,0),"GuideScan2":(120,71,0)},
"ref_multiple_v_single":{"N":160,"kRISP-meR":(18,51,91),"CRISPOR":(1,159,0),"GuideScan2":(118,42,0)},
"present_v_absent":     {"N":713,"kRISP-meR":(199,117,397),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
}
_RICE_FB={
"multiple_v_single":    {"N":133,"kRISP-meR":(4,61,68),"CRISPOR":(64,69,0),"GuideScan2":(107,26,0)},
"ref_present_v_absent": {"N":73, "kRISP-meR":(2,4,67), "CRISPOR":(35,38,0),"GuideScan2":(56,17,0)},
"ref_multiple_v_single":{"N":109,"kRISP-meR":(27,37,45),"CRISPOR":(10,99,0),"GuideScan2":(85,24,0)},
"present_v_absent":     {"N":237,"kRISP-meR":(75,33,129),"CRISPOR":(140,97,0),"GuideScan2":(0,0,237)},
}
HUMAN=_load([os.path.expanduser("~/human_assemblies_kmer_count/july-10/headtohead_classN_human.tsv"),
             os.path.join(os.path.dirname(__file__),"headtohead_classN_human.tsv"),
             "/mnt/windows-data/Thesis/data/headtohead_classN_human.tsv"], _HUMAN_FB)
RICE=_load([os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh/headtohead_classN_n22.tsv"),
            os.path.join(os.path.dirname(__file__),"headtohead_classN_n22.tsv"),
            "/mnt/windows-data/Thesis/data/headtohead_classN_n22.tsv"], _RICE_FB)
TITLE={"multiple_v_single":"Multiple in individual / Single in reference\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in reference\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in reference\n(should accept)",
"present_v_absent":"Present in individual / Absent in reference\n(should accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]

def pct(D,cls,t):
    N=D[cls]["N"]; a,r,g=D[cls][t]; return (a/N*100,r/N*100,g/N*100) if N else (0,0,0)

fig,axes=plt.subplots(2,4,figsize=(15.5,7.4))
POS={ ("H","multiple_v_single"):(0,0),("H","ref_present_v_absent"):(0,1),
      ("H","ref_multiple_v_single"):(1,0),("H","present_v_absent"):(1,1),
      ("R","multiple_v_single"):(0,2),("R","ref_present_v_absent"):(0,3),
      ("R","ref_multiple_v_single"):(1,2),("R","present_v_absent"):(1,3) }
for tag,D in [("H",HUMAN),("R",RICE)]:
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
        ax.set_title(TITLE[cls],fontsize=8,color="#2b2f3a")
        for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides"); axes[1,0].set_ylabel("% of class guides")

fig.text(0.28,0.965,"Human  (CHM13 vs GRCh38)",ha="center",fontsize=12,color="#2b2f3a",weight="bold")
fig.text(0.75,0.965,"Rice  (N22 vs japonica)",ha="center",fontsize=12,color="#2b2f3a",weight="bold")
fig.add_artist(plt.Line2D([0.515,0.515],[0.05,0.93],color="#cccccc",lw=1.0))

h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="lower center",ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,-0.005))
fig.tight_layout(rect=[0,0.03,1,0.95])
fig.savefig(f"{OUT}/fig_h2h_combined.png",dpi=200,bbox_inches="tight"); plt.close()
print(f"wrote {OUT}/fig_h2h_combined.png")
