#!/usr/bin/env python3
# AA2_plot_headtohead_human_with_aa.py
# HUMAN head-to-head, 4 classes, with AlleleAnalyzer LOCKED IN as a 4th bar.
#
# The 3 specificity tools (kRISP-meR / CRISPOR / GuideScan2) are scored on the
# CLASS-GUIDE denominator: % of class guides accepted / rejected / not-generated.
# AlleleAnalyzer is NOT a per-guide specificity tool -- it emits allele-
# discriminating guide PAIRS per het variant. Its honest head-to-head axis is
# LOCUS DISCRIMINABILITY: % of the class's target loci for which AA can build at
# least one allele-specific guide. So AA's bar is a SINGLE segment on a DIFFERENT
# denominator (loci, not class guides); it is drawn set apart and labelled, never
# stacked into the accept/reject semantics of the other three.
#
# Outputs into  <repo>/allele-results/ :
#   fig_h2h_human_with_alleleanalyzer.png
#   headtohead_human_with_aa.tsv   (the 3-tool + AA numbers actually plotted)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

OUT = "/mnt/windows-data/Thesis/allele-results"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#86B3C4"; REJ="#E3A15B"; NOG="#CDD2DA"; AAC="#8E7CC3"   # AA = distinct purple

tools=["kRISP-meR","CRISPOR","GuideScan2"]

# --- 3-tool head-to-head (BOTH-STRAND), CLASS-GUIDE denominator ---------------
# (a, r, not_generated) counts per class; matches data/headtohead_classN_human_bothstrand.tsv
HUMAN={
"multiple_v_single":    {"N":116,"kRISP-meR":(3,86,27),"CRISPOR":(59,57,0),"GuideScan2":(75,41,0)},
"ref_present_v_absent": {"N":191,"kRISP-meR":(1,16,174),"CRISPOR":(40,151,0),"GuideScan2":(120,71,0)},
"ref_multiple_v_single":{"N":160,"kRISP-meR":(45,101,14),"CRISPOR":(1,159,0),"GuideScan2":(118,42,0)},
"present_v_absent":     {"N":713,"kRISP-meR":(447,259,7),"CRISPOR":(561,152,0),"GuideScan2":(0,0,713)},
}
# --- AlleleAnalyzer discriminability: (n_window_loci, n_covered) ---------------
# from alleleanalyzer_discriminability_human.tsv (VM), sanity-checked cov==AA_ENST.
AA={
"multiple_v_single":    (50,24),   # 48%
"ref_present_v_absent": (50,14),   # 28%
"ref_multiple_v_single":(50,48),   # 96%
"present_v_absent":     (50,1),    # 2%
}
TRUTH={"multiple_v_single":"reject","ref_present_v_absent":"reject",
       "ref_multiple_v_single":"accept","present_v_absent":"accept"}
TITLE={"multiple_v_single":"Multiple in individual / Single in reference\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in reference\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in reference\n(should accept)",
"present_v_absent":"Present in individual / Absent in reference\n(should accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]

def pct3(cls,t):
    N=HUMAN[cls]["N"]; a,r,g=HUMAN[cls][t]; return (a/N*100,r/N*100,g/N*100) if N else (0,0,0)
def pctAA(cls):
    n,c=AA[cls]; return 100.0*c/n if n else 0.0

fig,axes=plt.subplots(2,2,figsize=(10.2,8.4))
POS={"multiple_v_single":(0,0),"ref_present_v_absent":(0,1),
     "ref_multiple_v_single":(1,0),"present_v_absent":(1,1)}
xlabels=tools+["AlleleAnalyzer"]
for cls in ORDER:
    r,c=POS[cls]; ax=axes[r,c]
    x=np.arange(4)
    A=[pct3(cls,t)[0] for t in tools]; R=[pct3(cls,t)[1] for t in tools]; G=[pct3(cls,t)[2] for t in tools]
    ax.bar(x[:3],A,0.62,color=ACC,edgecolor="white",lw=0.6,label="generated & accepted")
    ax.bar(x[:3],R,0.62,bottom=A,color=REJ,edgecolor="white",lw=0.6,label="generated & rejected")
    ax.bar(x[:3],G,0.62,bottom=[A[i]+R[i] for i in range(3)],color=NOG,edgecolor="white",lw=0.6,label="not generated")
    for i in range(3):
        for val,bot in [(A[i],0),(R[i],A[i]),(G[i],A[i]+R[i])]:
            if val>=8: ax.text(i,bot+val/2,f"{int(round(val))}",ha="center",va="center",fontsize=7,color="#333")
    # AA single-segment bar (different denominator: % of loci discriminable)
    aav=pctAA(cls)
    ax.bar([3],[aav],0.62,color=AAC,edgecolor="white",lw=0.6,label="AA: locus discriminable")
    ax.text(3,aav+2 if aav<92 else aav-6,f"{int(round(aav))}%",ha="center",
            va="bottom" if aav<92 else "top",fontsize=7.5,color="#4b3b7a",weight="bold")
    # separator between the 3 specificity tools and AA (different axis semantics)
    ax.axvline(2.5,color="#bbbbbb",lw=0.8,ls=(0,(3,3)))
    ax.set_xticks(x); ax.set_xticklabels(xlabels,fontsize=7.2,rotation=12)
    ax.set_ylim(0,100); ax.set_yticks([0,25,50,75,100]); ax.tick_params(labelsize=7.5)
    ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0,0].set_ylabel("% of class guides   |   AA: % of loci")
axes[1,0].set_ylabel("% of class guides   |   AA: % of loci")

fig.text(0.5,0.975,"Human head-to-head  (CHM13 vs GRCh38)  +  AlleleAnalyzer discriminability",
         ha="center",fontsize=12.5,color="#2b2f3a",weight="bold")

h,l=[],[]
for ax in [axes[0,0]]:
    hh,ll=ax.get_legend_handles_labels(); h+=hh; l+=ll
fig.legend(h,l,loc="lower center",ncol=4,frameon=False,fontsize=8.6,bbox_to_anchor=(0.5,-0.01))
fig.tight_layout(rect=[0,0.035,1,0.945])
png=f"{OUT}/fig_h2h_human_with_alleleanalyzer.png"
fig.savefig(png,dpi=200,bbox_inches="tight"); plt.close()
print(f"wrote {png}")

# --- dump the plotted numbers as a TSV ---------------------------------------
tsv=f"{OUT}/headtohead_human_with_aa.tsv"
with open(tsv,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["class","truth","tool","metric","accept_or_covered","reject","not_generated","denominator","pct_accept_or_covered"])
    for cls in ORDER:
        N=HUMAN[cls]["N"]
        for t in tools:
            a,r,g=HUMAN[cls][t]
            w.writerow([cls,TRUTH[cls],t,"class-guide accept/reject",a,r,g,N,round(100*a/N,1) if N else 0])
        n,cov=AA[cls]
        w.writerow([cls,TRUTH[cls],"AlleleAnalyzer","locus discriminability",cov,"","",n,round(100*cov/n,1) if n else 0])
print(f"wrote {tsv}")
