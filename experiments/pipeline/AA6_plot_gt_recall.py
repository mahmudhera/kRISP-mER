#!/usr/bin/env python3
# AA6_plot_gt_recall.py
# kRISP-meR vs AlleleAnalyzer RECALL of the individual's real targetable guides.
#
# Ground truth (AA5): per target window, the NGG guides present in the CHM13 assembly
# but ABSENT from GRCh38 = real, individual-specific targetable guides. Recall = of
# that set, how many each tool designed (matched on canonical 20-mer protospacer).
#
# SCOPE (decided): only the two SINGLE-IN-REFERENCE classes are reported, where the
# individual's locus is unambiguous and the window->CHM13 mapping is clean:
#   multiple_v_single (173 GT) and present_v_absent (468 GT).
# The two ref_* (multiple/present-in-reference) classes are EXCLUDED: the reference has
# several copies so "the individual's copy" is ambiguous, and the GT window cannot be
# guaranteed to be the same genomic copy kRISP-meR targeted -> not measurable this way.
#
# Left panel: the two clean classes. Right panel: pooled (100 targets, 641 GT guides).
# Writes allele-results/fig_gt_recall_human.png
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os

HERE="/mnt/windows-data/Thesis/allele-results"; os.makedirs(HERE,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
KR="#4E8FA6"; AAC="#C0576B"

# class -> (|GT|, krisp_recovered, aa_recovered)
CLEAN={
 "Multiple in individual /\nSingle in reference":(173,161,1),
 "Present in individual /\nAbsent in reference":(468,459,0),
}
POOL=(641,620,1)   # sum of the two clean classes

fig,(axL,axR)=plt.subplots(1,2,figsize=(11.8,5.6),gridspec_kw={"width_ratios":[2,1]})

# ---- left: two clean classes ----
labels=list(CLEAN); x=np.arange(len(labels)); w=0.38
kp=[100*CLEAN[c][1]/CLEAN[c][0] for c in labels]
ap=[100*CLEAN[c][2]/CLEAN[c][0] for c in labels]
axL.bar(x-w/2,kp,w,color=KR,edgecolor="white",lw=0.7,label="kRISP-meR (reads-based)")
axL.bar(x+w/2,ap,w,color=AAC,edgecolor="white",lw=0.7,label="AlleleAnalyzer (reads→GRCh38→VCF)")
for i,c in enumerate(labels):
    axL.text(i-w/2,kp[i]+1.5,f"{kp[i]:.0f}%",ha="center",va="bottom",fontsize=9,color="#333")
    axL.text(i+w/2,ap[i]+1.5,f"{ap[i]:.0f}%",ha="center",va="bottom",fontsize=9,color="#333")
    axL.text(i,-19,f"GT n={CLEAN[c][0]}",ha="center",va="top",fontsize=8,color="#666")
axL.set_xticks(x); axL.set_xticklabels(labels,fontsize=8.6)
axL.set_ylim(0,105); axL.set_yticks([0,25,50,75,100])
axL.set_ylabel("% of individual-specific ground-truth guides recovered")
axL.set_title("By class (single-in-reference only)",fontsize=10.5,color="#2b2f3a")
for s in ("top","right"): axL.spines[s].set_visible(False)
axL.legend(loc="upper center",bbox_to_anchor=(0.5,-0.18),ncol=2,frameon=False,fontsize=8.8)

# ---- right: pooled ----
n,k,a=POOL
xk=np.arange(2)
vals=[100*k/n,100*a/n]
axR.bar(xk,vals,0.55,color=[KR,AAC],edgecolor="white",lw=0.7)
for i,v in enumerate(vals):
    axR.text(i,v+1.5,f"{v:.0f}%",ha="center",va="bottom",fontsize=11,weight="bold",color="#333")
axR.set_xticks(xk); axR.set_xticklabels(["kRISP-meR","AlleleAnalyzer"],fontsize=9)
axR.set_ylim(0,105); axR.set_yticks([0,25,50,75,100])
axR.text(0.5,-8,f"pooled: 100 targets, {n} GT guides",ha="center",va="top",fontsize=8,color="#666",transform=axR.get_xaxis_transform())
axR.set_title("Pooled",fontsize=10.5,color="#2b2f3a")
for s in ("top","right"): axR.spines[s].set_visible(False)

fig.suptitle("Recovery of the individual's real targetable guides\n"
             "(ground truth = NGG guides in the CHM13 assembly but absent from GRCh38)",
             fontsize=12.5,color="#2b2f3a",weight="bold")
fig.subplots_adjust(top=0.84, bottom=0.24, left=0.09, right=0.97, wspace=0.28)
png=f"{HERE}/fig_gt_recall_human.png"; fig.savefig(png,dpi=200); plt.close()
print(f"wrote {png}")
