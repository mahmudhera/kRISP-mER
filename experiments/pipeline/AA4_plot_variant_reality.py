#!/usr/bin/env python3
# AA4_plot_variant_reality.py
# The AlleleAnalyzer reference-blindness figure.
#
# AlleleAnalyzer is reference-VCF driven: per heterozygous SNP it designs BOTH a
# gRNA_ref (targets the reference allele) and a gRNA_alt (targets the individual's
# alt allele). Using the individual's own genome (CHM13 k=23 jellyfish) as ground
# truth, we ask, per class, how often each designed guide's allele ACTUALLY EXISTS
# in the individual.
#
# Controls (validated in AA3): the ref guide is present in GRCh38 ~96% (reconstruction
# correct) and the alt guide is mostly ABSENT from GRCh38 (it carries the variant).
#
# FINDING = ref_in_chm13: how often the REFERENCE allele AA anchors on is present in
# the individual. It is LOW, and tracks the class definition -- lowest for
# present_v_absent (sequence individual-only) where AA's reference guide is phantom
# ~90% of the time. This quantifies the cost of AA's reference dependence, using the
# same individual-genome ground truth kRISP-meR is built on.
#
# Outputs into allele-results/:
#   fig_aa_variant_reality_human.png
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

HERE="/mnt/windows-data/Thesis/allele-results"
TSV=f"{HERE}/aa_variant_reality_human.tsv"
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})

REFC="#C0576B"   # reference-allele bar (the finding: often absent in individual)
ALTC="#5B8C9E"   # alt-allele bar (individual's own allele)
CTRL="#B9BEC7"   # GRCh38 control marker

TITLE={"multiple_v_single":"Multiple in individual /\nSingle in reference (reject)",
"ref_present_v_absent":"Absent in individual /\nPresent in reference (reject)",
"ref_multiple_v_single":"Single in individual /\nMultiple in reference (accept)",
"present_v_absent":"Present in individual /\nAbsent in reference (accept)"}
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]

D={}
for r in csv.DictReader(open(TSV),delimiter="\t"):
    D[r["class"]]=r

fig,ax=plt.subplots(figsize=(10.6,5.8))
x=np.arange(len(ORDER)); w=0.38
ref_chm=[float(D[c]["ref_in_chm13_pct"]) for c in ORDER]
alt_chm=[float(D[c]["alt_in_chm13_pct"]) for c in ORDER]
ref_grc=[float(D[c]["ref_in_grch38_pct"]) for c in ORDER]

b1=ax.bar(x-w/2, ref_chm, w, color=REFC, edgecolor="white", lw=0.7,
          label="AA REFERENCE-allele guide present in individual (CHM13)")
b2=ax.bar(x+w/2, alt_chm, w, color=ALTC, edgecolor="white", lw=0.7,
          label="AA alt-allele guide present in individual (CHM13)")
# GRCh38 control tick for the ref guide (should sit ~96%): faint marker + hairline
for i,v in enumerate(ref_grc):
    ax.plot([x[i]-w/2-w*0.42, x[i]-w/2+w*0.42],[v,v],color="#6b6f77",lw=1.4,zorder=5)
ax.plot([],[],color="#6b6f77",lw=1.4,label="(control) ref guide present in GRCh38 ≈ reconstruction OK")

for rects in (b1,b2):
    for r in rects:
        h=r.get_height()
        ax.text(r.get_x()+r.get_width()/2, h+1.5, f"{h:.0f}%", ha="center", va="bottom",
                fontsize=8.5, color="#333")

ax.set_xticks(x); ax.set_xticklabels([TITLE[c] for c in ORDER], fontsize=8.4)
ax.set_ylim(0,105); ax.set_yticks([0,25,50,75,100])
ax.set_ylabel("% of AA SNP guides whose allele exists in the individual")
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.set_title("AlleleAnalyzer anchors guides on a reference allele the individual often lacks\n"
             "(reference-driven design vs. CHM13 individual-genome ground truth; SNP guides)",
             fontsize=11.5, color="#2b2f3a", weight="bold", pad=12)
ax.legend(loc="lower center", bbox_to_anchor=(0.5,-0.30), ncol=1, frameon=False, fontsize=8.8)
fig.subplots_adjust(bottom=0.30)
png=f"{HERE}/fig_aa_variant_reality_human.png"
fig.savefig(png,dpi=200,bbox_inches="tight"); plt.close()
print(f"wrote {png}")
