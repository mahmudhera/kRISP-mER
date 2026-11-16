#!/usr/bin/env python3
# plot_cfd_distributions.py — TASK 2 figures.
# For each (organism, anchor): 4 classes x 2 groups (accepted vs rejected-side),
# box + strip of CRISPOR cfdSpecScore (0-100, high = specific/safe).
# rejected-side = kRISP-meR (generated&rejected) UNION (not-generated).
# Reads Thesis/data/cfd_dist_<org>_<anchor>.tsv ; writes plots/cfd-distributions/.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

OUT="/mnt/windows-data/Thesis/plots/cfd-distributions"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#5B8C6E"; REJ="#C6603F"   # accepted (green), rejected-side (red-brown)
ORDER=["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TITLE={"multiple_v_single":"Multiple in individual / Single in ref\n(should reject)",
"ref_present_v_absent":"Absent in individual / Present in ref\n(should reject)",
"ref_multiple_v_single":"Single in individual / Multiple in ref\n(should accept)",
"present_v_absent":"Present in individual / Absent in ref\n(should accept)"}

COMBOS=[("rice","n22","Rice — targets from N22 (individual)"),
        ("rice","japonica","Rice — targets from japonica (reference)"),
        ("human","chm13","Human — targets from CHM13 (individual)"),
        ("human","grch38","Human — targets from GRCh38 (reference)")]

def load(org,anchor):
    p=f"/mnt/windows-data/Thesis/data/cfd_dist_{org}_{anchor}.tsv"
    if not os.path.exists(p): return None
    D={c:{"accepted":[],"rejected_side":[]} for c in ORDER}
    for r in csv.DictReader(open(p),delimiter="\t"):
        D[r["class"]][r["group"]].append(float(r["cfd_spec"]))
    return D

for org,anchor,title in COMBOS:
    D=load(org,anchor)
    if D is None: print(f"[skip] {org} {anchor} (no tsv)"); continue
    fig,axes=plt.subplots(1,4,figsize=(15,4.3))
    for ax,cls in zip(axes,ORDER):
        acc=D[cls]["accepted"]; rej=D[cls]["rejected_side"]
        data=[acc,rej]; labels=[f"accepted\n(n={len(acc)})",f"rejected+\nnot-gen\n(n={len(rej)})"]
        cols=[ACC,REJ]
        bp=ax.boxplot(data,positions=[0,1],widths=0.5,patch_artist=True,showfliers=False,
                      medianprops=dict(color="#222",lw=1.4))
        for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor(c)
        for i,(vals,c) in enumerate(zip(data,cols)):
            if vals:
                jit=np.random.default_rng(0).normal(i,0.06,size=len(vals))
                ax.scatter(jit,vals,s=9,color=c,alpha=0.5,edgecolors="none",zorder=3)
        ax.set_xticks([0,1]); ax.set_xticklabels(labels,fontsize=7.5)
        ax.set_ylim(-3,103); ax.set_yticks([0,25,50,75,100])
        ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
        ax.axhline(50,ls="--",lw=0.8,color="#999",zorder=1)  # CFD reject threshold
        for s in ("top","right"): ax.spines[s].set_visible(False)
    axes[0].set_ylabel("CRISPOR CFD specificity  (0–100, high = specific)")
    fig.suptitle(title+"   —   CFD specificity by kRISP-meR verdict",
                 fontsize=12,color="#2b2f3a",weight="bold",y=1.02)
    fig.text(0.5,-0.04,"dashed line = CFD reject threshold (spec<50).  "
             "accepted = kRISP-meR Ecuts<=1.5 ; rejected+not-gen = the rest.",
             ha="center",fontsize=8,color="#555")
    fig.tight_layout()
    p=f"{OUT}/cfd_{org}_{anchor}.png"
    fig.savefig(p,dpi=200,bbox_inches="tight"); plt.close()
    print(f"wrote {p}")
