#!/usr/bin/env python3
# plot_reads_cfd_cuts.py -- TASK 2 (variant b) figures, 3-tool version.
# READS-BASED CFD off-target load (manuscript quantity): CFD-weighted expected
# genomic cuts = (1/lambda) * sum reads(neighbour)*CFD over Hamming<=2 neighbours.
# HIGH = DANGEROUS (opposite polarity to cfdSpecScore).
# Per (organism, anchor): 4 classes, each with 3 tools x {accepted, rejected+not-gen}
# = 6 boxes; box + strip, log-y. The cfd_cuts value is a property of the guide
# (same across tools); only the accept/reject split differs per tool.
# Reads Thesis/data/reads_cfd_dist_<org>_<anchor>.tsv ; writes plots/cfd-distributions/.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv
from collections import defaultdict

OUT="/mnt/windows-data/Thesis/plots/cfd-distributions"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#5B8C6E"; REJ="#C6603F"   # accepted (green), rejected+not-gen (red-brown)
TOOLS=["kRISP-meR","CRISPOR","GuideScan2"]
GROUPS=[("accepted",ACC),("rejected_side",REJ)]
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
    p=f"/mnt/windows-data/Thesis/data/reads_cfd_dist_{org}_{anchor}.tsv"
    if not os.path.exists(p): return None
    D={c:{t:{"accepted":[],"rejected_side":[]} for t in TOOLS} for c in ORDER}
    for r in csv.DictReader(open(p),delimiter="\t"):
        if r["class"] in D and r["tool"] in TOOLS:
            D[r["class"]][r["tool"]][r["group"]].append(float(r["cfd_cuts"]))
    return D

EPS=1e-2   # log-axis floor so zeros are visible
# x positions: 3 tools, each a tight (accepted, rejected) pair
POS={}; xt=[]; xl=[]
x=0
for t in TOOLS:
    POS[(t,"accepted")]=x; POS[(t,"rejected_side")]=x+0.8
    xt.append(x+0.4); xl.append(t)
    x+=2.4

YCAP=15.0   # linear y-ceiling; guides above are drawn at the cap with a ^ marker
for org,anchor,title in COMBOS:
    D=load(org,anchor)
    if D is None: print(f"[skip] {org} {anchor} (no tsv)"); continue
    fig,axes=plt.subplots(1,4,figsize=(18,4.6),sharey=True)
    for ax,cls in zip(axes,ORDER):
        for t in TOOLS:
            for grp,col in GROUPS:
                vals=D[cls][t][grp]
                pos=POS[(t,grp)]
                if vals:
                    # box on the clipped values so it stays inside the linear frame
                    vclip=[min(v,YCAP) for v in vals]
                    bp=ax.boxplot([vclip],positions=[pos],widths=0.62,patch_artist=True,
                                  showfliers=False,medianprops=dict(color="#222",lw=1.3))
                    bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.35); bp["boxes"][0].set_edgecolor(col)
                    rng=np.random.default_rng(0)
                    inside=[v for v in vals if v<=YCAP]; over=[v for v in vals if v>YCAP]
                    if inside:
                        ax.scatter(rng.normal(pos,0.07,size=len(inside)),inside,
                                   s=7,color=col,alpha=0.5,edgecolors="none",zorder=3)
                    if over:   # clipped outliers pinned at the ceiling with a caret
                        ax.scatter(rng.normal(pos,0.07,size=len(over)),
                                   [YCAP*0.985]*len(over),s=16,marker="^",
                                   color=col,alpha=0.8,edgecolors="none",zorder=4)
                        ax.text(pos,YCAP*1.01,f"+{len(over)}",ha="center",va="bottom",
                                fontsize=5.5,color=col)
                    # mean marker (accept/reject gap lives in the mean); clip its position too
                    mval=float(np.mean(vals))
                    ax.scatter([pos],[min(mval,YCAP)],marker="D",s=26,facecolor="white",
                               edgecolor=col,linewidths=1.4,zorder=5)
                elif grp=="accepted":
                    ax.text(pos,0.4,"designs\nnone",ha="center",va="bottom",
                            fontsize=6,color="#999",style="italic")
        ax.set_xticks(xt); ax.set_xticklabels(xl,fontsize=8)
        ax.set_ylim(0,YCAP*1.06)
        ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
        ax.axhline(1.5,ls="--",lw=0.8,color="#999",zorder=1)
        for s in ("top","right"): ax.spines[s].set_visible(False)
    axes[0].set_ylabel(f"Reads-based CFD-weighted expected cuts  (high = dangerous; capped at {YCAP:g})")
    # legend
    import matplotlib.patches as mp
    from matplotlib.lines import Line2D
    mean_h=Line2D([0],[0],marker="D",color="none",markerfacecolor="white",
                  markeredgecolor="#555",markersize=7,label="mean")
    fig.legend([mp.Patch(color=ACC,alpha=0.5),mp.Patch(color=REJ,alpha=0.5),mean_h],
               ["accepted (by that tool)","rejected + not-generated","mean (♦)"],
               loc="upper right",fontsize=8,frameon=False,bbox_to_anchor=(0.998,1.04))
    fig.suptitle(title+"   —   reads-CFD off-target load, per tool by its own verdict",
                 fontsize=12,color="#2b2f3a",weight="bold",y=1.03,x=0.42)
    fig.text(0.5,-0.05,f"CFD-weighted expected genomic cuts = (1/λ)·Σ reads(m)·CFD(g,m) over Hamming≤2 neighbours m (both strands).  "
             f"dashed = 1.5 cuts.  ▲+n = n guides above the {YCAP:g} cap.  reject: kRISP-meR Ecuts>1.5 | CRISPOR MITspec<50 | GuideScan2 spec<0.2.",
             ha="center",fontsize=7.5,color="#555")
    fig.tight_layout(rect=(0,0,1,0.97))
    p=f"{OUT}/reads_cfd_{org}_{anchor}.png"
    fig.savefig(p,dpi=200,bbox_inches="tight"); plt.close()
    print(f"wrote {p}")
