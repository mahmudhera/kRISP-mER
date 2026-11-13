#!/usr/bin/env python3
# plot_asm_cfd_cuts.py -- TASK 2 (variant c) figures, ASSEMBLY-based CFD load.
# CFD off-target load from the INDIVIDUAL'S ASSEMBLY k-mer counts (the source that
# DEFINES the 4 classes) -- NO reads, NO /lambda. Value = "expected genomic cuts":
# single-copy safe guide = 1.0, duplicated = 2.0+, dangerous = tens+. HIGH=DANGEROUS.
#
# ACCEPTED and REJECTED are drawn as SEPARATE figures (very different y-ranges):
#   *_ACCEPTED.png  -> tight linear 0-3 (the single-copy=1 vs duplicated=2 story)
#   *_REJECTED.png  -> wide linear, capped, log-friendly ceiling (the danger tail)
# Each: 4 classes x 3 tools. Reads Thesis/data/asm_cfd_dist_<org>_<anchor>.tsv.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, csv

OUT="/mnt/windows-data/Thesis/plots/cfd-distributions"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})
ACC="#5B8C6E"; REJ="#C6603F"
TOOLS=["kRISP-meR","CRISPOR","GuideScan2"]
TOOLCOL={"kRISP-meR":"#3F6D54","CRISPOR":"#4E79A7","GuideScan2":"#B07AA1"}
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
    p=f"/mnt/windows-data/Thesis/data/asm_cfd_dist_{org}_{anchor}.tsv"
    if not os.path.exists(p): return None
    D={c:{t:{"accepted":[],"rejected_side":[]} for t in TOOLS} for c in ORDER}
    for r in csv.DictReader(open(p),delimiter="\t"):
        if r["class"] in D and r["tool"] in TOOLS:
            D[r["class"]][r["tool"]][r["group"]].append(float(r["cfd_cuts"]))
    return D

# tool x positions within a class panel
POS={t:i for i,t in enumerate(TOOLS)}

def draw(D,title,group,ycap,fname,ylabel,reflines,note):
    fig,axes=plt.subplots(1,4,figsize=(15,4.4),sharey=True)
    for ax,cls in zip(axes,ORDER):
        for t in TOOLS:
            vals=D[cls][t][group]; pos=POS[t]; c=TOOLCOL[t]
            if not vals:
                if group=="accepted":
                    ax.text(pos,ycap*0.06,"designs\nnone",ha="center",va="bottom",
                            fontsize=6.5,color="#999",style="italic")
                continue
            vclip=[min(v,ycap) for v in vals]
            bp=ax.boxplot([vclip],positions=[pos],widths=0.55,patch_artist=True,
                          showfliers=False,medianprops=dict(color="#222",lw=1.4))
            bp["boxes"][0].set_facecolor(c); bp["boxes"][0].set_alpha(0.30); bp["boxes"][0].set_edgecolor(c)
            rng=np.random.default_rng(0)
            inside=[v for v in vals if v<=ycap]; over=[v for v in vals if v>ycap]
            if inside:
                ax.scatter(rng.normal(pos,0.06,size=len(inside)),inside,s=8,color=c,
                           alpha=0.5,edgecolors="none",zorder=3)
            if over:
                ax.scatter(rng.normal(pos,0.06,size=len(over)),[ycap*0.985]*len(over),
                           s=16,marker="^",color=c,alpha=0.85,edgecolors="none",zorder=4)
                ax.text(pos,ycap*1.008,f"+{len(over)}",ha="center",va="bottom",fontsize=6,color=c)
            ax.scatter([pos],[min(float(np.mean(vals)),ycap)],marker="D",s=28,
                       facecolor="white",edgecolor=c,linewidths=1.5,zorder=5)
            # n= just under the axis line, below the tick labels (annotate in axis coords)
            ax.annotate(f"n={len(vals)}",xy=(pos,0),xycoords=("data","axes fraction"),
                        xytext=(0,-24),textcoords="offset points",
                        ha="center",va="top",fontsize=6.5,color="#666")
        ax.set_xticks(range(len(TOOLS))); ax.set_xticklabels(TOOLS,fontsize=8)
        ax.set_ylim(0,ycap*1.05)
        ax.set_title(TITLE[cls],fontsize=8.5,color="#2b2f3a")
        for yv,st,cc in reflines: ax.axhline(yv,ls=st,lw=0.9,color=cc,zorder=1)
        for s in ("top","right"): ax.spines[s].set_visible(False)
    axes[0].set_ylabel(ylabel)
    tag="ACCEPTED (by each tool)" if group=="accepted" else "REJECTED + not-generated"
    fig.suptitle(f"{title}   —   assembly-CFD load, {tag}",
                 fontsize=12,color="#2b2f3a",weight="bold",y=1.03)
    fig.text(0.5,-0.06,note,ha="center",fontsize=7.3,color="#555")
    fig.tight_layout(rect=(0,0,1,0.97))
    fig.savefig(fname,dpi=200,bbox_inches="tight"); plt.close()
    print(f"wrote {fname}")

REJ_CAP=40.0
for org,anchor,title in COMBOS:
    D=load(org,anchor)
    if D is None: print(f"[skip] {org} {anchor} (no tsv)"); continue
    base=f"{OUT}/asm_cfd_{org}_{anchor}"
    # ACCEPTED — 0-8: single-copy(1) vs duplicated(2) story, with headroom for mild outliers
    draw(D,title,"accepted",ycap=8.0,fname=f"{base}_ACCEPTED.png",
         ylabel="Assembly CFD expected cuts  (accepted guides)",
         reflines=[(1.0,":","#3a7"),(2.0,"--","#999")],
         note="Guides each tool ACCEPTED.  dotted=1 (single copy, safe), dashed=2 (duplicated).  "
              "▲+n = accepted guides above the 8-cut cap.  reject rules: kRISP-meR Ecuts>1.5 | CRISPOR MITspec<50 | GuideScan2 spec<0.2.")
    # REJECTED — wide, capped: the danger tail
    draw(D,title,"rejected_side",ycap=REJ_CAP,fname=f"{base}_REJECTED.png",
         ylabel=f"Assembly CFD expected cuts  (rejected + not-gen; cap {REJ_CAP:g})",
         reflines=[(1.0,":","#3a7"),(2.0,"--","#999")],
         note="Guides each tool REJECTED or did not generate.  dotted=1, dashed=2.  "
              f"▲+n = guides above the {REJ_CAP:g}-cut cap (load reaches hundreds).  same reject rules as the accepted panel.")
