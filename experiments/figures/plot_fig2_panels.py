#!/usr/bin/env python3
# plot_fig2_panels.py  ORG        (run locally: python3 plot_fig2_panels.py human)
#
# Regenerates Figure-2 panels with LARGE fonts + tight margins (Hera's request),
# reproducing the ORIGINAL metrics exactly. Reads data from ./<ORG>/ and writes
# <ORG>_scatter.png and <ORG>_bar.png into ./<ORG>/.
#   scatter : all_scores  ->  x=Genome-based Score, y=Krispmer Score (Reads Only)
#   bar     : ot_summary_{crispor,guidescan,krispmer}_<org>, EXACT bucketing from
#             plot_offtargets_human.py, pooled over ALL targets.
# Violin comes from Udoy's script (GuideScan values corrupted in local CSVs).
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

ORG = sys.argv[1] if len(sys.argv) > 1 else "human"
D   = os.path.join(os.path.dirname(os.path.abspath(__file__)), ORG)
ORG_LABEL = {"human":"Human","mouse":"Mouse","yeast":"Yeast"}.get(ORG, ORG.capitalize())

# ---- BIG-FONT styling; Liberation Sans = metric-compatible Arial clone --------
plt.rcParams.update({
    "font.family":"sans-serif",
    "font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":20, "axes.titlesize":24, "axes.labelsize":22,
    "xtick.labelsize":18, "ytick.labelsize":18,
    "legend.fontsize":17, "legend.title_fontsize":18,
    "axes.linewidth":1.4, "figure.facecolor":"white", "axes.facecolor":"white",
})
BLUE="#1f77b4"

def save_both(fig, stem):
    """Save a panel as vector PDF (zoom-proof, for the manuscript) AND 300-dpi PNG."""
    pdf=os.path.join(D,f"{stem}.pdf"); png=os.path.join(D,f"{stem}.png")
    fig.savefig(pdf,bbox_inches="tight")           # vector -- resolution independent
    fig.savefig(png,dpi=300,bbox_inches="tight")   # raster preview
    plt.close(fig)
    print(f"wrote {pdf}\nwrote {png}")

# ============================ SCATTER ======================================
def scatter():
    # scatter data: a single all_scores, OR (yeast) all_scores_coding + all_scores_NC merged.
    srcs=[c for c in ["all_scores","all_scores_coding","all_scores_NC"] if os.path.exists(os.path.join(D,c))]
    if "all_scores" in srcs: srcs=["all_scores"]     # prefer the merged file if present
    g=[]; k=[]
    for s in srcs:
        with open(os.path.join(D,s)) as fh:
            for ln in fh:
                p=ln.split()
                if len(p)>=2:
                    try: g.append(float(p[0])); k.append(float(p[1]))
                    except: pass
    print(f"  scatter sources: {srcs}")
    fig,ax=plt.subplots(figsize=(6.4,6.4))
    ax.scatter(g,k,s=45,alpha=0.15,edgecolors=BLUE,facecolors="none",linewidths=1.2)
    lo=min(min(g),min(k)); hi=max(max(g),max(k))
    ax.plot([lo,hi],[lo,hi],ls="--",lw=2.6,color=BLUE)
    ax.set_xlabel("Genome-based Score")
    ax.set_ylabel("Krispmer Score (Reads Only)")
    ax.tick_params(width=1.4,length=6)
    fig.tight_layout(pad=0.6)
    save_both(fig, f"{ORG}_scatter")
    print(f"  (scatter n={len(g)})")

# ============================ BAR (two-bar grouped) =========================
# Reproduces the published d-panel: per tool, TWO bars side by side:
#   Bar 1 (stacked): gRNA counts by worst off-target mismatch level
#                    [0 mismatch OT, 1 mismatch OT, 2/more mismatch OT]  (mutually exclusive)
#   Bar 2 (single):  # gRNAs with CFD-safe (Expected_Cuts <= 1.5), + % annotation
# Guide set = the CFD summary file (authoritative; carries Expected_Cuts).
# NOTE: GuideScan bar-1 total here = 2007 (this data) vs 1951 in the published figure
#       (the published GuideScan set was slightly more filtered); documented, not fudged.
CFD_CUT=1.5
def _find(*cands):
    """return the first existing filename (handles .csv vs no-extension per organism)."""
    for c in cands:
        p=os.path.join(D,c)
        if os.path.exists(p): return c
    return cands[0]  # fall through (will error visibly if truly missing)

def bar():
    # GuideScan MUST use the *_filtered file (that is what the published figure used;
    # human filtered = 1951 rows, matching the paper). CRISPOR/kRISP use the plain summary.
    ot_files={
        "CRISPOR":   _find(f"ot_summary_crispor_{ORG}.csv", f"ot_summary_crispor_{ORG}"),
        "GuideScan2":_find(f"ot_summary_guidescan_{ORG}_filtered.csv", f"ot_summary_guidescan_{ORG}_filtered"),
        "kRISP-meR": _find(f"ot_summary_krispmer_{ORG}.csv", f"ot_summary_krispmer_{ORG}"),
    }
    cfd_files={
        "CRISPOR":   _find(f"cfd_score_summary_crispor_{ORG}.csv", f"cfd_score_summary_crispor_{ORG}_1_50.csv", f"cfd_score_summary_crispor_{ORG}"),
        "GuideScan2":_find(f"cfd_score_summary_guidescan_{ORG}.csv", f"cfd_score_summary_guidescan_{ORG}_1_50.csv", f"cfd_score_summary_guidescan_{ORG}"),
        "kRISP-meR": _find(f"cfd_score_summary_krispmer_{ORG}.csv", f"cfd_score_summary_krispmer_{ORG}_1_50.csv", f"cfd_score_summary_krispmer_{ORG}"),
    }
    order=["CRISPOR","GuideScan2","kRISP-meR"]

    mm={}       # tool -> [0mm, 1mm, 2+mm]  (mutually exclusive by worst OT)
    cfdsafe={}  # tool -> (n_safe, pct)
    for tool in order:
        # bar 1: mismatch categories from ot_summary. VERIFIED against published segments:
        #   "0 mismatch OT"      = a>0            (has a 0-mismatch off-target)
        #   "1 mismatch OT"      = a==0 & b>0     (worst OT is 1 mismatch)
        #   "2/more mismatch OT" = everything else (2+mm OR perfectly clean) -> total = all rows
        of=os.path.join(D,ot_files[tool]); rec=[0,0,0]
        if os.path.exists(of):
            df=pd.read_csv(of,sep=" ",header=None,names=["target","grna","a","b","c"])
            n=len(df)
            rec[0]=int((df["a"]>0).sum())
            rec[1]=int(((df["a"]==0)&(df["b"]>0)).sum())
            rec[2]=n-rec[0]-rec[1]
        mm[tool]=rec
        # bar 2: CFD-safe from cfd summary. DEDUP by guide first -- some GuideScan CFD
        # files repeat each guide once per off-target site (human GS: 4 dups; yeast GS:
        # ~5x, which inflated the published 5012). Each guide has a single Expected_Cuts
        # across its rows (verified), so keep-first == min == max; deduping is lossless.
        cf=os.path.join(D,cfd_files[tool])
        if os.path.exists(cf):
            cd=pd.read_csv(cf).drop_duplicates("gRNA")
            n=len(cd); safe=int((cd["Expected_Cuts"]<=CFD_CUT).sum())
            cfdsafe[tool]=(safe, 100*safe/n if n else 0)
        else:
            cfdsafe[tool]=(0,0)
        print(f"  {tool}: mismatch[0mm,1mm,2+mm]={rec} (sum={sum(rec)})  CFD<=1.5={cfdsafe[tool][0]} ({cfdsafe[tool][1]:.1f}%)")

    # ---- draw ----
    cat_labels=["0 mismatch OT","1 mismatch OT","2/more mismatch OT"]
    mm_colours=["#5B2A4A","#C0576B","#E08152"]   # dark / red / orange (published palette)
    cfd_colour="#F0C9A8"                          # peach
    fig,ax=plt.subplots(figsize=(7.6,6.4))
    x=np.arange(len(order)); w=0.36
    # bar 1 (stacked mismatch) at x-w/2
    bottoms=np.zeros(len(order))
    for ci in range(3):
        vals=np.array([mm[t][ci] for t in order],dtype=float)
        ax.bar(x-w/2,vals,w,bottom=bottoms,color=mm_colours[ci],edgecolor="black",lw=0.6,label=cat_labels[ci])
        bottoms+=vals
    # total labels above bar 1
    for xi,t in enumerate(order):
        tot=sum(mm[t]); ax.text(xi-w/2, tot+ax.get_ylim()[1]*0.01, f"{tot}", ha="center", va="bottom", fontsize=14)
    # bar 2 (CFD-safe) at x+w/2
    safe_vals=np.array([cfdsafe[t][0] for t in order],dtype=float)
    ax.bar(x+w/2,safe_vals,w,color=cfd_colour,edgecolor="black",lw=0.6,label="CFD ≤ 1.5")
    for xi,t in enumerate(order):
        s,p=cfdsafe[t]; ax.text(xi+w/2, s+ax.get_ylim()[1]*0.01, f"{s} ({p:.1f}%)", ha="center", va="bottom", fontsize=12)

    ax.set_xticks(x); ax.set_xticklabels(order)
    ax.set_ylabel("Number of gRNAs")   # panel title omitted; described in the caption / column headers
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(width=1.4,length=6)
    ax.legend(loc="upper right",frameon=True,fontsize=14)
    # headroom so the total/CFD% labels above bars don't collide with the top spine
    ax.set_ylim(0, ax.get_ylim()[1]*1.12)
    fig.tight_layout(pad=0.6)
    save_both(fig, f"{ORG}_bar")

# ============================ VIOLIN (rows g,h,i) ===========================
# Ported from Udoy's trial.py: per guide total_ot = ot0+ot1+ot2; y = log1p(total_ot);
# seaborn violin per tool. GuideScan uses the *_filtered file (consistent with the bar).
def violin():
    import seaborn as sns
    ot_files={
        "CRISPOR":   _find(f"ot_summary_crispor_{ORG}.csv", f"ot_summary_crispor_{ORG}"),
        "GuideScan2":_find(f"ot_summary_guidescan_{ORG}_filtered.csv", f"ot_summary_guidescan_{ORG}_filtered"),
        "kRISP-meR": _find(f"ot_summary_krispmer_{ORG}.csv", f"ot_summary_krispmer_{ORG}"),
    }
    order=["CRISPOR","GuideScan2","kRISP-meR"]
    palette={"kRISP-meR":"#A8DADC","CRISPOR":"#B5C99A","GuideScan2":"#F4B5C1"}  # trial.py colours
    frames=[]
    for tool in order:
        f=os.path.join(D,ot_files[tool])
        if not os.path.exists(f): print(f"  [violin] missing {ot_files[tool]}"); continue
        df=pd.read_csv(f,sep=r"\s+",header=None,names=["target","grna","ot0","ot1","ot2"])
        df["tool"]=tool; frames.append(df)
    alldf=pd.concat(frames,ignore_index=True)
    alldf["total_ot"]=alldf["ot0"]+alldf["ot1"]+alldf["ot2"]
    alldf["log_total_ot"]=np.log1p(alldf["total_ot"])

    fig,ax=plt.subplots(figsize=(7.2,6.0))
    sns.violinplot(data=alldf,x="tool",y="log_total_ot",order=order,hue="tool",
                   hue_order=order,palette=palette,inner=None,cut=0,linewidth=1.2,
                   legend=False,ax=ax)
    ax.set_ylabel("log(1 + total off-targets)")
    ax.set_xlabel("")   # panel title omitted; described in the caption
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(width=1.4,length=6)
    fig.tight_layout(pad=0.6)
    save_both(fig, f"{ORG}_violin")

if __name__=="__main__":
    print(f"[{ORG}] data dir: {D}")
    scatter()
    bar()
    try: violin()
    except Exception as e: print(f"  [violin] skipped: {e}")
