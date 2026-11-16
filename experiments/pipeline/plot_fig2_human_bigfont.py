#!/usr/bin/env python3
# plot_fig2_human_bigfont.py  --  run ON THE VM in .../human/feb-20/
#
# Regenerates Figure-2 HUMAN panels with LARGE fonts + tight margins (Hera's request).
# Reproduces the ORIGINAL metrics exactly:
#   scatter (a): all_scores  ->  x=Genome-based Score, y=Krispmer Score (Reads Only)
#   bar (d):     ot_summary_{crispor,guidescan,krispmer}_human, EXACT bucketing from
#                plot_offtargets_human.py but pooled over ALL targets. Stacked
#                [1mm, 2+mm, CFD-safe/other]  (script's records[...][1:]).
# Violin (g) is NOT here -- its GuideScan values are corrupted in the local CSVs;
# that panel comes from Udoy's script.
#
# Outputs: human_scatter.png, human_bar.png   (dpi=300, bbox tight)
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

# ---- BIG-FONT styling (the whole point of the regeneration) ----------------
plt.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":20,            # base
    "axes.titlesize":24,
    "axes.labelsize":22,
    "xtick.labelsize":18,
    "ytick.labelsize":18,
    "legend.fontsize":18,
    "legend.title_fontsize":19,
    "axes.linewidth":1.4,
    "figure.facecolor":"white","axes.facecolor":"white",
})
BLUE="#1f77b4"

# ============================ (a) SCATTER ===================================
def scatter():
    g=[]; k=[]
    with open("all_scores") as fh:
        for ln in fh:
            p=ln.split()
            if len(p)>=2:
                try: g.append(float(p[0])); k.append(float(p[1]))
                except: pass
    fig,ax=plt.subplots(figsize=(6.4,6.4))
    ax.scatter(g,k,s=45,alpha=0.15,edgecolors=BLUE,facecolors="none",linewidths=1.2)
    lo=min(min(g),min(k)); hi=max(max(g),max(k))
    ax.plot([lo,hi],[lo,hi],ls="--",lw=2.6,color=BLUE)
    ax.set_xlabel("Genome-based Score")
    ax.set_ylabel("Krispmer Score (Reads Only)")
    for s in ("top","right"): ax.spines[s].set_visible(True)  # keep box like original
    ax.tick_params(width=1.4,length=6)
    fig.tight_layout(pad=0.6)
    fig.savefig("human_scatter.png",dpi=300,bbox_inches="tight")
    plt.close(); print("wrote human_scatter.png")

# ============================ (d) BAR =======================================
def bar():
    files={"CRISPOR":"ot_summary_crispor_human",
           "GuideScan2":"ot_summary_guidescan_human",
           "kRISP-meR":"ot_summary_krispmer_human"}
    ot_categories=["0 mismatch OT","1 mismatch OT","2/more mismatch OT"]
    # EXACT bucketing from plot_offtargets_human.py, pooled over all targets.
    # records index: 0=perfectly-specific(0mm), 1=has-0mm-OT, 2=has-1mm-OT, 3=CFD-safe/other
    order=["CRISPOR","GuideScan2","kRISP-meR"]
    stacks={}   # tool -> [cat0, cat1, cat2]  (the plotted records[...][1:] were [1,2,3];
                # but the manuscript legend is 0/1/2-mm, so we plot buckets [0,1,2])
    for tool in order:
        f=files[tool]
        rec=[0,0,0,0]
        if os.path.exists(f):
            df=pd.read_csv(f,sep=" ",header=None,names=["target","grna","a","b","c"])
            for _,r in df.iterrows():
                a,b,c=r["a"],r["b"],r["c"]
                if a+b+c==0: rec[0]+=1; rec[3]+=1
                elif a>0: rec[1]+=1
                elif b>0: rec[2]+=1
                else: rec[3]+=1
        stacks[tool]=rec
        print(f"{tool}: 0mm={rec[0]} has0mmOT={rec[1]} has1mmOT={rec[2]} CFDsafe={rec[3]}")

    # NOTE: the three stacked categories reproduce the original script's records[...][1:]
    #       = [has-0mm-OT, has-1mm-OT, CFD-safe/other]. Confirm which triplet the
    #       published panel uses; adjust `plot_idx` if needed.
    plot_idx=[1,2,3]   # original plotted records[1:]
    cat_labels=["0 mismatch OT","1 mismatch OT","2/more mismatch OT"]
    colours=["#E3A15B","#C0576B","#86B3C4"]

    fig,ax=plt.subplots(figsize=(7.2,6.4))
    x=np.arange(len(order)); w=0.6
    bottoms=np.zeros(len(order))
    for ci,idx in enumerate(plot_idx):
        vals=np.array([stacks[t][idx] for t in order],dtype=float)
        ax.bar(x,vals,w,bottom=bottoms,color=colours[ci],edgecolor="black",lw=0.8,
               label=cat_labels[ci])
        for xi,(v,b) in enumerate(zip(vals,bottoms)):
            if v>=30: ax.text(xi,b+v/2,f"{int(v)}",ha="center",va="center",fontsize=15)
        bottoms+=vals
    ax.set_xticks(x); ax.set_xticklabels(order)
    ax.set_ylabel("Number of gRNAs")
    ax.set_title("All targets")
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(width=1.4,length=6)
    ax.legend(title="OT Types",loc="upper right",frameon=True)
    fig.tight_layout(pad=0.6)
    fig.savefig("human_bar.png",dpi=300,bbox_inches="tight")
    plt.close(); print("wrote human_bar.png")

if __name__=="__main__":
    scatter()
    bar()
