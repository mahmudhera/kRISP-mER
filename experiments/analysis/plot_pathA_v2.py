#!/usr/bin/env python3
# plot_pathA_v2.py  [A|B|both]
#
# Journal-styled replacements for plot_pathA.py's figure.
#
# What was wrong with v1, and what changes here:
#   - v1 carried an embedded title and italic subtitle. Journal figures put that
#     in the caption; inside the axes it reads as a slide.
#   - v1 shaded the above-diagonal region pink and annotated it with an editorial
#     conclusion ("amplified in the individual -> under-counted by the reference").
#     A figure should show the effect, not instruct the reader.
#   - v1 drew a CURVED grey line from the HBB label to its point. A curve through
#     a scatter reads as a fit or trend; this one was only a pointer.
#   - Colours are now the manuscript's Figure 2 palette so the paper looks like
#     one document.
#
# A: the same scatter, cleaned up. Keeps the CIRCLE-seq read-count sizing, which
#    is what makes the HBB point land.
# B: paired dot plot of the elevated sites only. Reference and K562 copy number
#    as two joined points per site, sorted by fold change. Makes the discrepancy
#    the subject rather than the correlation.
#
# Reads data/k562_offtarget_copynumber.tsv, writes into data/ alongside it.
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WHICH = (sys.argv[1] if len(sys.argv) > 1 else "both").upper()
HERE  = os.path.dirname(os.path.abspath(__file__))
ROOT  = os.path.dirname(HERE)
TSV   = os.path.join(ROOT, "data", "k562_offtarget_copynumber.tsv")
OUT   = os.path.join(ROOT, "data")

# Manuscript palette (see figure-regen/v2/plot_fig2_panels_v2.py)
ELEV   = "#D96F6F"   # individual-elevated
CONC   = "#B9C0CB"   # concordant
ACCENT = "#3B7EA1"
INK    = "#2b2b2b"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans", "Arial"],
    "font.size": 17, "axes.labelsize": 18,
    "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15,
    "axes.linewidth": 1.3, "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": "#4a4a4a", "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": "#4a4a4a", "ytick.color": "#4a4a4a",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

def load():
    d = pd.read_csv(TSV, sep="\t")
    # same preparation as v1: keep the strongest cut per (site, off-target),
    # restrict to off-targets actually present in the reference
    d = (d.sort_values("circle_reads", ascending=False)
           .drop_duplicates(subset=["targetsite", "offtarget_seq"], keep="first"))
    d = d[(d.grch38_copy >= 1) & (d.grch38_copy <= 6)].copy()
    return d

def save(fig, stem):
    for ext in ("pdf", "png"):
        p = os.path.join(OUT, f"{stem}.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.pdf / .png")

# ------------------------------------------------------------------ A
def figure_a(d):
    YMAX, XMAX = 12.0, 6.0
    d = d.copy()
    d["yv"] = d.k562_copy.clip(upper=YMAX)
    rng = np.random.default_rng(7)
    d["xv"] = d.grch38_copy + rng.uniform(-0.13, 0.13, len(d))
    d["ms"] = 26 + d.circle_reads.clip(upper=350) / 2.2

    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    ax.plot([0, YMAX], [0, YMAX], ls=(0, (6, 5)), lw=1.6, color="#888888", zorder=1)
    ax.text(0.97, 0.03, "equal copy number", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=14, color="#777777", style="italic")

    for cls, colour, z in [("concordant", CONC, 2), ("individual_elevated", ELEV, 3)]:
        s = d[d["class"] == cls]
        lbl = ("concordant" if cls == "concordant" else "elevated in K562")
        ax.scatter(s.xv, s.yv, s=s.ms, c=colour, alpha=0.75,
                   edgecolors="white", linewidths=0.7, zorder=z,
                   label=f"{lbl}  ($n$ = {len(s)})")

    # HBB: the strongest measured cut of that guide. Straight leader, no curve.
    hbb = d.loc[d.circle_reads.idxmax()]
    ax.annotate(f"{int(hbb.circle_reads):,} CIRCLE-seq reads\n"
                f"GRCh38 {hbb.grch38_copy:.0f}$\\times$, K562 {hbb.k562_copy:.1f}$\\times$",
                xy=(hbb.grch38_copy, min(hbb.k562_copy, YMAX)),
                xytext=(2.35, 7.4), fontsize=14, color=INK,
                ha="left", va="center", linespacing=1.25,
                arrowprops=dict(arrowstyle="-", lw=1.1, color="#666666",
                                shrinkA=0, shrinkB=9), zorder=6)

    ax.set_xlim(0.45, XMAX + 0.45); ax.set_ylim(-0.4, YMAX + 0.4)
    ax.set_xlabel("Copy number in GRCh38")
    ax.set_ylabel("Copy number in K562 (from reads)")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.tick_params(width=1.3, length=6)
    ax.grid(True, ls=":", lw=0.7, color="#cfcfcf", zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, handletextpad=0.4, borderpad=0.2)

    fig.tight_layout(pad=0.5)
    save(fig, "k562_copynumber_scatter")

# ------------------------------------------------------------------ B
def figure_b(d):
    e = d[d["class"] == "individual_elevated"].copy()
    e["fold"] = e.k562_copy / e.grch38_copy.clip(lower=1)
    e = e.sort_values("fold").reset_index(drop=True)
    y = np.arange(len(e))

    fig, ax = plt.subplots(figsize=(7.2, max(6.0, 0.20 * len(e))))
    ax.hlines(y, e.grch38_copy, e.k562_copy, color="#cbcbcb", lw=1.6, zorder=1)
    ax.scatter(e.grch38_copy, y, s=42, c=CONC, edgecolors="white", linewidths=0.7,
               zorder=3, label="GRCh38 (reference)")
    ax.scatter(e.k562_copy, y, s=30 + e.circle_reads.clip(upper=350) / 3.0,
               c=ELEV, edgecolors="white", linewidths=0.7, zorder=3,
               label="K562 (from reads)")

    hi = e.circle_reads.idxmax()
    ax.annotate(f"strongest measured cut\n({int(e.loc[hi,'circle_reads']):,} reads)",
                xy=(e.loc[hi, "k562_copy"], hi),
                xytext=(e.k562_copy.max() * 0.62, len(e) * 0.33),
                fontsize=14, color=INK, ha="left", va="center", linespacing=1.25,
                arrowprops=dict(arrowstyle="-", lw=1.1, color="#666666",
                                shrinkA=0, shrinkB=8), zorder=6)

    ax.set_yticks([])
    ax.set_xlabel("Copy number")
    ax.set_ylabel(f"Experimentally confirmed off-targets ($n$ = {len(e)})", labelpad=12)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.tick_params(width=1.3, length=6)
    ax.grid(True, axis="x", ls=":", lw=0.7, color="#cfcfcf", zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False, handletextpad=0.4)

    fig.tight_layout(pad=0.5)
    save(fig, "k562_copynumber_paired")

if __name__ == "__main__":
    if not os.path.exists(TSV):
        sys.exit(f"missing {TSV}\n"
                 "Copy k562_offtarget_copynumber.tsv from the VM "
                 "(~/pathA_circleseq/) into data/ first.")
    d = load()
    print(f"{len(d)} off-targets present in the reference (GRCh38 copy 1-6)")
    print(d["class"].value_counts().to_string())
    if WHICH in ("A", "BOTH"): figure_a(d)
    if WHICH in ("B", "BOTH"): figure_b(d)
