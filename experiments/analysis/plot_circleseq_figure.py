#!/usr/bin/env python3
"""
plot_circleseq_figure.py — CIRCLE-seq validation figure for the kRISP-meR manuscript.

Produces the paired dot plot: for each experimentally confirmed off-target that is
amplified in K562 relative to GRCh38, the reference copy number and the read-derived
K562 copy number are drawn as two joined points. The gap between them is the
quantity the figure is about — how far the reference under-counts a real, measured
cleavage site.

Design decisions, so they are not re-litigated later:

  * Every visual channel encodes data. Colour = which guide the off-target belongs
    to (six guides); marker area = CIRCLE-seq cleavage read count; row order =
    absolute copy-number gain. Nothing is decorative.
  * Rows are ordered by ABSOLUTE gain (k562 - grch38), not fold change. Fold change
    puts a 3x -> 5.8x site below a 1x -> 2.1x one, which makes the ladder appear to
    invert near the bottom and reads as an error.
  * No title inside the axes, no shaded regions, no editorial annotation. Those
    belong in the caption. An earlier version carried all three and read as a
    conference slide rather than a journal figure.
  * The legend is placed outside the data area. In the first draft it overlapped
    the lowest three rows.
  * Palette is colour-blind safe (Okabe-Ito derived) and prints legibly in
    greyscale, since the guides are also separated by row position.

Usage:
    python3 plot_circleseq_figure.py            # writes into ../data/
    python3 plot_circleseq_figure.py --outdir X

Input: data/k562_offtarget_copynumber.tsv
    columns used: targetsite, offtarget_seq, circle_reads, n_mm,
                  grch38_copy, k562_copy, class
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Okabe-Ito derived; distinguishable in greyscale by lightness as well as hue.
GUIDE_COLOURS = {
    "HBB":    "#C2453B",   # deep red   -- carries the headline site
    "EMX1":   "#1F6F8B",   # teal blue
    "FANCF":  "#D98C2B",   # amber
    "RNF2":   "#6B8E4E",   # olive
    "Site_1": "#7A5C9E",   # violet
    "Site_2": "#3E7A5E",   # green
    "Site_3": "#B0682F",   # burnt orange
}
REF_GREY = "#9AA3AF"
LINK     = "#CBD1D9"
INK      = "#23262B"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 16,
    "axes.labelsize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 13,
    "legend.fontsize": 13,
    "axes.linewidth": 1.1,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#4a4a4a",
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": "#4a4a4a",
    "ytick.color": "#4a4a4a",
    "pdf.fonttype": 42,     # editable text in the vector output
    "ps.fonttype": 42,
})


def load(tsv):
    """Replicate the original analysis's preparation exactly."""
    d = pd.read_csv(tsv, sep="\t")
    d = (d.sort_values("circle_reads", ascending=False)
           .drop_duplicates(subset=["targetsite", "offtarget_seq"], keep="first"))
    d = d[(d.grch38_copy >= 1) & (d.grch38_copy <= 6)]
    return d


def marker_area(reads):
    """CIRCLE-seq read counts span 4 to 1880, so scale on sqrt to keep the
    largest marker readable without drowning the smallest."""
    r = np.clip(np.asarray(reads, dtype=float), None, 900.0)
    return 30 + 13.0 * np.sqrt(r)


def build(d, outdir):
    e = d[d["class"] == "individual_elevated"].copy()
    e["gain"] = e.k562_copy - e.grch38_copy
    e = e.sort_values("gain").reset_index(drop=True)
    y = np.arange(len(e))

    fig, ax = plt.subplots(figsize=(9.4, 8.2))

    # connector first, so the markers sit on top of it
    ax.hlines(y, e.grch38_copy, e.k562_copy, color=LINK, lw=2.0, zorder=1,
              capstyle="round")

    # reference copy number: one neutral colour, it is the baseline
    ax.scatter(e.grch38_copy, y, s=46, c=REF_GREY, edgecolors="white",
               linewidths=0.9, zorder=3)

    # K562 copy number: coloured by guide, sized by measured cleavage strength
    for guide, sub in e.groupby("targetsite"):
        idx = sub.index.to_numpy()
        ax.scatter(sub.k562_copy, idx,
                   s=marker_area(sub.circle_reads),
                   c=GUIDE_COLOURS.get(guide, "#666666"),
                   edgecolors="white", linewidths=0.9, zorder=4, alpha=0.95)

    # the headline site, annotated factually
    hi = e.circle_reads.idxmax()
    # NB: "most strongly cleaved" is scoped to the sites PLOTTED here. The HBB
    # guide has an off-target with 2,356 reads, but it is absent from GRCh38
    # (class reference_blind) and so is not in this panel.
    ax.annotate(
        f"{int(e.loc[hi, 'circle_reads']):,} CIRCLE-seq reads\n"
        f"(most strongly cleaved site shown;\n{e.loc[hi, 'targetsite']} guide)",
        xy=(e.loc[hi, "k562_copy"], hi),
        xytext=(5.55, hi - 7.2),
        fontsize=14.5, color=INK, ha="left", va="center", linespacing=1.3,
        arrowprops=dict(arrowstyle="-", lw=1.0, color="#6b6b6b",
                        shrinkA=0, shrinkB=11,
                        connectionstyle="angle3,angleA=0,angleB=70"),
        zorder=7)

    ax.set_xlim(0.55, 9.4)
    ax.set_ylim(-1.2, len(e) + 0.6)
    ax.set_yticks([])
    ax.set_xlabel("Copy number")
    ax.set_ylabel(f"Confirmed off-target sites  ($n$ = {len(e)})", labelpad=14)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(width=1.1, length=6)
    ax.grid(True, axis="x", ls=":", lw=0.8, color="#d5d5d5", zorder=0)
    ax.set_axisbelow(True)

    # ---- legends, both outside the data ----
    order = [g for g in ["HBB", "EMX1", "FANCF", "RNF2", "Site_1", "Site_2", "Site_3"]
             if g in set(e.targetsite)]
    guide_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=9,
               markerfacecolor=GUIDE_COLOURS[g], markeredgecolor="white",
               markeredgewidth=0.9, label=g.replace("_", " "))
        for g in order]
    guide_handles.insert(0, Line2D([], [], marker="o", linestyle="none",
                                   markersize=8, markerfacecolor=REF_GREY,
                                   markeredgecolor="white", markeredgewidth=0.9,
                                   label="GRCh38 (reference)"))
    leg1 = ax.legend(handles=guide_handles, loc="upper left",
                     bbox_to_anchor=(0.0, -0.105), ncol=4, frameon=False,
                     handletextpad=0.5, columnspacing=2.2,
                     title="K562 copy number, coloured by guide")
    leg1.get_title().set_fontsize(14)
    leg1.get_title().set_color("#4a4a4a")
    ax.add_artist(leg1)

    # size key: three round numbers spanning the observed range
    for reads, label in [(10, "10"), (200, "200"), (1500, "1,500")]:
        ax.scatter([], [], s=marker_area(reads), c="#b9b9b9",
                   edgecolors="white", linewidths=0.9, label=label)
    handles, labels = ax.get_legend_handles_labels()
    leg2 = ax.legend(handles[-3:], labels[-3:], loc="upper left",
                     bbox_to_anchor=(0.0, -0.275), frameon=False, ncol=3,
                     handletextpad=0.9, columnspacing=3.0,
                     title="marker area: CIRCLE-seq cleavage reads")
    leg2.get_title().set_fontsize(14)
    leg2.get_title().set_color("#4a4a4a")

    fig.tight_layout(pad=0.6)
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"circleseq_copynumber.{ext}")
        fig.savefig(path, dpi=400, bbox_inches="tight")
        print(f"  wrote {path}")
    plt.close(fig)

    # numbers worth having alongside the figure
    print(f"\n  sites shown           : {len(e)}")
    print(f"  guides represented    : {', '.join(order)}")
    print(f"  copy-number gain      : {e.gain.min():.2f} to {e.gain.max():.2f}")
    print(f"  fold change (K562/ref): {(e.k562_copy / e.grch38_copy).min():.2f}"
          f" to {(e.k562_copy / e.grch38_copy).max():.2f}")
    print(f"  mismatch distances    : {sorted(e.n_mm.unique())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", default=os.path.join(ROOT, "data",
                                                  "k562_offtarget_copynumber.tsv"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data"))
    a = ap.parse_args()
    if not os.path.exists(a.tsv):
        raise SystemExit(f"missing input: {a.tsv}")
    d = load(a.tsv)
    print(f"{len(d)} off-targets present in the reference (GRCh38 copy 1-6)")
    print(d["class"].value_counts().to_string(), "\n")
    build(d, a.outdir)


if __name__ == "__main__":
    main()
