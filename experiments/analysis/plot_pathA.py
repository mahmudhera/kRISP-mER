#!/usr/bin/env python3
# plot_pathA.py  --  Path A figure (soft, minimal). Run in ~/pathA_circleseq.
# Off-targets present in the reference (grch38 copy 1-6): concordant on the
# diagonal, individual-elevated above it (amplified in K562).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.facecolor":"#ffffff","axes.facecolor":"#fbfbfd","axes.edgecolor":"#d7dae0",
    "axes.linewidth":1.0,"axes.labelcolor":"#3a3f4b","text.color":"#3a3f4b",
    "xtick.color":"#6b7280","ytick.color":"#6b7280","font.size":11,
    "axes.grid":True,"grid.color":"#eceef3","grid.linewidth":1.0,
})

d = pd.read_csv("k562_offtarget_copynumber.tsv", sep="\t")
d = d.sort_values("circle_reads", ascending=False)\
     .drop_duplicates(subset=["targetsite","offtarget_seq"], keep="first")
d = d[(d.grch38_copy >= 1) & (d.grch38_copy <= 6)].copy()

YMAX, XMAX = 12.0, 6.0
d["yv"] = d.k562_copy.clip(upper=YMAX)
rng = np.random.default_rng(7)
d["xv"] = d.grch38_copy + rng.uniform(-0.13, 0.13, len(d))
d["ms"] = 22 + d.circle_reads.clip(upper=350)/2.4

palette = {"concordant":"#c3ccd9", "individual_elevated":"#ec7a8f"}
labels  = {"concordant":"concordant", "individual_elevated":"individual-elevated"}
z = {"concordant":1, "individual_elevated":3}

fig, ax = plt.subplots(figsize=(8.2, 6.6))
ax.fill_between([0,YMAX],[0,YMAX],[YMAX,YMAX], color="#ec7a8f", alpha=0.05, zorder=0)
ax.plot([0,YMAX],[0,YMAX], ls=(0,(6,5)), lw=1.3, color="#aab2bf", zorder=1)

for cls in ["concordant","individual_elevated"]:
    s = d[d["class"]==cls]
    ax.scatter(s.xv, s.yv, s=s.ms, c=palette[cls], alpha=0.72,
               edgecolors="white", linewidths=0.7, zorder=z[cls],
               label=f"{labels[cls]}  (n={len(s)})")

ax.set_xlim(0.4, XMAX+0.45); ax.set_ylim(-0.45, YMAX+0.45)
ax.set_xlabel("GRCh38 copy number  (reference)", fontsize=12)
ax.set_ylabel("K562 copy number  (individual, reads / lambda)", fontsize=12)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.set_axisbelow(True)
leg = ax.legend(loc="lower right", frameon=True, fontsize=10, borderpad=0.9, labelspacing=0.7)
leg.get_frame().set_facecolor("#ffffff"); leg.get_frame().set_edgecolor("#e2e5ea")

fig.tight_layout()
fig.savefig("k562_vs_grch38_copynumber.png", dpi=220, bbox_inches="tight")
print("wrote k562_vs_grch38_copynumber.png")
print("\nshown:", len(d), "off-targets present in reference (grch38 copy 1-6)")
print(d["class"].value_counts())
