#!/usr/bin/env python3
# panel_legends.py -- two standalone legend PDFs for the discordance figure.
#   legend_hh.pdf  : generated & accepted / generated & rejected / not generated
#   legend_cfd.pdf : single-copy (=1) / duplicated (=2) reference lines
# Figure-2 fonts, soothing palette. Place once in the LaTeX table.
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

OUT=os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans","Arial"],
    "font.size":18,"text.color":"#2b2b2b"})
ACC="#6BA292"; REJ="#E4A66B"; NOG="#D4D2CC"

# ---- head-to-head legend ----
h1=[Patch(facecolor=ACC,edgecolor="white",label="generated & accepted"),
    Patch(facecolor=REJ,edgecolor="white",label="generated & rejected"),
    Patch(facecolor=NOG,edgecolor="white",label="not generated")]
fig=plt.figure(figsize=(9.5,0.4))
fig.legend(handles=h1,loc="center",ncol=3,frameon=False,fontsize=17,handlelength=1.4,columnspacing=2.0)
for ext in ("pdf","png"):
    p=f"{OUT}/legend_hh.{ext}"; fig.savefig(p,dpi=300,bbox_inches="tight",pad_inches=0.02); print(f"wrote {p}")
plt.close()

# ---- CFD legend ----
h2=[Line2D([0],[0],color="#3a8f7a",ls=":",lw=2.2,label="single-copy (expected cuts = 1)"),
    Line2D([0],[0],color="#888",ls="--",lw=2.2,label="duplicated (expected cuts = 2)"),
    Line2D([0],[0],marker="D",color="none",markerfacecolor="white",markeredgecolor="#4a4a4a",
           markersize=11,markeredgewidth=1.6,label="mean")]
fig=plt.figure(figsize=(9.5,0.4))
fig.legend(handles=h2,loc="center",ncol=3,frameon=False,fontsize=17,handlelength=2.2,columnspacing=2.0)
for ext in ("pdf","png"):
    p=f"{OUT}/legend_cfd.{ext}"; fig.savefig(p,dpi=300,bbox_inches="tight",pad_inches=0.02); print(f"wrote {p}")
plt.close()
