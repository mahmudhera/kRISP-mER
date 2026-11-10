#!/usr/bin/env python3
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np, os

OUT = "/mnt/windows-data/Thesis/plots"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#333333",
    "axes.labelcolor":"#222222","text.color":"#222222","xtick.color":"#333333","ytick.color":"#333333",
    "font.size":10,"font.family":"DejaVu Sans","axes.linewidth":1.0})

KR, CR, GS = "#9DC3CB", "#A6C08D", "#DFA9B2"
EDGE = "#4d4d4d"
COL = {"kRISP-meR":KR,"CRISPOR":CR,"GuideScan2":GS}
CLS   = ["multiple_v_single","ref_multiple_v_single","present_v_absent","ref_present_v_absent"]
SHORT = ["Multiple in individual\nSingle in reference",
         "Single in individual\nMultiple in reference",
         "Present in individual\nAbsent in reference",
         "Absent in individual\nPresent in reference"]
tools = ["kRISP-meR","CRISPOR","GuideScan2"]
AGREE = {"multiple_v_single":(96,49,41),"ref_multiple_v_single":(26,0,72),
         "present_v_absent":(62,78,None),"ref_present_v_absent":(100,79,37)}
TRUTH = {"multiple_v_single":"reject","ref_multiple_v_single":"accept",
         "present_v_absent":"accept","ref_present_v_absent":"reject"}

teal = LinearSegmentedColormap.from_list("teal",["#f4f8f9","#9DC3CB","#3f7d88","#245a63"])
orng = LinearSegmentedColormap.from_list("orng",["#fdf5ec","#f3c48c","#e0913f","#a85f18"])

# ================= FIG 1: agreement grouped bar =================
fig,ax=plt.subplots(figsize=(8.4,4.2)); x=np.arange(len(CLS)); w=0.26
for i,t in enumerate(tools):
    for j,c in enumerate(CLS):
        v=AGREE[c][i]; xi=x[j]+(i-1)*w
        if v is None:
            ax.text(xi,2,"no gRNA found",ha="center",va="bottom",fontsize=6.5,color="#9aa1ac",rotation=90); continue
        ax.bar(xi,v,w,color=COL[t],edgecolor=EDGE,linewidth=0.7,label=t if j==0 else None,zorder=3)
        ax.text(xi,v+1.2,f"{v}",ha="center",va="bottom",fontsize=7,color="#555")
ax.set_xticks(x); ax.set_xticklabels(SHORT,fontsize=8)
ax.set_ylim(0,100); ax.set_ylabel("% agreement with true genome",fontsize=9.5)
ax.set_yticks([0,25,50,75,100]); ax.grid(axis="y",color="#eef0f4",lw=0.8,zorder=0)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=3,frameon=False,
          fontsize=9,handlelength=1.1,columnspacing=1.6)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_agreement.png",dpi=220); plt.close()

# ================= FIG 2: correct-action heatmap (teal=rejected, orange=accepted) =================
img=np.ones((3,len(CLS),4))
for j,c in enumerate(CLS):
    cmap = teal if TRUTH[c]=="reject" else orng
    for i in range(3):
        v=AGREE[c][i]
        img[i,j]=(0.965,0.965,0.975,1) if v is None else cmap(v/100.0)
fig,ax=plt.subplots(figsize=(8.0,4.0))
ax.imshow(img,aspect="auto")
for i in range(3):
    for j,c in enumerate(CLS):
        v=AGREE[c][i]
        if v is None:
            ax.text(j,i,"no gRNA\nfound",ha="center",va="center",fontsize=7.5,color="#999"); continue
        ax.text(j,i,f"{int(v)}",ha="center",va="center",fontsize=10,
                color="white" if v>58 else "#222")
HM=[SHORT[j]+f"\n(should {TRUTH[c]})" for j,c in enumerate(CLS)]
ax.set_xticks(range(len(CLS))); ax.set_xticklabels(HM,fontsize=7)
ax.set_yticks(range(3)); ax.set_yticklabels(tools,fontsize=9.5)
ax.set_xticks(np.arange(-.5,len(CLS),1),minor=True); ax.set_yticks(np.arange(-.5,3,1),minor=True)
ax.grid(which="minor",color="white",lw=2.5); ax.tick_params(which="minor",length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.legend(handles=[Patch(facecolor=teal(0.78),label="% rejected"),
                   Patch(facecolor=orng(0.72),label="% accepted")],
          loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2,frameon=False,
          fontsize=9,handlelength=1.2,columnspacing=2.0)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_heatmap.png",dpi=220); plt.close()

# ================= FIG 3: kRISP-meR median Expected_cuts =================
MED=[("Present in individual\nAbsent in reference",1.09),
     ("Single in individual\nMultiple in reference",2.12),
     ("Multiple in individual\nSingle in reference",2.95)]
fig,ax=plt.subplots(figsize=(6.2,3.8)); xs=range(len(MED))
ax.bar(xs,[m for _,m in MED],0.6,color=KR,edgecolor=EDGE,linewidth=0.7,zorder=3)
for i,(_,m) in enumerate(MED): ax.text(i,m+0.05,f"{m:.2f}",ha="center",va="bottom",fontsize=8,color="#555")
ax.set_xticks(xs); ax.set_xticklabels([n for n,_ in MED],fontsize=7.5)
ax.set_ylabel("median Expected_cuts",fontsize=9.5); ax.set_ylim(0,3.4)
ax.grid(axis="y",color="#eef0f4",lw=0.8,zorder=0)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_ecuts.png",dpi=220); plt.close()
print("done")
