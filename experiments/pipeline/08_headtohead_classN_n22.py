#!/usr/bin/env python3
# 08_headtohead_classN_n22.py — FRESH N22 3-tool head-to-head, CLASS-GUIDE denominator.
#
# CORRECTED denominator: N = ALL class-defining guides present in that class's 50
# target windows (canonical guides in the window that are members of the class
# bucket), INDEPENDENT of which tools designed them. A guide a tool did not design
# counts as "not-generated" against this fixed N. (The earlier version used the
# union of the three tools' designed guides as N, which let each tool's design
# coverage distort the denominator -- wrong.)
#
# For each tool over that fixed class-guide set: accept / reject / not-generated.
import glob, csv, os
import pandas as pd

FRESH   = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
BUCKET_DIR = os.path.expanduser("~/krispmer/new_datasets/exps/n22")
CLASSES = ["multiple_v_single","ref_multiple_v_single","present_v_absent","ref_present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_multiple_v_single":"accept",
           "present_v_absent":"accept","ref_present_v_absent":"reject"}
K_CUT,MIT_CUT,GS_CUT = 1.5, 50, 0.2
def canon(k): return min(k,k.translate(str.maketrans("ACGT","TGCA"))[::-1])

# ---- the fixed class-guide denominator: class guides present in the windows ----
def class_guides_in_windows(cls):
    bucket=set()
    with open(f"{BUCKET_DIR}/{cls}.txt") as fh:
        for ln in fh: bucket.add(canon(ln.split()[0].upper()))
    cg=set()
    for f in glob.glob(f"{FRESH}/targets/{cls}/*.fasta"):
        seq="".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
        for i in range(len(seq)-22):
            g=canon(seq[i:i+23])
            if g in bucket: cg.add(g)
    return cg

# ---- per-tool verdict maps (guide -> accept/reject) ----
def krisp(cls):
    d={}
    for f in glob.glob(f"{FRESH}/scores/{cls}/*.csv"):
        try: df=pd.read_csv(f,comment="#",skip_blank_lines=True)
        except: continue
        if "tgt_in_plus" not in df.columns: continue
        for _,r in df.iterrows():
            try:
                g=canon(str(r["tgt_in_plus"]).upper()); ec=float(r["Expected_cuts_in_genome"])
                d[g]="reject" if ec>K_CUT else "accept"
            except: pass
    return d

def crispor(cls):
    d={}
    for f in glob.glob(f"{FRESH}/crispor_class_scores/{cls}/*_crispor.tsv"):
        with open(f) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip(): continue
                p=ln.rstrip("\n").split("\t")
                if len(p)>=5 and len(p[2].strip())>=20:
                    try:
                        g=canon(p[2].strip().upper()); mit=float(p[3])
                        d[g]="reject" if mit<MIT_CUT else "accept"
                    except: pass
    return d

def guidescan(cls):
    d={}
    fn=f"{FRESH}/guidescan_class_scores/out/{cls}_ont.csv"
    if not os.path.exists(fn): return d
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if row["match_distance"].strip()!="0": continue
            g=canon(row["match_sequence"].upper()); sp=float(row["specificity"])
            if g not in d: d[g]="reject" if sp<GS_CUT else "accept"
    return d

print(f"reject rules: kRISPmer Ecuts>{K_CUT} | CRISPOR MITspec<{MIT_CUT} | GuideScan2 spec<{GS_CUT}")
print("denominator = ALL class guides present in the class's target windows (fixed)\n")
rows=[]
for cls in CLASSES:
    classset=class_guides_in_windows(cls)          # <-- fixed class-guide N
    kd,cd,gd=krisp(cls),crispor(cls),guidescan(cls)
    N=len(classset)
    print(f"### {cls}   (truth: {TRUTH[cls]})   class-guide N = {N}")
    print(f"{'tool':12s} {'accept':>13s} {'reject':>13s} {'not-designed':>15s}")
    for name,d in [("kRISP-meR",kd),("CRISPOR",cd),("GuideScan2",gd)]:
        a=sum(1 for g in classset if d.get(g)=="accept")
        r=sum(1 for g in classset if d.get(g)=="reject")
        nd=N-a-r
        pa,pr,pnd=(100*a//N,100*r//N,100*nd//N) if N else (0,0,0)
        print(f"{name:12s} {a:5d} ({pa:3d}%) {r:5d} ({pr:3d}%) {nd:5d} ({pnd:3d}%)")
        rows.append({"class":cls,"truth":TRUTH[cls],"N":N,"tool":name,
                     "accept":a,"reject":r,"not_designed":nd})
    print()
pd.DataFrame(rows).to_csv(f"{FRESH}/headtohead_classN_n22.tsv",sep="\t",index=False)
print(f"wrote {FRESH}/headtohead_classN_n22.tsv")
