#!/usr/bin/env python3
# headtohead_common.py   (run in ~/human_assemblies_kmer_count/july-10)
# COMMON-DENOMINATOR head-to-head.
# Fixed class-guide set per class = union of the 3 tools' class guides
# (~= all class guides present in the 50 windows). For each tool, over that
# SAME set, count accept / reject / not-designed.
import glob, csv
import pandas as pd

CLASSES=["multiple_v_single","ref_multiple_v_single","present_v_absent","ref_present_v_absent"]
TRUTH={"multiple_v_single":"reject","ref_multiple_v_single":"accept",
       "present_v_absent":"accept","ref_present_v_absent":"reject"}
K_CUT,MIT_CUT,GS_CUT=1.5,50,0.2
def canon(k): return min(k,k.translate(str.maketrans("ACGT","TGCA"))[::-1])

def krisp(cls):
    d={}
    for f in glob.glob(f"scores/{cls}/*.csv"):
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
    for f in glob.glob(f"crispor_class_scores/{cls}/*_crispor.tsv"):
        with open(f) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip(): continue
                p=ln.rstrip("\n").split("\t")
                if len(p)>=5 and len(p[2].strip())>=20:
                    try:
                        g=canon(p[2].strip().upper()); mit=float(p[3])   # mitSpecScore
                        d[g]="reject" if mit<MIT_CUT else "accept"
                    except: pass
    return d
def guidescan(cls):
    d={}
    with open(f"guidescan_class_scores/out/{cls}_ont.csv") as fh:
        for row in csv.DictReader(fh):
            g=canon(row["match_sequence"].upper()); sp=float(row["specificity"])
            if g not in d: d[g]="reject" if sp<GS_CUT else "accept"
    return d

print(f"reject rules: kRISPmer Ecuts>{K_CUT} | CRISPOR MITspec<{MIT_CUT} | GuideScan2 spec<{GS_CUT}")
print("denominator = FIXED class-guide set per class (union of the 3 tools)\n")
rows=[]
for cls in CLASSES:
    kd,cd,gd=krisp(cls),crispor(cls),guidescan(cls)
    cand=set(kd)|set(cd)|set(gd)
    classset=set()
    with open(f"{cls}.txt") as fh:                       # stream bucket (present_v_absent ~60M lines)
        for ln in fh:
            g=canon(ln.split()[0].upper())
            if g in cand: classset.add(g)
    N=len(classset)
    print(f"### {cls}   (truth: {TRUTH[cls]})   fixed N = {N}")
    print(f"{'tool':12s} {'accept':>13s} {'reject':>13s} {'not-designed':>15s}")
    for name,d in [("kRISP-meR",kd),("CRISPOR",cd),("GuideScan2",gd)]:
        a=sum(1 for g in classset if d.get(g)=="accept")
        r=sum(1 for g in classset if d.get(g)=="reject")
        nd=N-a-r
        print(f"{name:12s} {a:5d} ({100*a//N:3d}%) {r:5d} ({100*r//N:3d}%) {nd:5d} ({100*nd//N:3d}%)")
        rows.append({"class":cls,"truth":TRUTH[cls],"N":N,"tool":name,
                     "accept":a,"reject":r,"not_designed":nd})
    print()
pd.DataFrame(rows).to_csv("headtohead_common.tsv",sep="\t",index=False)
print("wrote headtohead_common.tsv")
