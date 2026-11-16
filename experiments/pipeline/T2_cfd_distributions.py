#!/usr/bin/env python3
# T2_cfd_distributions.py  ORG  ANCHOR
# ----------------------------------------------------------------------------
# TASK 2: per-class CFD-specificity distributions, split by kRISP-meR's verdict.
# For every class, two groups of guides:
#    ACCEPTED       = kRISP-meR generated & accepted (Ecuts <= 1.5)
#    REJECTED-side  = kRISP-meR (generated & rejected)  UNION  (not-generated)
# For each guide we take CRISPOR's cfdSpecScore (col cfdSpecScore, 0-100, high =
# specific/safe). CRISPOR designs ~all window guides, so it covers the
# kRISP-meR-not-generated guides too. Output a TSV of (class, group, guide,
# cfd_spec) for all 4 classes -> plotted separately.
#
# Usage on VM:
#   python3 T2_cfd_distributions.py rice n22      # rice, N22-anchored
#   python3 T2_cfd_distributions.py rice japonica # rice, japonica-anchored
#   python3 T2_cfd_distributions.py human chm13   # human, CHM13-anchored
#   python3 T2_cfd_distributions.py human grch38  # human, GRCh38-anchored (existing)
import os, sys, glob, csv
import pandas as pd

ORG    = sys.argv[1] if len(sys.argv) > 1 else "rice"
ANCHOR = sys.argv[2] if len(sys.argv) > 2 else "n22"
CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
           "ref_multiple_v_single":"accept","present_v_absent":"accept"}
K_CUT = 1.5
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))

# ---- locate the run dirs for (org, anchor) --------------------------------
if ORG == "rice":
    BASE = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
    BUCKET_DIR = os.path.expanduser("~/krispmer/new_datasets/exps/n22")
    if ANCHOR == "n22":       KDIR, CDIR, TDIR = "scores_bothstrand", "crispor_class_scores", "targets"
    else:                     KDIR, CDIR, TDIR = "ref_scores_bothstrand", "ref_crispor_class_scores", "ref_targets"
else:  # human
    BASE = os.path.expanduser("~/human_assemblies_kmer_count/july-10")
    BUCKET_DIR = BASE
    if ANCHOR == "chm13":     KDIR, CDIR, TDIR = "ind_scores_bothstrand", "ind_crispor_class_scores", "ind_targets"
    else:                     KDIR, CDIR, TDIR = "scores_bothstrand", "crispor_class_scores", "targets"

def class_guides_in_windows(cls):
    bucket=set(canon(l.split()[0].upper()) for l in open(f"{BUCKET_DIR}/{cls}.txt") if len(l.split()[0])==23)
    cg=set()
    for f in glob.glob(f"{BASE}/{TDIR}/{cls}/*.fasta")+glob.glob(f"{BASE}/{TDIR}/{cls}/*.fa"):
        seq="".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
        for i in range(len(seq)-22):
            c=canon(seq[i:i+23])
            if c in bucket: cg.add(c)
    return cg

def krisp_verdict(cls):
    d={}
    for f in glob.glob(f"{BASE}/{KDIR}/{cls}/*.csv"):
        try: df=pd.read_csv(f,comment="#",skip_blank_lines=True)
        except: continue
        if "tgt_in_plus" not in df.columns: continue
        for _,r in df.iterrows():
            try:
                g=canon(str(r["tgt_in_plus"]).upper()); ec=float(r["Expected_cuts_in_genome"])
                d[g]="reject" if ec>K_CUT else "accept"
            except: pass
    return d   # guides NOT here = not-generated

def crispor_cfd(cls):
    d={}
    for f in glob.glob(f"{BASE}/{CDIR}/{cls}/*_crispor.tsv"):
        with open(f) as fh:
            rd=csv.reader(fh,delimiter="\t"); hdr=None
            for row in rd:
                if row and row[0].startswith("#seqId"): hdr=row; continue
                if hdr is None or len(row)<5: continue
                try:
                    seq=row[2].strip().upper()
                    if len(seq)<20: continue
                    cfd=float(row[4])                 # cfdSpecScore
                    d[canon(seq)]=cfd
                except: pass
    return d

rows=[]
for cls in CLASSES:
    classset=class_guides_in_windows(cls)
    kd=krisp_verdict(cls); cf=crispor_cfd(cls)
    for g in classset:
        v=kd.get(g,"not_generated")
        group="accepted" if v=="accept" else "rejected_side"   # reject OR not_generated
        cfd=cf.get(g)
        if cfd is None: continue     # no CRISPOR CFD for this guide (rare)
        rows.append({"class":cls,"truth":TRUTH[cls],"group":group,
                     "krisp_verdict":v,"cfd_spec":cfd})

out=f"{BASE}/cfd_dist_{ORG}_{ANCHOR}.tsv"
pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
print(f"wrote {out}  ({len(rows)} guide-CFD rows)")
# quick summary
df=pd.DataFrame(rows)
for cls in CLASSES:
    sub=df[df["class"]==cls]
    for grp in ["accepted","rejected_side"]:
        s=sub[sub["group"]==grp]["cfd_spec"]
        if len(s):
            print(f"  {cls:22s} {grp:14s} n={len(s):4d}  cfd_spec median={s.median():5.1f} mean={s.mean():5.1f}")
        else:
            print(f"  {cls:22s} {grp:14s} n=   0")
