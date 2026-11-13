#!/usr/bin/env python3
# headtohead_crispor_krispmer.py  (run in ~/human_assemblies_kmer_count/july-10)
# Per-class head-to-head: kRISP-meR Expected_cuts (CHM13 reads) vs CRISPOR
# cfdSpecScore + offtargetCount (hg38), matched per class guide.
import glob, os, statistics
import pandas as pd

CLASSES = ["multiple_v_single","present_v_absent","ref_multiple_v_single","ref_present_v_absent"]
KDIR, CDIR = "scores", "crispor_class_scores"

def canon(k):
    return min(k, k.translate(str.maketrans("ACGT","TGCA"))[::-1])

def load_krisp(f):
    try: df = pd.read_csv(f, comment="#", skip_blank_lines=True)
    except: return {}
    if "tgt_in_plus" not in df.columns: return {}
    d={}
    for _,r in df.iterrows():
        try: d[canon(str(r["tgt_in_plus"]).upper())] = float(r["Expected_cuts_in_genome"])
        except: pass
    return d

def load_crispor(f):
    d={}
    try:
        with open(f) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip(): continue
                p = ln.rstrip("\n").split("\t")
                if len(p) < 6: continue
                seq = p[2].strip().upper()
                if len(seq) < 20: continue
                try: d[canon(seq)] = (float(p[4]), int(p[5]))   # cfdSpec, offtargetCount
                except: continue
    except: pass
    return d

# pass 1: pair guides present in BOTH tools' output for the same target
paired = {c:[] for c in CLASSES}; allg = {c:set() for c in CLASSES}
for c in CLASSES:
    for kf in sorted(glob.glob(f"{KDIR}/{c}/*.csv")):
        name = os.path.basename(kf)[:-4]
        cf = f"{CDIR}/{c}/{name}_crispor.tsv"
        if not os.path.exists(cf): continue
        kd, cd = load_krisp(kf), load_crispor(cf)
        for g in set(kd) & set(cd):
            paired[c].append((g, kd[g], cd[g][0], cd[g][1])); allg[c].add(g)

# pass 2: keep only actual class guides (stream the bucket)
isclass={}
for c in CLASSES:
    want=allg[c]; found=set()
    with open(f"{c}.txt") as fh:
        for ln in fh:
            cg=canon(ln.split()[0].upper())
            if cg in want: found.add(cg)
    isclass[c]=found

# report
print(f"{'class':24s} {'n':>4s} {'kRISP Ecuts':>12s} {'CRISPOR cfdSpec':>16s} {'CRISPOR offtgt':>15s}")
print("-"*76)
rows=[]
for c in CLASSES:
    sel=[(g,ec,cfd,ot) for (g,ec,cfd,ot) in paired[c] if g in isclass[c]]
    if not sel:
        print(f"{c:24s} {0:4d}   (no paired class guides)"); continue
    ecs=[x[1] for x in sel]; cfds=[x[2] for x in sel]; ots=[x[3] for x in sel]
    print(f"{c:24s} {len(sel):4d} {statistics.median(ecs):12.2f} "
          f"{statistics.median(cfds):16.1f} {statistics.median(ots):15.1f}")
    for g,ec,cfd,ot in sel:
        rows.append({"class":c,"guide":g,"krispmer_Ecuts":ec,
                     "crispor_cfdSpec":cfd,"crispor_offtargets":ot})
pd.DataFrame(rows).to_csv("headtohead_class_guides.tsv", sep="\t", index=False)
print("\nwrote headtohead_class_guides.tsv  (per-guide, for plotting)")
print("\nreminder: CRISPOR cfdSpec HIGH = 'safe' (inverse of Expected_cuts).")
