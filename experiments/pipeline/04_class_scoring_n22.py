#!/usr/bin/env python3
# 04_class_scoring_n22.py  — run after kRISP-meR over the FRESH targets
# ----------------------------------------------------------------------------
# Per-class kRISP-meR scoring for the clean N22 rebuild.
# kRISP-meR scores EVERY NGG guide it finds in each 150bp window (one CSV row
# each). We isolate the CLASS guides by intersecting the scored guides with the
# class bucket (canonical match), then report class/total, median & mean
# Expected_cuts_in_genome. Same method as the human Table 3.
# ----------------------------------------------------------------------------
import glob, os, statistics, pandas as pd
FRESH   = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
BD      = os.path.expanduser("~/krispmer/new_datasets/exps/n22")
CLASSES = ["present_v_absent","ref_multiple_v_single","multiple_v_single","ref_present_v_absent"]
def canon(k): return min(k,k.translate(str.maketrans("ACGT","TGCA"))[::-1])

print(f"{'class':24s} {'class/total':>12s} {'median':>8s} {'mean':>8s}")
for cls in CLASSES:
    ge = {}                                   # canonical guide -> Expected_cuts
    for f in glob.glob(f"{FRESH}/scores/{cls}/*.csv"):
        try: df = pd.read_csv(f, comment="#", skip_blank_lines=True)
        except: continue
        if "tgt_in_plus" not in df.columns: continue
        for _, r in df.iterrows():
            try: ge[canon(str(r["tgt_in_plus"]).upper())] = float(r["Expected_cuts_in_genome"])
            except: pass
    cand = set(ge)
    classset = set()
    with open(f"{BD}/{cls}.txt") as fh:
        for ln in fh:
            g = canon(ln.split()[0].upper())
            if g in cand: classset.add(g)
    ecs = [ge[g] for g in classset]
    if ecs:
        print(f"{cls:24s} {f'{len(ecs)}/{len(cand)}':>12s} "
              f"{statistics.median(ecs):8.2f} {sum(ecs)/len(ecs):8.2f}")
    else:
        print(f"{cls:24s} {f'0/{len(cand)}':>12s} {'--':>8s} {'--':>8s}")
