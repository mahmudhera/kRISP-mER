#!/usr/bin/env python3
# 04b_class_scoring_human_bothstrand.py — Table 3 (human) on the BOTH-STRAND run.
# Reads scores_bothstrand/ (kRISP-meR run WITH -n). Reports class/total scored,
# median & mean Expected_cuts_in_genome per class.
import glob, os, statistics, pandas as pd
H       = os.path.expanduser("~/human_assemblies_kmer_count/july-10")
BD      = H
CLASSES = ["present_v_absent","ref_multiple_v_single","multiple_v_single","ref_present_v_absent"]
def canon(k): return min(k,k.translate(str.maketrans("ACGT","TGCA"))[::-1])
print(f"{'class':24s} {'class/total':>12s} {'median':>8s} {'mean':>8s}")
for cls in CLASSES:
    ge = {}
    for f in glob.glob(f"{H}/ind_scores_bothstrand/{cls}/*.csv"):
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
