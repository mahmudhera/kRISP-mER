#!/usr/bin/env python3
# verify_notgen_readsupport.py  [class]   (default: present_v_absent)
# ----------------------------------------------------------------------------
# Verifies the "not-generated = no read support" mechanism (krispmer.py:291/742).
# For a class, split its window class-guides into GENERATED (in the scores CSVs)
# vs NOT-GENERATED, then check each guide two ways:
#   (A) present in the N22 ASSEMBLY  -> lexo_mer_counts_PAM_sorted_n22.txt
#   (B) present in the N22 READS     -> n22-real/mer_counts.jf  (jellyfish query)
# CLAIM: not-generated guides are ASSEMBLY-present but READS-absent (count 0),
#        while generated guides are READS-present (count > 0).
# Run on the VM (needs the jf + assembly dump). Prints a clean contingency table.
import glob, os, subprocess, tempfile, sys, statistics
import pandas as pd

cls = sys.argv[1] if len(sys.argv) > 1 else "present_v_absent"
FRESH = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
BD    = os.path.expanduser("~/krispmer/new_datasets/exps/n22")
ASM   = f"{BD}/lexo_mer_counts_PAM_sorted_n22.txt"        # N22 assembly PAM 23-mers
JF    = os.path.expanduser("~/krispmer/new_datasets/n22-real/mer_counts.jf")  # N22 reads

def rc(s):    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))

# ---- class guides present in this class's target windows -------------------
bucket = set(canon(l.split()[0].upper()) for l in open(f"{BD}/{cls}.txt"))
win = set()
for f in glob.glob(f"{FRESH}/targets/{cls}/*.fasta"):
    seq = "".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
    for i in range(len(seq)-22):
        c = canon(seq[i:i+23])
        if c in bucket: win.add(c)

# ---- which of them kRISP-meR actually wrote to output (generated) ----------
gen = set()
for f in glob.glob(f"{FRESH}/scores/{cls}/*.csv"):
    try: df = pd.read_csv(f, comment="#", skip_blank_lines=True)
    except: continue
    if "tgt_in_plus" not in df.columns: continue
    for g in df["tgt_in_plus"].dropna(): gen.add(canon(str(g).upper()))

generated  = win & gen
notgen     = win - gen

# ---- (A) assembly membership: stream the assembly dump once ----------------
def in_assembly(guides):
    want = set(guides); found = set()
    with open(ASM) as fh:
        for ln in fh:
            g = canon(ln.split()[0].upper())
            if g in want: found.add(g)
    return found

# ---- (B) read counts ------------------------------------------------------
# IMPORTANT: `jellyfish query` takes k-mers as POSITIONAL ARGS (mers:string+),
# NOT from stdin (unless -i). It also does NOT canonicalize the query, so we
# pass BOTH strands and take the max. Chunk the args to stay under ARG_MAX.
def read_counts(guides):
    guides = list(guides)
    mers = []
    for g in guides: mers += [g, rc(g)]
    seen = {}
    CHUNK = 2000
    for j in range(0, len(mers), CHUNK):
        batch = mers[j:j+CHUNK]
        out = subprocess.run(["jellyfish","query",JF, *batch],
                             capture_output=True, text=True).stdout
        for ln in out.split("\n"):
            p = ln.split()
            if len(p) == 2:
                c = canon(p[0].upper()); seen[c] = max(seen.get(c,0), int(p[1]))
    return {g: seen.get(g, 0) for g in guides}

print(f"class = {cls}")
print(f"class guides in windows : {len(win)}  (generated {len(generated)}, not-generated {len(notgen)})\n")

asm_all = in_assembly(win)
rc_all  = read_counts(win)

def report(name, S):
    n = len(S)
    if not n: print(f"  {name}: (empty)"); return
    in_asm  = sum(1 for g in S if g in asm_all)
    reads0  = sum(1 for g in S if rc_all[g] == 0)
    reads_pos = n - reads0
    vals = sorted(rc_all[g] for g in S)
    print(f"  {name:14s} n={n:3d} | in N22 assembly: {in_asm}/{n} ({100*in_asm//n}%) "
          f"| reads==0: {reads0}/{n} ({100*reads0//n}%) | reads>0: {reads_pos} "
          f"| read median={statistics.median(vals)} max={vals[-1]}")

print("membership check (assembly presence vs read presence):")
report("generated",     generated)
report("NOT-generated", notgen)
print("\nCLAIM: not-generated -> ~100% in assembly AND ~100% reads==0 ;"
      " generated -> reads>0")
