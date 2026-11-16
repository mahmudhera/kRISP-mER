#!/usr/bin/env python3
# 03_verify_targets_n22.py
# ----------------------------------------------------------------------------
# GATE: after building fresh targets, confirm each class's windows actually
# contain that class's OWN bucket guides (strong diagonal) BEFORE spending a
# kRISP-meR run. This is the exact matrix that exposed the earlier 0/50 break.
#
# Pass criterion: the diagonal (window contains its own class guide) is strongly
# non-zero AND dominates the OTHER class buckets. We EXCLUDE the present_v_absent
# column from "other" because every window is cut from N22 exons, so it is
# unavoidably full of incidental N22-present guides that fall into the huge
# present_v_absent bucket (~6.85M guides = the N22-only universe). That column is
# background, not target contamination; scoring (04) isolates the exact class
# guide per window via bucket intersection, so this background never affects the
# actual Ecuts numbers -- it only inflates this coarse matrix.
MIN_DIAG = 20          # need a real signal of the class guide in the windows
# ----------------------------------------------------------------------------
import os, glob
FRESH  = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
BD     = os.path.expanduser("~/krispmer/new_datasets/exps/n22")
CLS    = ["present_v_absent","ref_present_v_absent","multiple_v_single","ref_multiple_v_single"]
def canon(k): return min(k,k.translate(str.maketrans("ACGT","TGCA"))[::-1])

def bucket(n):
    with open(f"{BD}/{n}.txt") as fh:
        return {canon(l.split()[0].upper()) for l in fh}
B = {n:bucket(n) for n in CLS}

print("target-window guides  x  bucket   (diagonal should dominate each row)\n")
hdr = "targets \\ bucket".ljust(24) + "".join(c[:10].rjust(12) for c in CLS)
print(hdr); print("-"*len(hdr))
ok = True
for cls in CLS:
    tg=set()
    for f in glob.glob(f"{FRESH}/targets/{cls}/*.fasta"):
        seq="".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
        for i in range(len(seq)-22): tg.add(canon(seq[i:i+23]))
    row = {b: len(tg & B[b]) for b in CLS}
    diag = row[cls]
    # "other" = the other CLASS buckets, excluding the present_v_absent background
    others = [row[b] for b in CLS if b != cls and b != "present_v_absent"]
    max_other = max(others) if others else 0
    diag_ok = (diag >= MIN_DIAG and diag >= max_other)
    ok &= diag_ok
    mark = "  OK" if diag_ok else "  <-- FAIL"
    print(cls.ljust(24) + "".join(str(row[b]).rjust(12) for b in CLS) + mark)
print(f"\n(diagonal must be >= {MIN_DIAG} and >= the other CLASS buckets; "
      f"present_v_absent column is expected background and excluded)")
print("GATE:", "PASS — every class window contains its own class guide, safe to run kRISP-meR" if ok
      else "FAIL — a class's windows lack its own guides; do NOT run kRISP-meR")
