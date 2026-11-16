#!/usr/bin/env python3
# T1_enrichment_human.py  [ANCHOR=grch38]
# ----------------------------------------------------------------------------
# TASK 1 rigor: enrichment of FAMOUS/important genes among the guides that
# reference tools MIS-ACCEPT (dangerous) vs the guides kRISP-meR ACCEPTS (safe).
#
# Two guide groups (per class and pooled):
#   MISACC = reference tool (CRISPOR or GuideScan2) accepts, assembly-CFD >= DUP_CUT
#            (duplicated / high off-target in the individual), kRISP-meR does NOT accept.
#   KSAFE  = kRISP-meR accepts (Ecuts<=1.5).
# For each group: fraction of guides whose gene is in the famous-gene universe.
# Fisher's exact 2x2 (in-famous vs not) tests whether MISACC is enriched vs KSAFE.
#
# guide->ENST->symbol via the same grch38 SAM + enst2sym map as T1_famous_genes_human.
# Env: DATA, SAM_DIR, ASM_DIR (as in T1_famous_genes_human.py).  Run ON VM.
import os, sys, csv, re
from collections import defaultdict

ANCHOR = sys.argv[1] if len(sys.argv) > 1 else "grch38"
DATA    = os.path.expanduser(os.environ.get("DATA", "/mnt/windows-data/Thesis/data"))
SAM_DIR = os.path.expanduser(os.environ.get("SAM_DIR", f"{DATA}/human_sams"))
ASM_DIR = os.path.expanduser(os.environ.get("ASM_DIR", DATA))
CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
DUP_CUT = 2.0
REF_TOOLS = ["CRISPOR","GuideScan2"]

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))
def _find(name):
    for p in (f"{DATA}/{name}", f"{DATA}/famous_genes/{name}"):
        if os.path.exists(p): return p
    sys.exit(f"[ERR] cannot find {name}")

famous=set()
for ln in open(_find("famous_genes_human.tsv")):
    p=ln.rstrip("\n").split("\t")
    if len(p)>=1 and p[0] not in ("symbol",""): famous.add(p[0].upper())
enst2sym={}
for ln in open(_find("enst2sym.tsv")):
    p=ln.rstrip("\n").split("\t")
    if len(p)==2: enst2sym[p[0].split(".")[0]]=p[1]
ENST_RE=re.compile(r"(ENST\d+)")
print(f"famous genes: {len(famous)} ; ENST->symbol: {len(enst2sym)}")

def load_asm(anchor):
    p=f"{ASM_DIR}/asm_cfd_dist_human_{anchor}.tsv"
    G=defaultdict(lambda:{"cfd":None,"verd":{},"seq":None,"cls":None})
    for r in csv.DictReader(open(p),delimiter="\t"):
        g=canon(r["guide"].upper()); key=(r["class"],g)
        G[key]["cfd"]=float(r["cfd_cuts"]); G[key]["seq"]=r["guide"].upper()
        G[key]["cls"]=r["class"]; G[key]["verd"][r["tool"]]=r["verdict"]
    return G

def groups(anchor):
    """per class -> (misacc_set, ksafe_set) of canon guides."""
    G=load_asm(anchor)
    out={c:(set(),set()) for c in CLASSES}
    for (cls,g),d in G.items():
        if cls not in out: continue
        kacc = d["verd"].get("kRISP-meR")=="accept"
        dangerous = d["cfd"] is not None and d["cfd"]>=DUP_CUT
        refacc = any(d["verd"].get(t)=="accept" for t in REF_TOOLS)
        if kacc: out[cls][1].add(g)                       # KSAFE
        if dangerous and refacc and not kacc: out[cls][0].add(g)  # MISACC
    return out

def guides_in_famous(cls, wanted):
    """stream the SAM once; return subset of `wanted` whose gene is famous."""
    sam=f"{SAM_DIR}/{cls}_exons_max_0.sam"
    if not os.path.exists(sam): print(f"  [WARN] no SAM {sam}"); return set(),{}
    g2fam=set(); g2sym=defaultdict(set)
    for ln in open(sam):
        if ln.startswith("@"): continue
        p=ln.split("\t")
        if len(p)<10: continue
        s=p[9].strip().upper()
        if len(s)!=23: continue
        c=canon(s)
        if c not in wanted: continue
        m=ENST_RE.search(p[2])
        if not m: continue
        sym=enst2sym.get(m.group(1))
        if sym:
            g2sym[c].add(sym)
            if sym in famous: g2fam.add(c)
    return g2fam, g2sym

def fisher(a,b,c,d):
    """Fisher exact right-tail p for 2x2 [[a,b],[c,d]] (a=misacc&fam, b=misacc&notfam,
       c=ksafe&fam, d=ksafe&notfam). Uses scipy if present else a hypergeom fallback."""
    try:
        from scipy.stats import fisher_exact
        _,p=fisher_exact([[a,b],[c,d]],alternative="greater"); return p
    except Exception:
        from math import comb
        n=a+b+c+d; row1=a+b; col1=a+c
        # right tail: sum P(X>=a)
        lo=max(0,col1-(n-row1)); hi=min(row1,col1)
        def P(x): return comb(row1,x)*comb(n-row1,col1-x)/comb(n,col1)
        return sum(P(x) for x in range(a,hi+1))

GR=groups(ANCHOR)
print(f"\n{'class':22s} {'MISACC in-fam/N':>16s} {'KSAFE in-fam/N':>16s} {'MISACC%':>8s} {'KSAFE%':>8s} {'ratio':>6s} {'p(Fisher)':>10s}")
tot=[0,0,0,0]
for cls in CLASSES:
    mis,ksafe=GR[cls]
    fam_mis,_=guides_in_famous(cls, mis)
    fam_saf,_=guides_in_famous(cls, ksafe)
    a=len(fam_mis); b=len(mis)-a; c=len(fam_saf); d=len(ksafe)-c
    tot[0]+=a; tot[1]+=b; tot[2]+=c; tot[3]+=d
    pm = 100*a/len(mis) if mis else 0
    ps = 100*c/len(ksafe) if ksafe else 0
    ratio = (pm/ps) if ps>0 else float('inf')
    p=fisher(a,b,c,d) if (mis and ksafe) else float('nan')
    print(f"{cls:22s} {a:6d}/{len(mis):<9d} {c:6d}/{len(ksafe):<9d} {pm:7.1f}% {ps:7.1f}% {ratio:6.2f} {p:10.4g}")
a,b,c,d=tot
pm=100*a/(a+b) if (a+b) else 0; ps=100*c/(c+d) if (c+d) else 0
ratio=(pm/ps) if ps>0 else float('inf'); p=fisher(a,b,c,d)
print(f"{'POOLED':22s} {a:6d}/{a+b:<9d} {c:6d}/{c+d:<9d} {pm:7.1f}% {ps:7.1f}% {ratio:6.2f} {p:10.4g}")
print("\nMISACC = ref tool accepts, asm-CFD>=2 (dangerous), kRISP rejects.  KSAFE = kRISP accepts.")
print("ratio>1 & small p => reference tools' mis-accepted guides are ENRICHED for famous/important genes.")
