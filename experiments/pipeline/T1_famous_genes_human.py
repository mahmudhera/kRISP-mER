#!/usr/bin/env python3
# T1_famous_genes_human.py  [ANCHOR=grch38]
# ----------------------------------------------------------------------------
# TASK 1 (human): show that reference tools recommend DANGEROUS guides inside
# FAMOUS (COSMIC cancer-census) genes.
#
# A "hit" = a COSMIC gene that contains a class guide which a REFERENCE TOOL
# ACCEPTS but is actually dangerous:  assembly-CFD >= DUP_CUT (duplicated in the
# individual) AND kRISP-meR does NOT accept it.  (definition chosen with Navid:
# "gene contains a mis-accepted guide", tied to Task 2c.)
#
# Chain, per class:
#   mis-accepted guide (from asm_cfd_dist_human_<anchor>.tsv + verdict)
#     -> ENST         (grch38 exon SAM: col10 seq -> regex ENST out of col3 rname)
#     -> gene symbol  (enst2sym.tsv, Ensembl GRCh38.110)
#     -> in COSMIC?   (cosmic_cgc_symbols.txt, 723 genes)
#
# Runs LOCALLY (needs the asm_cfd TSVs in data/, the SAMs on the VM or copied,
# enst2sym.tsv + cosmic list in data/). SAM is streamed once, filtered to the
# mis-accepted guide set only. Writes data/task1_human_<anchor>.tsv + a summary.
import os, sys, csv, re
from collections import defaultdict

ANCHOR = sys.argv[1] if len(sys.argv) > 1 else "grch38"
# Paths resolve for BOTH local and VM runs via env overrides:
#   DATA   = dir holding cosmic_cgc_symbols.txt, enst2sym.tsv, asm_cfd_dist_human_*.tsv
#   SAM_DIR= dir holding <class>_exons_max_0.sam
# On the VM run e.g.:
#   DATA=~/task1_human SAM_DIR=~/human_assemblies_kmer_count/july-10 \
#     ASM_DIR=~/human_assemblies_kmer_count/july-10 python3 T1_famous_genes_human.py grch38
DATA    = os.path.expanduser(os.environ.get("DATA", "/mnt/windows-data/Thesis/data"))
SAM_DIR = os.path.expanduser(os.environ.get("SAM_DIR", f"{DATA}/human_sams"))
ASM_DIR = os.path.expanduser(os.environ.get("ASM_DIR", DATA))   # where asm_cfd_dist_*.tsv live

CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
           "ref_multiple_v_single":"accept","present_v_absent":"accept"}
DUP_CUT = 2.0          # assembly-CFD >= 2 => duplicated in individual (dangerous)
REF_TOOLS = ["CRISPOR","GuideScan2"]

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))

# ---- famous-gene universe (labeled) + ENST->symbol map --------------------
# reference files may sit directly under DATA (VM: ~/task1_human/) or in a
# famous_genes/ subdir (local: Thesis/data/famous_genes/). Try both.
def _find(name):
    for p in (f"{DATA}/{name}", f"{DATA}/famous_genes/{name}"):
        if os.path.exists(p): return p
    sys.exit(f"[ERR] cannot find {name} under {DATA} or {DATA}/famous_genes")
# famous_genes_human.tsv: symbol<TAB>sources (COSMIC/ACMG/ClinVar/ClinGen/FDA/BROCA)
famous = {}   # SYMBOL -> "src1;src2"
for ln in open(_find("famous_genes_human.tsv")):
    p=ln.rstrip("\n").split("\t")
    if len(p)==2 and p[0]!="symbol": famous[p[0].upper()] = p[1]
enst2sym = {}
for ln in open(_find("enst2sym.tsv")):
    p=ln.rstrip("\n").split("\t")
    if len(p)==2: enst2sym[p[0].split(".")[0]] = p[1]   # strip version
print(f"famous genes: {len(famous)} ; ENST->symbol: {len(enst2sym)}")

ENST_RE = re.compile(r"(ENST\d+)")

# ---- per-class assembly-CFD + verdicts ------------------------------------
def load_asm(anchor):
    """canon guide -> {'cfd':x, 'tool_verdict':{tool:accept/reject/not_generated}, 'seq':w, 'class':c}"""
    p=f"{ASM_DIR}/asm_cfd_dist_human_{anchor}.tsv"
    G=defaultdict(lambda: {"cfd":None,"verd":{}, "seq":None})
    for r in csv.DictReader(open(p),delimiter="\t"):
        g=canon(r["guide"].upper()); key=(r["class"],g)
        G[key]["cfd"]=float(r["cfd_cuts"]); G[key]["seq"]=r["guide"].upper()
        G[key]["verd"][r["tool"]]=r["verdict"]
    return G

def mis_accepted(anchor):
    """per class: set of guides a reference tool ACCEPTS but are dangerous
       (assembly-CFD >= DUP_CUT) AND kRISP-meR does not accept."""
    G=load_asm(anchor)
    out={c:{} for c in CLASSES}   # class -> {canon_guide: {seq,cfd,tools:[...]}}
    for (cls,g),d in G.items():
        if cls not in out: continue
        if d["cfd"] is None or d["cfd"]<DUP_CUT: continue          # not dangerous
        if d["verd"].get("kRISP-meR")=="accept": continue          # kRISP accepts -> not a kRISP miss
        tools=[t for t in REF_TOOLS if d["verd"].get(t)=="accept"] # ref tools that accepted it
        if not tools: continue
        out[cls][g]={"seq":d["seq"],"cfd":d["cfd"],"tools":tools}
    return out

# ---- stream the grch38 exon SAM: guide seq -> set of ENST -----------------
def guide_to_enst(cls, wanted_seqs):
    """wanted_seqs = set of canon guide seqs; return {canon_seq: set(ENST)} by
       streaming <class>_exons_max_0.sam (grch38 arm)."""
    sam=f"{SAM_DIR}/{cls}_exons_max_0.sam"
    if not os.path.exists(sam):
        print(f"  [WARN] SAM not found: {sam}"); return {}
    hit=defaultdict(set)
    with open(sam) as fh:
        for ln in fh:
            if ln.startswith("@"): continue
            p=ln.split("\t")
            if len(p)<10: continue
            seq=p[9].strip().upper()
            if len(seq)!=23: continue
            c=canon(seq)
            if c not in wanted_seqs: continue
            m=ENST_RE.search(p[2])
            if m: hit[c].add(m.group(1))
    return hit

# ---- run -------------------------------------------------------------------
MA=mis_accepted(ANCHOR)
rows=[]; summary=[]
for cls in CLASSES:
    guides=MA[cls]
    wanted=set(guides.keys())
    print(f"### {cls}: {len(wanted)} mis-accepted dangerous guides -> streaming SAM ...", flush=True)
    g2e=guide_to_enst(cls, wanted) if wanted else {}
    famous_hits=set(); ghits=0
    for g,info in guides.items():
        ensts=g2e.get(g,set())
        syms={enst2sym.get(e) for e in ensts if enst2sym.get(e)}
        fam ={s for s in syms if s in famous}       # symbol in the famous-gene universe
        if fam:
            ghits+=1; famous_hits|=fam
            for s in sorted(fam):
                rows.append({"class":cls,"truth":TRUTH[cls],"famous_gene":s,
                             "famous_sources":famous[s],
                             "guide":info["seq"],"asm_cfd":round(info["cfd"],3),
                             "mis_accepted_by":";".join(info["tools"]),
                             "all_symbols":";".join(sorted(syms))})
    summary.append((cls,len(wanted),ghits,len(famous_hits),sorted(famous_hits)))
    print(f"    {ghits} guides land in {len(famous_hits)} FAMOUS genes: {sorted(famous_hits)[:14]}"
          + (" ..." if len(famous_hits)>14 else ""))

out=f"{DATA}/task1_human_{ANCHOR}.tsv"
import pandas as pd
pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
print(f"\nwrote {out}  ({len(rows)} guide-in-famous-gene rows)")
print("\n=== SUMMARY: FAMOUS/important genes with a reference-tool-mis-accepted dangerous guide ===")
print(f"{'class':22s} {'truth':7s} {'mis-acc guides':>14s} {'->in famous':>11s} {'#famous genes':>13s}")
for cls,nma,gh,nfam,genes in summary:
    print(f"{cls:22s} {TRUTH[cls]:7s} {nma:14d} {gh:11d} {nfam:13d}")
