#!/usr/bin/env python3
# T1_famous_genes_rice.py
# ----------------------------------------------------------------------------
# TASK 1 (rice), REVERSE route -- no OsN22<->japonica ID map / no BLAST needed.
# Rice class genes are OsN22 (aus cultivar); trait-gene DBs (funRiceGenes) use
# japonica MSU (LOC_Os) IDs -- there is NO crosswalk. So instead of
# guide->OsN22->trait, we go trait-gene -> its japonica CDS -> guides in it ->
# do any match a MIS-ACCEPTED dangerous guide?  The guide SEQUENCE is the bridge
# (assembly-independent), so no ID mapping is required.
#
# Pipeline:
#   1. funRiceGenes trait genes -> MSU locus (funrice_msu2sym.tsv)
#   2. that locus' CDS sequence  (Osativa_323_v7.0.cds_primaryTranscriptOnly.fa,
#      Phytozome MSU v7.0; header ">LOC_Os..g.. ... locus=LOC_Os..g..")
#   3. extract every NGG 23-mer guide from the CDS (both strands)
#   4. intersect (canonical seq) with MIS-ACCEPTED dangerous guides from the rice
#      N22 assembly-CFD run: asm-CFD>=2 (duplicated in individual) AND a reference
#      tool (CRISPOR/GuideScan2) ACCEPTS it AND kRISP-meR does NOT accept.
#   5. for hits: N22 (individual) vs japonica (reference) 23-mer copy count.
#
# Env (run ON VM):
#   DATA   = dir with funrice_msu2sym.tsv (scp'd)              [~/task1_rice]
#   CDS    = Osativa MSU CDS fasta
#   ASM    = asm_cfd_dist_rice_n22.tsv
#   N22_JF = N22 assembly jf (individual)      GRCH-analogue = japonica
#   JAP_JF = japonica assembly jf (reference)
import os, sys, csv, re, subprocess, tempfile
from collections import defaultdict

DATA   = os.path.expanduser(os.environ.get("DATA", "/mnt/windows-data/Thesis/data"))
CDS    = os.path.expanduser(os.environ.get("CDS",
          "~/krispmer/new_datasets/assemblies/reference/annotation/Osativa_323_v7.0.cds_primaryTranscriptOnly.fa"))
ASM    = os.path.expanduser(os.environ.get("ASM",
          "~/krispmer/new_datasets/exps/n22/fresh/asm_cfd_dist_rice_n22.tsv"))
N22_JF = os.path.expanduser(os.environ.get("N22_JF",
          "~/krispmer/new_datasets/oryza_sativa_n22/genome_counts.jf"))
JAP_JF = os.path.expanduser(os.environ.get("JAP_JF",
          "~/krispmer/new_datasets/assemblies/reference/mer_counts.jf"))
DUP_CUT = 2.0
REF_TOOLS = ["CRISPOR","GuideScan2"]

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))
def _find(name):
    for p in (f"{DATA}/{name}", f"{DATA}/famous_genes/{name}"):
        if os.path.exists(p): return p
    sys.exit(f"[ERR] cannot find {name} under {DATA}")

# ---- 1. funRiceGenes MSU-locus -> trait symbol ----------------------------
msu2sym={}
for ln in open(_find("funrice_msu2sym.tsv")):
    p=ln.rstrip("\n").split("\t")
    if len(p)==2 and p[0]!="msu_locus": msu2sym[p[0]]=p[1]
trait_loci=set(msu2sym)
print(f"funRiceGenes trait loci: {len(trait_loci)}")

# ---- 4a. mis-accepted dangerous guides (rice N22 asm-CFD) ------------------
def mis_accepted():
    G=defaultdict(lambda:{"cfd":None,"verd":{},"seq":None})
    for r in csv.DictReader(open(ASM),delimiter="\t"):
        g=canon(r["guide"].upper())
        G[g]["cfd"]=float(r["cfd_cuts"]); G[g]["seq"]=r["guide"].upper()
        G[g]["verd"][r["tool"]]=r["verdict"]
    out={}
    for g,d in G.items():
        if d["cfd"] is None or d["cfd"]<DUP_CUT: continue
        if d["verd"].get("kRISP-meR")=="accept": continue
        tools=[t for t in REF_TOOLS if d["verd"].get(t)=="accept"]
        if not tools: continue
        out[g]={"seq":d["seq"],"cfd":d["cfd"],"tools":tools}
    return out
MA=mis_accepted()
mis_canon=set(MA)
print(f"mis-accepted dangerous guides (rice N22, all classes pooled): {len(MA)}")

# ---- 2+3. stream trait-gene CDS -> NGG guides -> match mis-accepted --------
LOC_RE=re.compile(r"(LOC_Os\d+g\d+)")
def ngg_guides(seq):
    """all canonical 23-mers ending in .GG (plus strand) or starting CC. (minus)"""
    seq=seq.upper(); out=set()
    n=len(seq)
    for i in range(n-22):
        w=seq[i:i+23]
        if 'N' in w: continue
        if w[-2:]=="GG": out.add(canon(w))            # plus-strand NGG
    for i in range(n-22):                              # minus strand: CCN..
        w=seq[i:i+23]
        if 'N' in w: continue
        if w[:2]=="CC": out.add(canon(rc(w)))
    return out

# read the CDS fasta, keep only trait-gene loci, collect guide->loci
guide2loci=defaultdict(set)
cur_loc=None; cur=[]
def flush(loc,chunks):
    if loc is None or loc not in trait_loci: return
    seq="".join(chunks)
    for cg in ngg_guides(seq):
        if cg in mis_canon: guide2loci[cg].add(loc)
with open(CDS) as fh:
    for ln in fh:
        if ln.startswith(">"):
            flush(cur_loc,cur)
            m=LOC_RE.search(ln); cur_loc=m.group(1) if m else None; cur=[]
        else:
            cur.append(ln.strip())
    flush(cur_loc,cur)
print(f"mis-accepted guides that fall in a trait-gene CDS: {len(guide2loci)}")

# ---- 5. N22 vs japonica copy counts for the hit guides --------------------
hits=sorted(guide2loci)
def jfq(jf,seqs):
    need=set()
    for cg in seqs:
        w=MA[cg]["seq"]; need.add(w); need.add(rc(w))
    with tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False) as fhh:
        tmp=fhh.name
        for i,k in enumerate(need): fhh.write(f">v{i}\n{k}\n")
    r=subprocess.run(["jellyfish","query",jf,"-s",tmp],capture_output=True,text=True); os.unlink(tmp)
    raw={}
    for line in r.stdout.splitlines():
        x=line.split()
        if len(x)==2: raw[x[0]]=int(x[1])
    return {cg: raw.get(MA[cg]["seq"],0)+raw.get(rc(MA[cg]["seq"]),0) for cg in seqs}
n22=jfq(N22_JF,hits) if hits else {}
jap=jfq(JAP_JF,hits) if hits else {}

rows=[]
for cg in hits:
    loci=sorted(guide2loci[cg])
    syms=sorted({msu2sym[l] for l in loci})
    for l in loci:
        rows.append({"trait_gene":msu2sym[l],"msu_locus":l,"guide":MA[cg]["seq"],
                     "N22_copies":n22.get(cg,0),"japonica_copies":jap.get(cg,0),
                     "asm_cfd":round(MA[cg]["cfd"],3),
                     "recommended_by":";".join(MA[cg]["tools"])})
import pandas as pd
out=f"{DATA}/task1_rice_n22.tsv"
pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
print(f"\nwrote {out}  ({len(rows)} guide-in-trait-gene rows)")
print(f'\n{"trait":14s} {"MSU":16s} {"N22":>3s} {"jap":>3s} {"asmCFD":>6s} {"by":18s} guide')
for r in sorted(rows,key=lambda r:-r["N22_copies"]):
    print(f'{r["trait_gene"][:14]:14s} {r["msu_locus"]:16s} {r["N22_copies"]:3d} {r["japonica_copies"]:3d} '
          f'{r["asm_cfd"]:6.2f} {r["recommended_by"]:18s} {r["guide"]}')
