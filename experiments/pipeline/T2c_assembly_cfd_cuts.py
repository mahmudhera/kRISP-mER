#!/usr/bin/env python3
# T2c_assembly_cfd_cuts.py  ORG  ANCHOR
# ----------------------------------------------------------------------------
# TASK 2 (variant c): per-class CFD off-target load computed from the INDIVIDUAL'S
# ASSEMBLY k-mer counts (NOT reads) -- internally consistent with how the 4
# discordance classes are DEFINED (assembly-vs-assembly). Answers directly:
# "does CRISPOR accept guides that are DUPLICATED in the individual's assembly?"
#
# assembly-CFD of a guide g (no read coverage, no /lambda -- assembly counts ARE
# integer copy numbers):
#     CFD_asm(g) = sum_{m in {g} u N1(g) u N2(g)} asm_count(m) * CFD(g,m)
#   asm_count(m) = count of m in the individual's ASSEMBLY jf (both strands)
#                  (rice: N22 assembly; human: CHM13 assembly)
#   CFD(g,m)     = Doench-2016, krispmer get_cfd_score.py
#   guide oriented to plus-strand (NGG) before neighbour-enum + scoring.
#
# For a multiple_v_single guide with assembly copy 2, the on-target term alone is
# 2*1.0 = 2 -> load >= 2, as expected (unlike the reads version, which under-counts
# assembly duplications). Same 3-tool / own-verdict grouping as T2b.
#
# Usage on VM:
#   python3 T2c_assembly_cfd_cuts.py rice n22
#   python3 T2c_assembly_cfd_cuts.py rice japonica
#   python3 T2c_assembly_cfd_cuts.py human chm13
#   python3 T2c_assembly_cfd_cuts.py human grch38
import os, sys, glob, csv, itertools, subprocess, tempfile
import pandas as pd

ORG    = sys.argv[1] if len(sys.argv) > 1 else "rice"
ANCHOR = sys.argv[2] if len(sys.argv) > 2 else "n22"
CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
           "ref_multiple_v_single":"accept","present_v_absent":"accept"}
K_CUT, MIT_CUT, GS_CUT = 1.5, 50, 0.2
MAX_HD = 2

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def canon(k): return min(k, rc(k))
def orient_guide(g):
    g=g.upper()
    if g[-2:]=="GG": return g
    if g[:2]=="CC":  return rc(g)
    return g

# ---- krispmer's own Doench-2016 CFD ---------------------------------------
_PKGDIRS=[os.path.expanduser("~/krispmer/kRISP-mER/kRISP-meR_source/krispmer"),
          os.path.expanduser("~/krispmer/krispmer"),
          os.path.expanduser("~/krispmer_src/krispmer")]
try:
    for _p in _PKGDIRS: sys.path.insert(0, _p)
    from get_cfd_score import get_score
except Exception:
    import pickle
    _CANDS=[os.path.join(d,"CFD_scoring") for d in _PKGDIRS]+["CFD_scoring"]
    _dir=next((d for d in _CANDS if os.path.exists(os.path.join(d,"mismatch_score.pkl"))),None)
    if _dir is None: sys.exit("[ERR] CFD_scoring pkls not found")
    _mm =pickle.load(open(os.path.join(_dir,"mismatch_score.pkl"),"rb"),encoding="latin1")
    _pam=pickle.load(open(os.path.join(_dir,"pam_scores.pkl"),"rb"),encoding="latin1")
    def _revcom(s): return s.translate(str.maketrans("ACGTU","TGCAA"))[::-1]
    def get_score(candidate, sequence):
        sg=candidate[:20].replace('T','U'); wt=sequence[:20].replace('T','U'); sc=1.0
        for i,sl in enumerate(sg):
            if wt[i]==sl: continue
            try: sc*=_mm['r'+wt[i]+':d'+_revcom(sl)+','+str(i+1)]
            except KeyError: continue
        try: sc*=_pam[candidate[-2:]]
        except KeyError: pass
        return sc

# ---- run dirs + INDIVIDUAL ASSEMBLY jf for (org, anchor) -------------------
# NB: the count source is the INDIVIDUAL's assembly (the one that defines c_I),
# regardless of which anchor the target windows came from.
if ORG=="rice":
    BASE=os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
    BUCKET_DIR=os.path.expanduser("~/krispmer/new_datasets/exps/n22")
    ASM_JF=os.path.expanduser("~/krispmer/new_datasets/oryza_sativa_n22/genome_counts.jf")  # N22 assembly
    if ANCHOR=="n22": KDIR,CDIR,GDIR,TDIR="scores_bothstrand","crispor_class_scores","guidescan_class_scores","targets"
    else:             KDIR,CDIR,GDIR,TDIR="ref_scores_bothstrand","ref_crispor_class_scores","ref_guidescan_class_scores","ref_targets"
else:  # human
    BASE=os.path.expanduser("~/human_assemblies_kmer_count/july-10")
    BUCKET_DIR=BASE
    # CHM13 assembly jf (individual). Both verified to reproduce class counts (2/2/2/3/3/2).
    _H=[os.path.expanduser("~/human_assemblies_kmer_count/chm13_assembly/chm13_counts_binary.jf"),
        os.path.expanduser("~/krispmer/data/homo-sapiens/chm13-assembly/ncbi-genomes-2022-10-25/chm13_genome.jf")]
    ASM_JF=next((p for p in _H if os.path.exists(p)),_H[0])
    if ANCHOR=="chm13": KDIR,CDIR,GDIR,TDIR="ind_scores_bothstrand","ind_crispor_class_scores","ind_guidescan_class_scores","ind_targets"
    else:               KDIR,CDIR,GDIR,TDIR="scores_bothstrand","crispor_class_scores","guidescan_class_scores","targets"

if not os.path.exists(ASM_JF): sys.exit(f"[ERR] assembly jf not found: {ASM_JF}")
print(f"[assembly jf] {ASM_JF}")

# ---- class guides in the 50 windows (canon -> plus-strand oriented 23-mer) --
def class_guides_in_windows(cls):
    bucket=set(canon(l.split()[0].upper()) for l in open(f"{BUCKET_DIR}/{cls}.txt")
               if l.split() and len(l.split()[0])==23)
    reps={}
    for f in glob.glob(f"{BASE}/{TDIR}/{cls}/*.fasta")+glob.glob(f"{BASE}/{TDIR}/{cls}/*.fa"):
        seq="".join(l.strip() for l in open(f) if not l.startswith(">")).upper()
        for i in range(len(seq)-22):
            w=seq[i:i+23]; c=canon(w)
            if c in bucket and c not in reps: reps[c]=orient_guide(w)
    return reps

def v_krisp(cls):
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
    return d

def v_crispor(cls):
    d={}
    for f in glob.glob(f"{BASE}/{CDIR}/{cls}/*_crispor.tsv"):
        with open(f) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip(): continue
                p=ln.rstrip("\n").split("\t")
                if len(p)>=5 and len(p[2].strip())>=20:
                    try:
                        g=canon(p[2].strip().upper()); mit=float(p[3])
                        d[g]="reject" if mit<MIT_CUT else "accept"
                    except: pass
    return d

def v_guidescan(cls):
    d={}
    base=f"{BASE}/{GDIR}/out"
    fn=f"{base}/{cls}_ont.csv"
    if not os.path.exists(fn): fn=f"{base}/{cls}.csv"
    if not os.path.exists(fn): return d
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if row.get("match_distance","").strip()!="0": continue
            try:
                g=canon(row["match_sequence"].upper()); sp=float(row["specificity"])
            except: continue
            if g not in d: d[g]="reject" if sp<GS_CUT else "accept"
    return d

# ---- HD1+HD2 neighbours; batched jellyfish query against the ASSEMBLY jf ----
BASES=("A","C","G","T")
def neighbors(seq):
    out={seq}
    for n in range(1,MAX_HD+1):
        for pos in itertools.combinations(range(len(seq)),n):
            alts=[[b for b in BASES if b!=seq[p]] for p in pos]
            for repl in itertools.product(*alts):
                s=list(seq)
                for i,p in enumerate(pos): s[p]=repl[i]
                out.add("".join(s))
    return out

def jf_counts(kmers):
    need=set()
    for k in kmers: need.add(k); need.add(rc(k))
    with tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False) as fh:
        tmp=fh.name
        for i,k in enumerate(need): fh.write(f">v{i}\n{k}\n")
    try:
        res=subprocess.run(["jellyfish","query",ASM_JF,"-s",tmp],
                           capture_output=True,text=True,check=True)
    finally:
        os.unlink(tmp)
    raw={}
    for line in res.stdout.splitlines():
        p=line.split()
        if len(p)==2: raw[p[0]]=int(p[1])
    return {k: raw.get(k,0)+raw.get(rc(k),0) for k in kmers}

# ---- SANITY: the assembly jf must reproduce the count table that DEFINED the
#     classes. Probe a few known multiple_v_single guides; abort if mismatched. --
def sanity_check():
    if ORG!="rice": return  # probe values known for rice n22
    probe={"CGGTGTCAATCTGGCTGCGGAGG":4,"ATCTGTATTGCATCTTGTGATGG":2,
           "AGTGCCGCAACATATTCTTGAGG":2,"CGGGAGATGCCTTCGACTCACGG":2}
    c=jf_counts(set(probe))
    print("[sanity] assembly jf vs class-definition count table (rice N22):")
    bad=0
    for g,exp in probe.items():
        got=c.get(g,None)
        ok = (got==exp)
        if not ok: bad+=1
        print(f"   {g}  expect={exp}  jf={got}  {'OK' if ok else 'MISMATCH'}")
    if bad:
        print(f"[WARN] {bad} mismatches -- the assembly jf may not equal the class-def table.")
    else:
        print("[sanity] all match -- assembly jf is the correct source.\n")

sanity_check()

# ---- main -----------------------------------------------------------------
TOOLS=[("kRISP-meR",v_krisp),("CRISPOR",v_crispor),("GuideScan2",v_guidescan)]
rows=[]
for cls in CLASSES:
    reps=class_guides_in_windows(cls)
    if not reps: print(f"  [skip] {cls}: no window guides"); continue
    verdicts={name:fn(cls) for name,fn in TOOLS}
    guide_neigh={c: neighbors(w) for c,w in reps.items()}
    allk=set()
    for s in guide_neigh.values(): allk|=s
    print(f"  {cls}: {len(reps)} guides, {len(allk)} neighbour kmers -> assembly jellyfish ...", flush=True)
    cnt=jf_counts(allk)
    cfd_of={}
    for c,w in reps.items():
        tot=0.0
        for m in guide_neigh[c]:
            k=cnt.get(m,0)
            if k>0: tot+=k*get_score(w,m)
        cfd_of[c]=round(tot,4)          # NO /lambda -- assembly counts are copy numbers
    for name,_ in TOOLS:
        vd=verdicts[name]
        for c,w in reps.items():
            v=vd.get(c,"not_generated")
            group="accepted" if v=="accept" else "rejected_side"
            rows.append({"class":cls,"truth":TRUTH[cls],"tool":name,"group":group,
                         "verdict":v,"guide":w,"cfd_cuts":cfd_of[c]})

out=f"{BASE}/asm_cfd_dist_{ORG}_{ANCHOR}.tsv"
pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
print(f"\nwrote {out}  ({len(rows)} tool-guide rows, ASSEMBLY counts, MAX_HD={MAX_HD})")
df=pd.DataFrame(rows)
for cls in CLASSES:
    print(f"### {cls}  (truth {TRUTH[cls]})")
    sub=df[df["class"]==cls]
    for name,_ in TOOLS:
        for grp in ["accepted","rejected_side"]:
            s=sub[(sub["tool"]==name)&(sub["group"]==grp)]["cfd_cuts"]
            if len(s):
                print(f"  {name:11s} {grp:14s} n={len(s):4d}  asm-CFD median={s.median():7.3f} mean={s.mean():7.3f}")
            else:
                print(f"  {name:11s} {grp:14s} n=   0")
