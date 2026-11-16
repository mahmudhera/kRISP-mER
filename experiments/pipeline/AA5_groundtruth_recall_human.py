#!/usr/bin/env python3
# AA5_groundtruth_recall_human.py  --  run ON THE VM.
#
# THE EXPERIMENT (Navid's spec):
#   For each class's 50 targets, build the GROUND-TRUTH set of individual-specific
#   targetable guides = NGG 23-mers that are PRESENT in the CHM13 assembly but
#   ABSENT from GRCh38 (definition A). Then ask, of that ground-truth set, how many
#   did kRISP-meR (reads-based) recover vs AlleleAnalyzer (reads->GRCh38->VCF based).
#   kRISP-meR should recover the individual-specific guides; AA should miss the ones
#   its reference-alignment pipeline never surfaced.
#
# GROUND TRUTH construction (per target, per class):
#   1. target fasta header -> GRCh38 ENST  (>hg38_..._ENST00000396890.6)
#   2. ENST -> CHM13 transcript span via CAT Liftoff GFF (source_transcript=ENST..)
#      col1 is chrN (UCSC); CHM13 fasta uses CP accessions -> map via chr number.
#   3. faidx CHM13 for that span -> CHM13 sequence at the individual's locus
#   4. enumerate every NGG 23-mer (both strands) in that CHM13 sequence
#   5. keep guide G iff  chm13_kmer[G] > 0  AND  grch38_kmer[G] == 0
#      = individual-specific targetable guide (ground truth)
#
# RECALL:
#   kRISP-meR recovered G  <-  G (canonical) in that target's kRISP-meR designed guides
#                             (scores_bothstrand/<cls>/target_N.csv, col tgt_in_plus)
#   AA recovered G          <-  G (canonical) in AA's reconstructed 23-mers for the class
#                             (alleleanalyzer.tsv gRNA_ref/gRNA_alt + genomic PAM)
#
# Output: gt_recall_human.tsv   (per class: |GT|, krisp n/%, aa n/%, both, neither)
import os, csv, re, subprocess, tempfile, glob
from collections import defaultdict

AA   = os.path.expanduser("~/krispmer/experiments/krispmer-expreminets/manuscript-experiments/human/july-13/specialized_targets")
H    = os.path.expanduser("~/human_assemblies_kmer_count")
CHM13_FA = f"{H}/chm13_assembly/cmh13v2.fasta"
GFF      = f"{H}/chm13_assembly/chm13v2.0_CAT_Liftoff.gff3"
CHM13_JF = f"{H}/chm13_assembly/chm13_counts_binary.jf"
GRCH38_JF= f"{H}/grch38.jf"
KRISP_DIR= f"{H}/july-10/scores_bothstrand"
OUTDIR   = os.path.expanduser("~")
CLASSES  = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH    = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
            "ref_multiple_v_single":"accept","present_v_absent":"accept"}
ENST_RE  = re.compile(r"ENST\d+")
_HOMO    = {b*20 for b in "ACGT"}

def rc(s): return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]
def canon(k): return min(k, rc(k))
def proto20(g):
    """reduce any guide to its canonical 20-nt protospacer for cross-tool matching.
    kRISP tgt_in_plus/minus, AA reconstructions, and GT guides are all 23-mers
    (20nt protospacer + NGG PAM), so drop the last 3 (PAM) then canonicalize."""
    g=g.upper()
    if len(g)==23: g=g[:20]
    elif len(g)>20: g=g[:20]
    return canon(g)

# WINDOW-PRESENT-IN-INDIVIDUAL threshold. The window->CHM13 match hamming distances
# are cleanly bimodal: real matches are 0-1 (exact, or the individual's own SNPs);
# non-matches jump to 26-42 (~20-28% divergence over 150bp) with an empty gap between.
# A non-match means the individual does NOT carry a clean copy of this reference locus
# (expected in the ref_* classes: single/absent in individual, multiple/present in ref).
# In that case the CORRECT ground truth is ZERO individual-specific guides -- there is
# nothing to target. So: ALL 50 targets/class stay in the denominator; a non-matching
# target contributes GT=0 (it is NOT removed). This keeps the same 200-target set the
# other experiments used, while refusing to enumerate guides off a wrong/paralog locus.
HAM_GATE = 5

def best_window_in(chm_seq, win):
    """Locate the 150bp GRCh38 window inside the CHM13 transcript sequence.
    Returns (chm_subsequence, best_hamming). Caller applies HAM_GATE. Tries both
    orientations (CAT span may be reverse vs the window)."""
    W=len(win); L=len(chm_seq)
    if L<W: return chm_seq, W        # transcript shorter than window: not a clean match
    best=None; bestd=W+1
    for s in (chm_seq, rc(chm_seq)):
        idx=s.find(win)
        if idx>=0: return s[idx:idx+W], 0
        for i in range(0, L-W+1):
            sub=s[i:i+W]
            d=sum(1 for a,b in zip(sub,win) if a!=b)
            if d<bestd: bestd=d; best=sub
    return (best or ""), bestd

# ---- chrN -> CHM13 CP accession (from fasta headers: ">CP.. chromosome N") -----
def load_chr_map():
    m={}
    with open(CHM13_FA) as fh:
        for ln in fh:
            if not ln.startswith(">"): continue
            acc=ln[1:].split()[0]
            mm=re.search(r"chromosome (\w+)", ln)
            if mm: m["chr"+mm.group(1)]=acc
            if "mitochondrion" in ln: m["chrM"]=acc
        # header scan reads whole file lazily? no -- fasta seq lines are huge; stop early
    return m
def load_chr_map_fast():
    m={}
    r=subprocess.run(["grep",">",CHM13_FA],capture_output=True,text=True)
    for ln in r.stdout.splitlines():
        acc=ln[1:].split()[0]
        mm=re.search(r"chromosome (\w+)", ln)
        if mm: m["chr"+mm.group(1)]=acc
        if "mitochondrion" in ln: m["chrM"]=acc
    return m
CHRMAP=load_chr_map_fast()
print(f"chr map: {len(CHRMAP)} sequences (e.g. chr16 -> {CHRMAP.get('chr16')})")

# ---- ENST -> CHM13 (chrN, start, stop) from CAT GFF (transcript lines) ----------
# an ENST may map to >1 CHM13 locus (paralog copies); keep ALL, we want every place
# the individual carries it.
def load_enst_spans(needed_enst):
    spans=defaultdict(list)   # ENST(no ver) -> list of (chrN,start,stop)
    with open(GFF) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<9 or p[2]!="transcript": continue
            m=re.search(r"source_transcript=(ENST\d+)", p[8])
            if not m: continue
            e=m.group(1)
            if e not in needed_enst: continue
            spans[e].append((p[0],int(p[3]),int(p[4])))
    return spans

# ---- faidx CHM13 region -> sequence --------------------------------------------
def faidx(acc,start,stop):
    r=subprocess.run(["samtools","faidx",CHM13_FA,f"{acc}:{start}-{stop}"],capture_output=True,text=True)
    return "".join(l.strip() for l in r.stdout.splitlines() if not l.startswith(">")).upper()

def ngg_guides(seq):
    seq=seq.upper(); out=set(); n=len(seq)
    for i in range(n-22):
        w=seq[i:i+23]
        if "N" in w: continue
        if w[-2:]=="GG": out.add(w)               # plus strand NGG
        if w[:2]=="CC": out.add(rc(w))            # minus strand CCN -> NGG
    return out

def jf_query(jf, seqs):
    seqs=[s for s in seqs if s]
    if not seqs: return {}
    with tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False) as fh:
        tmp=fh.name
        for i,s in enumerate(seqs): fh.write(f">v{i}\n{s}\n")
    r=subprocess.run(["jellyfish","query",jf,"-s",tmp],capture_output=True,text=True); os.unlink(tmp)
    d={}
    for ln in r.stdout.splitlines():
        x=ln.split()
        if len(x)==2: d[x[0]]=int(x[1])
    return d
def kcount(seq,tbl): return tbl.get(seq,0)+tbl.get(rc(seq),0)

# ---- kRISP-meR designed guides per target (canonical set) ----------------------
def krisp_guides(cls,tgt):
    """kRISP-meR designed guides for a target, as canonical 20-mer protospacers.
    tgt_in_plus AND tgt_in_minus are both 23-mers (20+PAM)."""
    f=f"{KRISP_DIR}/{cls}/{tgt}.csv"
    S=set()
    if not os.path.exists(f): return S
    try:
        import pandas as pd
        df=pd.read_csv(f,comment="#",skip_blank_lines=True)
        for col in ("tgt_in_plus","tgt_in_minus"):
            if col in df.columns:
                for v in df[col].dropna():
                    s=str(v).upper()
                    if len(s)>=20: S.add(proto20(s))
    except Exception: pass
    return S

# ---- AA reconstructed 23-mers per class (canonical set) -------------------------
HG38="~/krispmer/new_datasets/feb-23-rice-n22/crispor/crisporWebsite/genomes/hg38/hg38.fa"
HG38=os.path.expanduser(HG38)
def hg38_faidx(region):
    r=subprocess.run(["samtools","faidx",HG38,region],capture_output=True,text=True)
    return "".join(l.strip() for l in r.stdout.splitlines() if not l.startswith(">")).upper()
def aa_guides(cls):
    f=f"{AA}/{cls}/alleleanalyzer.tsv"; S=set()
    if not os.path.exists(f): return S
    for r in csv.DictReader(open(f),delimiter="\t"):
        for proto in (r["gRNA_ref"].upper(), r["gRNA_alt"].upper()):
            if proto in _HOMO: continue
            try:
                st,sp,strand=int(r["start"]),int(r["stop"]),r["strand"]
                if strand=="positive": g=proto+hg38_faidx(f'{r["chrom"]}:{sp+1}-{sp+3}')
                else: g=proto+rc(hg38_faidx(f'{r["chrom"]}:{st-2}-{st}'))
                if len(g)==23: S.add(proto20(g))
            except Exception: pass
    return S

# ================================ main ==========================================
rows=[]
print(f"\n{'class':22s} {'truth':6s} {'|GT|':>6s} {'kRISP':>12s} {'AA':>12s} {'both':>6s} {'neither':>7s}")
for cls in CLASSES:
    # target -> (ENST, 150bp GRCh38 window sequence)
    tgt_enst={}; tgt_win={}
    for f in glob.glob(f"{AA}/{cls}/target_*.fasta"):
        tgt=os.path.basename(f)[:-6]  # target_N
        lines=open(f).read().splitlines()
        hdr=lines[0]; win="".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()
        m=ENST_RE.search(hdr)
        if m: tgt_enst[tgt]=m.group(0); tgt_win[tgt]=win
    needed=set(tgt_enst.values())
    spans=load_enst_spans(needed)

    # build ground-truth guides per target:
    #   the 150bp GRCh38 window, located inside the CHM13 transcript span (aligned),
    #   -> enumerate NGG guides from the CHM13 window -> keep those present in CHM13
    #      k-mers but ABSENT from GRCh38 k-mers. Stored as canonical 20-mer protospacers.
    gt_by_tgt={}; all_gt=set(); n_gated=0; n_nospan=0
    for tgt,enst in tgt_enst.items():
        win=tgt_win[tgt]
        loci=spans.get(enst,[])
        if not loci: n_nospan+=1; gt_by_tgt[tgt]=set(); continue
        # among the ENST's CAT loci, take the CHM13 window with the BEST hamming match;
        # apply the paralog/wrong-locus gate before enumerating guides.
        best_win=""; best_d=len(win)+1
        for (chrN,st,sp) in loci:
            acc=CHRMAP.get(chrN)
            if not acc: continue
            w_chm,d=best_window_in(faidx(acc,st,sp), win)
            if d<best_d: best_d=d; best_win=w_chm
        if best_d>HAM_GATE or not best_win:   # individual lacks a clean copy -> GT=0 (target kept)
            n_gated+=1; gt_by_tgt[tgt]=set(); continue
        cand=ngg_guides(best_win)
        if not cand: gt_by_tgt[tgt]=set(); continue
        chm=jf_query(CHM13_JF,cand); grc=jf_query(GRCH38_JF,cand)
        gt={proto20(g) for g in cand if kcount(g,chm)>0 and kcount(g,grc)==0}
        gt_by_tgt[tgt]=gt; all_gt|=gt

    # recall
    aa_set=aa_guides(cls)     # class-level (AA output is pooled per class)
    krisp_hit=aa_hit=both=neither=0
    for tgt,gt in gt_by_tgt.items():
        kg=krisp_guides(cls,tgt)
        for g in gt:
            k = g in kg
            a = g in aa_set
            krisp_hit+=k; aa_hit+=a
            if k and a: both+=1
            if not k and not a: neither+=1
    N=len(all_gt)
    n_used=sum(1 for t in gt_by_tgt if gt_by_tgt[t])
    def pc(x): return f"{x} ({100*x//N if N else 0}%)"
    print(f"{cls:22s} {TRUTH[cls]:6s} {N:6d} {pc(krisp_hit):>12s} {pc(aa_hit):>12s} {both:6d} {neither:7d}"
          f"   [all 50 kept | individual-lacks-locus:{n_gated} no-CAT-span:{n_nospan} targets-with-GT:{n_used}]")
    rows.append({"class":cls,"truth":TRUTH[cls],"n_targets":len(tgt_enst),"ground_truth_guides":N,
                 "krispmer_recovered":krisp_hit,"krispmer_pct":round(100*krisp_hit/N,1) if N else 0,
                 "aa_recovered":aa_hit,"aa_pct":round(100*aa_hit/N,1) if N else 0,
                 "both":both,"neither":neither,
                 "targets_individual_lacks_locus":n_gated,"targets_no_cat_span":n_nospan,"targets_with_gt":n_used})

out=f"{OUTDIR}/gt_recall_human.tsv"
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"\nwrote {out}")
