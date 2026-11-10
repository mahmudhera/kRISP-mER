#!/usr/bin/env python3
# 02_build_targets_n22.py  CLASS [N]
# ----------------------------------------------------------------------------
# STEP 2 of the clean N22 rebuild — gene lists + centered target windows, all
# in the OsN22 Ensembl namespace (n22_exons.fa from step 1).
#
# For one CLASS:
#   1. align that class's bucket guides (<class>.txt, canonical 23-mers) to
#      n22_exons.fa with bowtie -v 0 (exact)  -> per-guide transcript+position
#   2. write genes_<class>_analysis.txt (unique parent-gene ids)  [record]
#   3. pick N transcripts that carry a class guide; for each, cut a 150bp window
#      centered on that guide -> targets/<class>/target_i.fasta
#
# Because guides, exon FASTA and SAM now share ONE namespace, a class guide is
# GUARANTEED to sit inside its own window. verify_targets.py checks this after.
#
# NOTE on canonical guides: buckets store canonical PAM 23-mers (min of k / rc).
# bowtie aligns the literal sequence, so we feed BOTH the stored strand and its
# reverse-complement, and record whichever maps. The extracted window always
# contains the guide 23-mer on one strand, which is all class_scoring needs
# (it re-canonicalizes every 23-mer in the window).
# ----------------------------------------------------------------------------
import os, sys, subprocess, random
from Bio import SeqIO

CLASS = sys.argv[1]
N     = int(sys.argv[2]) if len(sys.argv) > 2 else 50
WIN, HALF = 150, 75
random.seed(42)                      # reproducible target selection

FRESH   = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
BUCKET  = os.path.expanduser(f"~/krispmer/new_datasets/exps/n22/{CLASS}.txt")
EXONS   = f"{FRESH}/jap_exons.fa"
IDX     = f"{FRESH}/jap_exons_idx"
OUTDIR  = f"{FRESH}/ref_targets/{CLASS}"
SAMOUT  = f"{FRESH}/ref_{CLASS}_exons_max_0.sam"
GENES   = f"{FRESH}/ref_genes_{CLASS}_analysis.txt"
os.makedirs(OUTDIR, exist_ok=True)

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

# ---- 1. build a bowtie query FASTA of the class guides (both strands) -------
qry = f"{FRESH}/ref_{CLASS}_guides.fa"
ng = 0
with open(BUCKET) as fh, open(qry,"w") as out:
    for ln in fh:
        g = ln.split()[0].strip().upper()
        if len(g) != 23: continue
        out.write(f">{ng}|f\n{g}\n>{ng}|r\n{rc(g)}\n")
        ng += 1
print(f"[{CLASS}] bucket guides: {ng}")

# ---- 2. bowtie exact alignment to the transcript source --------------------
# -v 0 exact, -a all hits, -f FASTA query, --sam ; suppress if >... keep simple
print(f"[{CLASS}] bowtie -v 0 vs jap_exons ...")
with open(SAMOUT,"w") as so:
    subprocess.run(["bowtie","-v","0","-a","--sam","-f","-p",str(os.cpu_count() or 4),
                    IDX, qry], stdout=so, stderr=subprocess.DEVNULL, check=True)

# ---- 3. parse SAM: transcript -> first guide position; collect gene ids -----
# SAM cols: qname flag rname pos ... (0-based flag&16 = reverse)
pos_by_tx = {}                 # tx -> (pos0based, guide_seq_on_plus)
genes = set()
with open(SAMOUT) as fh:
    for ln in fh:
        if ln.startswith("@"): continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 10: continue
        flag = int(c[1])
        if flag & 4: continue                      # unmapped
        tx, pos, seq = c[2], int(c[3]), c[9]       # pos is 1-based leftmost
        if tx not in pos_by_tx:
            pos_by_tx[tx] = (pos-1, seq)           # store 0-based start
        # japonica RefSeq tx id = XM_015766610.3 -> gene = strip the .version
        # (there is no OsN22-style _gene suffix here; use the accession w/o version)
        gene = tx.rsplit(".",1)[0]
        genes.add(gene)

with open(GENES,"w") as g:
    for x in sorted(genes): g.write(x+"\n")
print(f"[{CLASS}] transcripts hit: {len(pos_by_tx)}  parent genes: {len(genes)}")

# ---- 4. load exon seqs for the hit transcripts, center 150bp windows -------
want = set(pos_by_tx)
seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(EXONS,"fasta") if r.id in want}

txs = [t for t in pos_by_tx if t in seqs]
random.shuffle(txs)
done = 0
for tx in txs:
    if done >= N: break
    start0, gseq = pos_by_tx[tx]
    s = seqs[tx]
    center = start0 + 11                      # middle of the 23-mer
    a = max(0, center - HALF)
    b = min(len(s), a + WIN)
    a = max(0, b - WIN)
    window = s[a:b]
    if len(window) < 40: continue             # too-short transcript, skip
    done += 1
    with open(f"{OUTDIR}/target_{done}.fasta","w") as o:
        o.write(f">{tx}|guidepos{start0-a}\n{window}\n")

print(f"[{CLASS}] wrote {done} target windows -> {OUTDIR}")
