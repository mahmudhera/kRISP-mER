#!/usr/bin/env python3
# T3_handpicked_ranks.py  ORG   (human | rice)
# -----------------------------------------------------------------------------
# For each HANDPICKED guide (fixed/task1_handpicked_examples_*.txt), report where
# each tool ranked that guide WITHIN ITS OWN TARGET REGION -- i.e. "oder respective
# target region e koy number e chilo".
#
# Per-region rank, each tool by its own score:
#   kRISP-meR  : Expected_cuts_in_genome  ascending  (lower = safer = rank 1)
#   CRISPOR    : MIT specificity (col 4)  descending (higher = better = rank 1)
#   GuideScan2 : specificity              descending (higher = better = rank 1)
#
# A "target region" = one 150bp window fasta in targets/<class>/<region>.fasta.
# A region's guide set = all NGG 23-mers found in that window. A tool's rank for a
# guide is its position among the region's guides that the tool actually scored.
#
# Run this ON THE VM (raw per-region tool outputs live there).
# -----------------------------------------------------------------------------
import glob, csv, os, sys
import pandas as pd

ORG = sys.argv[1] if len(sys.argv) > 1 else "human"

# ---- per-organism VM paths (mirror 08b_*_bothstrand.py) ----
if ORG == "human":
    ROOT   = os.path.expanduser("~/human_assemblies_kmer_count/july-10")
    HANDPK = os.path.expanduser("~/task1_handpicked_examples_human.txt")   # or wherever you keep it
elif ORG == "rice":
    ROOT   = os.path.expanduser("~/krispmer/new_datasets/exps/n22/fresh")
    HANDPK = os.path.expanduser("~/task1_handpicked_examples_rice.txt")
else:
    sys.exit("ORG must be human|rice")

CLASS  = "multiple_v_single"    # every handpicked example is this class
K_CUT, MIT_CUT, GS_CUT = 1.5, 50, 0.2

def canon(k):
    k = k.upper()
    return min(k, k.translate(str.maketrans("ACGT", "TGCA"))[::-1])

# ---- handpicked guides: list of (gene, orig_guide, canon_guide) ----
# A handpicked row is a line that contains BOTH a leading gene symbol and a 23-mer
# guide. We keep the original (as-written) guide for display and canon() for lookup.
def load_handpicked(path):
    out = []
    for ln in open(path):
        toks = ln.replace(",", " ").split()
        if not toks: continue
        # the guide: the single 23-mer ACGT token on the line
        guides = [t.upper() for t in toks if len(t) == 23 and set(t.upper()) <= set("ACGT")]
        if not guides: continue
        # the gene: first token, only if it's a plausible symbol (not a comment/rule)
        gene = toks[0]
        if gene.startswith(("#", "-", "gene", "trait", "asm", "Copy", "(")): continue
        for g in guides:
            out.append((gene, g, canon(g)))
    return out

def revcom(s): return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# ---- a region's guides, BOTH strands, as ORIENTED NGG 23-mers ----
# Returns {canon_key: oriented_ngg_23mer}. From the oriented 23-mer we can derive
# the GuideScan2 protospacer (first 20 nt) for that guide.
def region_guides_oriented(fasta):
    seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()
    d = {}
    for i in range(len(seq) - 22):
        sub = seq[i:i+23]
        if set(sub) > set("ACGT"): continue
        if sub[-2:] == "GG":                 # plus-strand NGG guide
            d[canon(sub)] = sub
        if sub[:2] == "CC":                  # minus-strand NGG guide
            d[canon(sub)] = revcom(sub)
    return d

# ---- a region's guide set: canon 23-mers in that window, BOTH strands ----
def region_guides(fasta):
    seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()
    gs = set()
    for i in range(len(seq) - 22):
        sub = seq[i:i+23]
        if set(sub) > set("ACGT"): continue
        # plus-strand NGG guide, OR minus-strand NGG guide (CCN on the plus strand)
        if sub[-2:] == "GG" or sub[:2] == "CC":
            gs.add(canon(sub))
    return gs

# ---- per-region tool score maps: {canon_guide: score} ----
def krisp_region(region_id):
    # exact region file: target_N.csv  (no prefix-collision with target_N0..N9)
    d = {}
    f = f"{ROOT}/scores_bothstrand/{CLASS}/{region_id}.csv"
    if os.path.exists(f):
        try: df = pd.read_csv(f, comment="#", skip_blank_lines=True)
        except: return d
        if "tgt_in_plus" not in df.columns: return d
        for _, r in df.iterrows():
            try: d[canon(str(r["tgt_in_plus"]))] = float(r["Expected_cuts_in_genome"])
            except: pass
    return d   # lower is better

def crispor_region(region_id):
    # exact region file: target_N_crispor.tsv  (NOT target_N_off.tsv)
    d = {}
    f = f"{ROOT}/crispor_class_scores/{CLASS}/{region_id}_crispor.tsv"
    if os.path.exists(f):
        for ln in open(f):
            if ln.startswith("#") or not ln.strip(): continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 5 and len(p[2].strip()) >= 20:
                try: d[canon(p[2].strip())] = float(p[3])
                except: pass
    return d   # higher is better

# GuideScan2 is one pooled file per class. A guide's identity is its QUERY sequence
# (col `sequence`, which carries an 'N' in the PAM); each guide appears in many rows
# (one per genomic match). specificity is a per-guide property, so we key on the
# 20-nt protospacer (sequence[:20], PAM-independent) and take the guide's own value.
_GS_CACHE = {}
def guidescan_all():
    if CLASS in _GS_CACHE: return _GS_CACHE[CLASS]
    d = {}   # {protospacer20 : specificity}
    fn = f"{ROOT}/guidescan_class_scores/out/{CLASS}_ont.csv"
    if os.path.exists(fn):
        for row in csv.DictReader(open(fn)):
            proto = row["sequence"].upper()[:20]      # query protospacer (before the N-PAM)
            if len(proto) < 20 or set(proto) > set("ACGT"): continue
            d[proto] = float(row["specificity"])       # same value across the guide's rows
    _GS_CACHE[CLASS] = d
    return d   # higher is better; keyed by protospacer

def rank_in(scoremap, region_set, guide, higher_better, accepts):
    """Ranks of `guide` among this tool's region guides, two denominators:
       (a) among ALL scored guides in the region,
       (b) among only the guides the tool ACCEPTS (`accepts` = fn score->bool).
    Returns rank_all, n_scored, score, rank_acc, n_accept.
    rank_acc is None when the guide itself is rejected (not in the accepted set)."""
    scored = {g: s for g, s in scoremap.items() if g in region_set}
    accepted = {g: s for g, s in scored.items() if accepts(s)}
    n_scored, n_accept = len(scored), len(accepted)
    if guide not in scored:
        return None, n_scored, None, None, n_accept
    score = scored[guide]
    order = sorted(scored.items(), key=lambda kv: kv[1], reverse=higher_better)
    rank_all = next(i for i, (g, _) in enumerate(order, 1) if g == guide)
    rank_acc = None
    if guide in accepted:
        order_a = sorted(accepted.items(), key=lambda kv: kv[1], reverse=higher_better)
        rank_acc = next(i for i, (g, _) in enumerate(order_a, 1) if g == guide)
    return rank_all, n_scored, score, rank_acc, n_accept

# -----------------------------------------------------------------------------
handpicked = load_handpicked(HANDPK)
gs_all = guidescan_all()                       # keyed by 20-nt protospacer
region_fastas = glob.glob(f"{ROOT}/targets/{CLASS}/*.fasta") + \
                glob.glob(f"{ROOT}/targets/{CLASS}/*.fa")
# per region: canon-key set (kRISP/CRISPOR) and canon->oriented map (for GS2 protospacers)
REGORI = {os.path.splitext(os.path.basename(fa))[0]: region_guides_oriented(fa) for fa in region_fastas}
REGSET = {rid: set(m.keys()) for rid, m in REGORI.items()}

rows = []
for gene, orig, cguide in handpicked:
    hit_regions = [rid for rid, rset in REGSET.items() if cguide in rset]
    if not hit_regions:
        rows.append({"gene": gene, "guide": orig, "region": "NOT_FOUND"})
        continue
    for rid in hit_regions:
        rset = REGSET[rid]
        # GuideScan2 works in protospacer space: region protospacer set + this guide's protospacer
        ori_map = REGORI[rid]
        proto_set = {ori_map[c][:20] for c in rset}
        tgt_proto = ori_map[cguide][:20]
        kr, kn, ks, kar, kna = rank_in(krisp_region(rid),   rset, cguide,
                                       higher_better=False, accepts=lambda s: s <= K_CUT)
        cr, cn, cs, car, cna = rank_in(crispor_region(rid), rset, cguide,
                                       higher_better=True,  accepts=lambda s: s >= MIT_CUT)
        gr, gn, gsc, gar, gna = rank_in(gs_all,             proto_set, tgt_proto,
                                       higher_better=True,  accepts=lambda s: s >= GS_CUT)
        def acc(r, n):  # "rank/accepted" string; '-' means the tool rejected this guide
            return f"{r}/{n}" if r is not None else f"-/{n}"
        rows.append({
            "gene": gene, "guide": orig, "region": rid, "region_guides": len(rset),
            "krispmer_rank": kr, "krispmer_of": kn, "krispmer_acc": acc(kar, kna), "krispmer_Ecuts": ks,
            "crispor_rank": cr,  "crispor_of": cn,  "crispor_acc": acc(car, cna),  "crispor_MITspec": cs,
            "guidescan_rank": gr,"guidescan_of": gn,"guidescan_acc": acc(gar, gna),"guidescan_spec": gsc,
        })

out = pd.DataFrame(rows)
outfn = f"{ROOT}/handpicked_ranks_{ORG}.tsv"
out.to_csv(outfn, sep="\t", index=False)
pd.set_option("display.width", 200, "display.max_columns", 30)
print(out.to_string(index=False))
print(f"\nwrote {outfn}")
