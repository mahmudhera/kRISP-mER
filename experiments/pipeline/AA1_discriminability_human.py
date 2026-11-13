#!/usr/bin/env python3
# AA1_discriminability_human.py  --  run ON THE VM.
# AlleleAnalyzer head-to-head metric = per-class LOCUS DISCRIMINABILITY.
#
# AlleleAnalyzer does NOT emit a specificity score on our class guides. It emits,
# per heterozygous variant, an allele-discriminating guide PAIR (gRNA_ref / gRNA_alt)
# with columns: chrom start stop ref alt variant_position_in_guide gRNA_ref gRNA_alt
#               variant_position strand cas_type locus id guide_id
# The natural head-to-head axis is therefore COVERAGE, not accept/reject:
#   of the ENST loci spanned by a class's 50 target windows, for how many can
#   AlleleAnalyzer build at least one allele-specific guide?
#
# Bridge = ENST id. Target fasta headers look like:
#   >hg38_wgEncodeGencodeBasicV38_ENST00000396890.6
# AA 'locus' column = ENST00000272091.8 (same ENST, version may differ) -> match on
# the version-STRIPPED accession (ENST\d+), which is how the rest of the pipeline
# bridges guides->ENST already.
#
# Output: alleleanalyzer_discriminability_human.tsv  (one row per class)
#   class  truth  n_window_loci  n_aa_loci  n_covered  pct_covered  n_guide_pairs
import os, re, glob, csv

AA = os.path.expanduser(os.environ.get("AA",
      "~/krispmer/experiments/krispmer-expreminets/manuscript-experiments/human/july-13/specialized_targets"))
OUT = os.path.expanduser(os.environ.get("OUT", f"{AA}/alleleanalyzer_discriminability_human.tsv"))

CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
           "ref_multiple_v_single":"accept","present_v_absent":"accept"}
ENST = re.compile(r"ENST\d+")

def window_loci(cls):
    """distinct version-stripped ENST across the 50 target_*.fasta headers."""
    s=set()
    for f in glob.glob(f"{AA}/{cls}/target_*.fasta"):
        with open(f) as fh:
            for ln in fh:
                if ln.startswith(">"):
                    m=ENST.search(ln)
                    if m: s.add(m.group(0))
    return s

def aa_loci(cls):
    """distinct version-stripped ENST AA built >=1 guide pair for, + pair count."""
    f=f"{AA}/{cls}/alleleanalyzer.tsv"
    loci=set(); pairs=0
    if not os.path.exists(f): return loci,0
    with open(f) as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        for r in rd:
            pairs+=1
            m=ENST.search(r.get("locus","") or "")
            if m: loci.add(m.group(0))
    return loci,pairs

rows=[]
print(f"{'class':22s} {'truth':6s} {'win_loci':>8s} {'aa_loci':>7s} {'covered':>7s} {'pct':>5s} {'pairs':>6s}")
for cls in CLASSES:
    W=window_loci(cls); A,pairs=aa_loci(cls)
    covered=W & A
    n=len(W); nc=len(covered)
    pct = 100.0*nc/n if n else 0.0
    print(f"{cls:22s} {TRUTH[cls]:6s} {n:8d} {len(A):7d} {nc:7d} {pct:5.1f} {pairs:6d}")
    rows.append({"class":cls,"truth":TRUTH[cls],"n_window_loci":n,
                 "n_aa_loci":len(A),"n_covered":nc,"pct_covered":round(pct,2),
                 "n_guide_pairs":pairs})

with open(OUT,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}")
