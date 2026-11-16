#!/usr/bin/env python3
# AA3_variant_reality_human.py  --  run ON THE VM.
#
# QUESTION (Navid's framing): AlleleAnalyzer is reference-VCF-driven. At our target
# loci, does the allele AA proposes to target ACTUALLY EXIST in the individual
# (CHM13)?  We have the individual as ground truth = CHM13 k=23 jellyfish counts.
# So for every AA-proposed allele-specific guide we ask: is its 23-mer (protospacer
# + NGG PAM) present in CHM13?  This measures reference-blindness directly.
#
# AA emits 20-mer protospacers (gRNA_ref / gRNA_alt), no PAM, but gives
# chrom/start/stop/strand -> reconstruct the 23-mer from hg38.fa:
#   positive strand: 23mer = protospacer + hg38[stop : stop+3]      (NGG)
#   negative strand: guide is on the minus strand; the 20-mer as written is already
#                    the protospacer 5'->3'. PAM (NGG on the guide's strand) sits at
#                    hg38[start-3 : start] on the plus strand as CCN -> revcomp = NGG.
#                    23mer = protospacer + revcomp(hg38[start-3:start]).
# jf is canonical (-C) so query strand doesn't matter; we still build the biological
# 23-mer and let jellyfish canonicalize.
#
# For gRNA_alt we must substitute the alt allele into the reconstructed sequence at
# the variant position (AA already did this in the 20-mer; the PAM base is genomic
# and unchanged unless the variant is IN the PAM window -- rare; handled by taking
# the PAM from genome, which reflects ref. That's fine: PAM identity NGG only needs
# positions 2-3 = GG, allele rarely changes that; we record if it does).
#
# Output: alleleanalyzer_variant_reality_human.tsv  (per class summary)
#   class truth n_rows alt_in_chm13 alt_absent ref_in_grch38 alt_frac_real
# plus a per-row dump alleleanalyzer_variant_reality_rows_human.tsv
import os, csv, subprocess, tempfile

AA   = os.path.expanduser(os.environ.get("AA",
        "~/krispmer/experiments/krispmer-expreminets/manuscript-experiments/human/july-13/specialized_targets"))
HG38 = os.path.expanduser(os.environ.get("HG38",
        "~/krispmer/new_datasets/feb-23-rice-n22/crispor/crisporWebsite/genomes/hg38/hg38.fa"))
CHM13_JF = os.path.expanduser(os.environ.get("CHM13_JF",
        "~/human_assemblies_kmer_count/chm13_assembly/chm13_counts_binary.jf"))
GRCH38_JF = os.path.expanduser(os.environ.get("GRCH38_JF",
        "~/human_assemblies_kmer_count/grch38.jf"))
OUTDIR = os.path.expanduser(os.environ.get("OUTDIR", AA))

CLASSES = ["multiple_v_single","ref_present_v_absent","ref_multiple_v_single","present_v_absent"]
TRUTH   = {"multiple_v_single":"reject","ref_present_v_absent":"reject",
           "ref_multiple_v_single":"accept","present_v_absent":"accept"}
# AA marks "no guide buildable on this allele" with a 20-base HOMOPOLYMER sentinel.
# It is NOT always poly-G: ref-field sentinels are poly-G; alt-field sentinels are
# poly-C (minus strand) and also poly-A/poly-T. So treat ANY 20-mer of a single base
# as a sentinel. (Verified: {poly-G ref, poly-C/A/T alt}; a few rows sentinel BOTH.)
_HOMO = {b*20 for b in "ACGT"}
def is_sentinel(s): return s in _HOMO

def rc(s): return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]

# ---- pysam if available, else samtools faidx for PAM base retrieval ----------
def faidx(region):
    """region like 'NC_000001.11:225984104-225984106' -> uppercase seq (samtools)."""
    r=subprocess.run(["samtools","faidx",HG38,region],capture_output=True,text=True)
    return "".join(l.strip() for l in r.stdout.splitlines() if not l.startswith(">")).upper()

def build23(chrom,start,stop,strand,protospacer):
    """protospacer is the 20-mer as AA wrote it. Return biological 23-mer w/ PAM."""
    if strand=="positive":
        pam=faidx(f"{chrom}:{stop+1}-{stop+3}")          # 1-based, base after stop
        return protospacer+pam
    else:
        ccn=faidx(f"{chrom}:{start-2}-{start}")           # 3 bases before start (plus strand)
        return protospacer+rc(ccn)

def jf_query(jf, seqs):
    if not seqs: return {}
    with tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False) as fh:
        tmp=fh.name
        for i,s in enumerate(sorted(seqs)): fh.write(f">v{i}\n{s}\n")
    r=subprocess.run(["jellyfish","query",jf,"-s",tmp],capture_output=True,text=True); os.unlink(tmp)
    out={}
    for ln in r.stdout.splitlines():
        p=ln.split()
        if len(p)==2: out[p[0]]=int(p[1])
    return out

summary=[]; allrows=[]
for cls in CLASSES:
    f=f"{AA}/{cls}/alleleanalyzer.tsv"
    if not os.path.exists(f):
        print(f"[skip] {cls}: no alleleanalyzer.tsv"); continue
    rows=list(csv.DictReader(open(f),delimiter="\t"))
    # AA emits a poly-G sentinel (GGGGGGGGGGGGGGGGGGGG) in the ref OR alt column when
    # NO guide is buildable on THAT allele -- i.e. that allele has no PAM (the variant
    # either sits outside the guide window or disrupts/creates the PAM). These are NOT
    # bad guides: a poly-G ref + real alt is an allele-SPECIFIC-by-PAM-creation guide
    # (only the alt allele carries a cut site). So we do NOT drop the row; we skip the
    # SENTINEL FIELD and only reconstruct/query the real guide(s). We also tally the
    # PAM-only categories separately (they are a legitimate AA output class).
    # 'both-polyG' never occurs (verified) but is guarded anyway.
    recon=[]
    for r in rows:
        chrom=r["chrom"]; strand=r["strand"]
        gref=r["gRNA_ref"].upper(); galt=r["gRNA_alt"].upper()
        ref_sentinel = is_sentinel(gref); alt_sentinel = is_sentinel(galt)
        a23=r23=None
        try:
            start=int(r["start"]); stop=int(r["stop"])
            if not alt_sentinel: a23=build23(chrom,start,stop,strand,galt)
            if not ref_sentinel: r23=build23(chrom,start,stop,strand,gref)
        except Exception:
            a23=r23=None
        recon.append((r,a23,r23,ref_sentinel,alt_sentinel))
    alt_set={a for _,a,_,_,_ in recon if a}; ref_set={x for _,_,x,_,_ in recon if x}
    chm13=jf_query(CHM13_JF, alt_set|ref_set)
    grch38=jf_query(GRCH38_JF, alt_set|ref_set)
    def cnt(seq,tbl):
        if not seq: return 0
        return tbl.get(seq,0)+tbl.get(rc(seq),0)   # jf canonical, but be safe
    def is_snp(r): return len(r["ref"])==1 and len(r["alt"])==1
    # --- SNP-only full 2x2: {ref,alt} guide x {CHM13,GRCh38} presence -----------
    # (indels can't be k-mer-checked by appending a ref PAM -> tallied separately)
    # denominators are per-guide (a row can contribute a ref guide and/or an alt guide,
    # whichever is non-sentinel). Controls: ref-in-GRCh38 ~ alt-not-in-GRCh38 validate
    # reconstruction; the finding is ref-in-CHM13 (how often AA's reference allele
    # actually exists in the individual).
    c=dict(n_alt=0,alt_in_chm13=0,alt_in_grch38=0,
           n_ref=0,ref_in_chm13=0,ref_in_grch38=0)
    indel_rows=0; both_sentinel=0
    pam_alt_only=0; pam_ref_only=0; both_guides=0
    for r,a23,r23,ref_sent,alt_sent in recon:
        if ref_sent and alt_sent: both_sentinel+=1; continue
        if not is_snp(r): indel_rows+=1; continue          # SNP-only analysis
        if ref_sent and not alt_sent: pam_alt_only+=1
        elif alt_sent and not ref_sent: pam_ref_only+=1
        else: both_guides+=1
        ac=aG=rc13=rG=None
        if a23 is not None:
            ac=cnt(a23,chm13); aG=cnt(a23,grch38)
            c["n_alt"]+=1; c["alt_in_chm13"]+= (ac>0); c["alt_in_grch38"]+= (aG>0)
        if r23 is not None:
            rc13=cnt(r23,chm13); rG=cnt(r23,grch38)
            c["n_ref"]+=1; c["ref_in_chm13"]+= (rc13>0); c["ref_in_grch38"]+= (rG>0)
        allrows.append({"class":cls,"locus":r["locus"],"guide_id":r["guide_id"],
                        "strand":r["strand"],"ref":r["ref"],"alt":r["alt"],
                        "variant_position_in_guide":r["variant_position_in_guide"],
                        "gRNA_ref_23":r23 or "","ref_in_chm13":"" if rc13 is None else int(rc13>0),
                        "ref_in_grch38":"" if rG is None else int(rG>0),
                        "gRNA_alt_23":a23 or "","alt_in_chm13":"" if ac is None else int(ac>0),
                        "alt_in_grch38":"" if aG is None else int(aG>0)})
    def pc(n,d): return round(100*n/d,1) if d else 0.0
    print(f"### {cls}  (truth {TRUTH[cls]})   SNP guides: ref={c['n_ref']}  alt={c['n_alt']}   "
          f"[indel rows set aside: {indel_rows}; both-sentinel dropped: {both_sentinel}]")
    print(f"    REF allele guide  ->  in GRCh38 {pc(c['ref_in_grch38'],c['n_ref']):5.1f}% (control ok)   "
          f"in CHM13 {pc(c['ref_in_chm13'],c['n_ref']):5.1f}%  <== reference allele present in individual")
    print(f"    ALT allele guide  ->  in GRCh38 {pc(c['alt_in_grch38'],c['n_alt']):5.1f}% (control: low)  "
          f"in CHM13 {pc(c['alt_in_chm13'],c['n_alt']):5.1f}%  (individual carries alt)")
    print(f"    PAM-specific: alt-only(ref=homopolymer)={pam_alt_only}  ref-only={pam_ref_only}  both-allele={both_guides}\n")
    summary.append({"class":cls,"truth":TRUTH[cls],
                    "n_snp_ref_guides":c["n_ref"],"ref_in_grch38_pct":pc(c["ref_in_grch38"],c["n_ref"]),
                    "ref_in_chm13_pct":pc(c["ref_in_chm13"],c["n_ref"]),
                    "n_snp_alt_guides":c["n_alt"],"alt_in_grch38_pct":pc(c["alt_in_grch38"],c["n_alt"]),
                    "alt_in_chm13_pct":pc(c["alt_in_chm13"],c["n_alt"]),
                    "ref_in_chm13_n":c["ref_in_chm13"],"alt_in_chm13_n":c["alt_in_chm13"],
                    "indel_rows_set_aside":indel_rows,"both_sentinel_dropped":both_sentinel,
                    "pam_alt_only":pam_alt_only,"pam_ref_only":pam_ref_only,"both_allele_guides":both_guides})

s_out=f"{OUTDIR}/alleleanalyzer_variant_reality_human.tsv"
with open(s_out,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
r_out=f"{OUTDIR}/alleleanalyzer_variant_reality_rows_human.tsv"
with open(r_out,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)
print(f"\nwrote {s_out}\nwrote {r_out}")
