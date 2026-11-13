# kRISP-meR class-analysis full pipeline (Stage A + Stage B)

Captured 2026-07-10. Built/run on human (chm13 vs grch38) in
`~/human_assemblies_kmer_count/July-10`. Goal: a single combined script,
to be adapted for **rice N22** at the end (see N22 ADAPTATION notes).

## Concept
Assembly-vs-assembly comparison. CASE = sample assembly, CONTROL = reference assembly.
Canonical PAM 23-mers (jellyfish -m 23 -C, PAM filter `^CC || GG$`). Counts = copy numbers.

## STAGE A — counting + classification (produces the 4 buckets)
1. `jellyfish count -m 23 -C -L 0 -U 999999 GENOME.fasta -o mer_counts.jf`
   `jellyfish dump -c mer_counts.jf > *_counts_readable`
   PAM filter: `awk '$1 ~ /^CC/ || $1 ~ /GG$/' *_counts_readable > *_counts_gRNAs.txt`
   `sort *_counts_gRNAs.txt > *_counts_gRNAs_sorted.txt`   # sorted = the hawk input (DATA lives here; unsorted may be empty)
2. hawk (krisp_hawk): reads case_sorted_files.txt + control_sorted_files.txt
   -> 4 buckets: present_v_absent.txt (case>0,ref=0), multiple_v_single.txt (case>1,ref=1),
      ref_present_v_absent.txt (case=0,ref>0), ref_multiple_v_single.txt (case=1,ref>1)

## JULY-10 SETUP (isolate, don't clobber old outputs)
```
cd ~/human_assemblies_kmer_count ; mkdir -p July-10 ; cd July-10
ln -sf ../exons.fa . ; ln -sf ../exons.*.bt2 .                 # read-only inputs -> symlink
ln -sf ../chm13_counts_gRNAs_sorted.txt  chm13_counts_gRNAs.txt   # NB: point at the *_sorted* files (unsorted are EMPTY)
ln -sf ../grch38_counts_gRNAs_sorted.txt grch38_counts_gRNAs.txt
ln -sf ../present_v_absent.txt . ; ln -sf ../multiple_v_single.txt .
ln -sf ../ref_present_v_absent.txt . ; ln -sf ../ref_multiple_v_single.txt .
# scripts -> copy (may edit)
```

## STAGE B step 1 — per-bucket: FASTA -> bowtie -> mismatch filter -> gene IDs  (runStageB.sh)
```
#!/bin/bash
# runStageB.sh BUCKET [all|best]
set -uo pipefail
BUCKET="$1"; MODE="${2:-best}"; REF=exons
BT2=/home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2
awk '{printf ">%s_%s_%s\n%s\n", NR,$2,$3,$1}' "${BUCKET}.txt" > "${BUCKET}.fasta"
AFLAG=""; [ "$MODE" = "all" ] && AFLAG="-a"
$BT2 -p 16 $AFLAG -x "$REF" -f "${BUCKET}.fasta" --no-unal -S "${BUCKET}_${REF}.sam"
grep 'XM:i:0' "${BUCKET}_${REF}.sam" > "${BUCKET}_${REF}_max_0.sam" || true
grep 'XM:i:0\|XM:i:1\|XM:i:2\|XM:i:3' "${BUCKET}_${REF}.sam" > "${BUCKET}_${REF}_max_3.sam" || true
for mm in 0 3; do
  awk '{print $3}' "${BUCKET}_${REF}_max_${mm}.sam" | awk -F'_' '{print $3}' \
    | sort | uniq > "geneIDs_${BUCKET}_max_${mm}_uniq.txt"     # <-- GENCODE header (human). N22 differs (see below)
done
```
Run: multiple_v_single & ref_multiple_v_single with `all` (-a); present_v_absent & ref_present_v_absent with `best`.
(present_v_absent = 59.7M guides -> background.)

## STAGE B step 2 — single-copy "unique-guide" transcript sets  (runUnique.sh)
```
#!/bin/bash
# runUnique.sh ASM   (chm13 | grch38 ; for N22: n22 | ref)
set -uo pipefail
ASM="$1"; REF=exons
BT2=/home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2
awk '($1 ~ /^CC/ || $1 ~ /GG$/) && $2 == 1 {print}' "${ASM}_counts_gRNAs.txt" > "${ASM}_unique_gRNAs.txt"
awk '{printf ">%s_%s\n%s\n", NR,$2,$1}' "${ASM}_unique_gRNAs.txt" > "${ASM}_unique_gRNAs.fasta"
$BT2 -p 16 -x "$REF" -f "${ASM}_unique_gRNAs.fasta" --no-unal -S "${ASM}_unique_gRNAs_${REF}.sam"
grep 'XM:i:0' "${ASM}_unique_gRNAs_${REF}.sam" > "${ASM}_unique_gRNAs_${REF}_max_0.sam" || true
awk '{print $3}' "${ASM}_unique_gRNAs_${REF}_max_0.sam" | awk -F'_' '{print $3}' | sort | uniq > "U_${ASM}.txt"
```
Run sequentially (each -p 16): `bash runUnique.sh chm13; bash runUnique.sh grch38` (~18 min total).

## STAGE B step 3 — comm cross-analysis (final target-gene lists)
```
export LC_ALL=C
comm -12 <(sort geneIDs_multiple_v_single_max_0_uniq.txt)     <(sort U_chm13.txt)  > genes_multiple_v_single_analysis.txt
comm -23 <(sort geneIDs_present_v_absent_max_3_uniq.txt)      <(sort U_grch38.txt) > genes_present_v_absent_analysis.txt
comm -12 <(sort geneIDs_ref_multiple_v_single_max_0_uniq.txt) <(sort U_grch38.txt) > genes_ref_multiple_v_single_analysis.txt
comm -23 <(sort geneIDs_ref_present_v_absent_max_0_uniq.txt)  <(sort U_chm13.txt)  > genes_ref_present_v_absent_analysis.txt
```
Logic: bucket-side mismatch = max_0 if the guides are PRESENT in the reference (multiple_v_single, ref_*),
max_3 only for present_v_absent (CHM13-only, no exact hg38 match). comm -12 = intersect (the "multiple" classes),
comm -23 = subtract (the "presence" classes). ref_* pairings are the symmetric (swap sample<->ref) inference.

## HUMAN RESULTS (this run)
Stage B per-bucket exon-align: multiple_v_single 34.75% (genes max0=29,889 / max3=56,180);
ref_multiple_v_single 19.07% (6,463 / 32,209); ref_present_v_absent 6.66% (27,598 / 44,061);
present_v_absent 7.13% (24,813 / 76,255).
Unique sets (fresh, best-hit): U_chm13=93,166 ; U_grch38=77,924  (old -a versions: 94,517 / 81,214).
comm results (final target-gene lists):
  genes_multiple_v_single_analysis.txt      = 29,565   (of 29,889 pre-comm; comm -12 vs U_chm13)
  genes_present_v_absent_analysis.txt       = 16,208   (of 76,255 max_3; comm -23 vs U_grch38)
  genes_ref_multiple_v_single_analysis.txt  =  6,129   (of 6,463; comm -12 vs U_grch38)
  genes_ref_present_v_absent_analysis.txt   =    214   (of 27,598; comm -23 vs U_chm13)
HEADLINE: present_v_absent (16,208) vs ref_present_v_absent (214) = ~75x asymmetry
  -> CHM13 opens unique guide access to ~16k genes GRCh38 design would miss; reverse only 214.
  (The comm -23 subtraction drives the result; comm -12 intersections barely filter because
   almost every gene has some single-copy guide.)

## N22 ADAPTATION (for the combined script the user will request)
- Labels: sample = **n22** (not chm13), reference = **ref** (not grch38).
- Counts files: case=`lexo_mer_counts_PAM_sorted_n22.txt`, control=`lexo_mer_counts_PAM_sorted.txt`
  (in exps/n22). For runUnique, rename to n22_counts_gRNAs.txt / ref_counts_gRNAs.txt (symlink the sorted data).
- **Exon header parse differs**: human = GENCODE `hg38_wgEncodeGencodeBasicV38_ENST..._...` -> `awk -F'_' '{print $3}'`.
  N22 = RefSeq `rna-XM_...` -> `sed 's/^rna-//'` (see n22/getGeneList_present_v_absent_all.sh).
  So the gene-ID extraction line must switch by dataset.
- exons index/fa: use the n22 dir's exons.* (rice exons).
- comm pairings map: chm13->n22, grch38->ref (U_n22, U_ref).
- Everything else identical.
