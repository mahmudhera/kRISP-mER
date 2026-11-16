#!/usr/bin/env bash
# R01h_build_chm13_exons.sh  — INDIVIDUAL-anchored step 1 for HUMAN.
# Build CHM13 exons from the CHM13v2.0 assembly + T2T CAT/Liftoff GENCODE GFF3,
# mirroring rice R01 (N22 assembly + N22 GFF). This is the genuine individual
# (CHM13) gene-model source, the human analogue of N22's OsN22RS2 exons.
#
# WRINKLE: the assembly FASTA uses GenBank accession seqids (>CP068277.2 ...
# "chromosome 1"), but the CAT/Liftoff GFF uses chr1..chrY. So we RENAME the
# FASTA seqids CP->chr (deterministically, from each header's own "chromosome N"
# text) before gffread. Output CHM13 exon FASTA + bowtie index under HW/.
set -euo pipefail

HW=~/human_assemblies_kmer_count/july-10          # where the human run lives
ASMDIR=~/human_assemblies_kmer_count/chm13_assembly
ASM="$ASMDIR/cmh13v2.fasta"
GFF="$ASMDIR/chm13v2.0_CAT_Liftoff.gff3"          # chr-named CAT/Liftoff GENCODE
cd "$HW"

[ -s "$ASM" ] || { echo "ERROR: missing CHM13 assembly $ASM" >&2; exit 1; }
[ -s "$GFF" ] || { echo "ERROR: missing CHM13 GFF $GFF"       >&2; exit 1; }

echo "[*] renaming assembly seqids CP-accession -> chr (from header text) ..."
# For each '>CP... chromosome N' header emit '>chrN'; X/Y handled; drop the rest.
# mito in CHM13v2 is 'mitochondrion' -> chrM.
awk '
  /^>/ {
    hdr=$0
    # extract the token after "chromosome "
    if (match(hdr, /chromosome ([0-9XY]+)/, m)) { print ">chr" m[1]; next }
    if (hdr ~ /mitochondrion/)                  { print ">chrM";     next }
    # unknown contig: keep a safe unique name so gffread never errors on it
    n++; print ">unplaced_" n; next
  }
  { print }
' "$ASM" > chm13v2_chrnamed.fasta

echo "[*] renamed contigs present:"; grep '^>' chm13v2_chrnamed.fasta | head -26

echo "[*] gffread: extracting spliced transcripts (chr-named) ..."
gffread -w chm13_exons_raw.fa -g chm13v2_chrnamed.fasta "$GFF"
# normalize headers: keep the transcript id token only (strip 'transcript:'/'rna-')
sed -e 's/^>rna-/>/' -e 's/^>transcript:/>/' chm13_exons_raw.fa > chm13_exons.fa
rm -f chm13_exons_raw.fa

echo "[*] building bowtie index over the CHM13 transcript source ..."
rm -f chm13_exons_idx*.ebwt          # clear any partial index from an interrupted run
bowtie-build --quiet chm13_exons.fa chm13_exons_idx

echo "transcripts : $(grep -c '^>' chm13_exons.fa)"
echo "sample hdrs :"; grep -m3 '^>' chm13_exons.fa
echo "[done] chm13_exons.fa + bowtie index (chm13_exons_idx.*) ready in $HW"
