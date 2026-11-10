#!/usr/bin/env bash
# 01_build_n22_exons.sh
# ----------------------------------------------------------------------------
# STEP 1 of the clean N22 rebuild.
# Extract spliced transcripts from the OsN22RS2 Ensembl assembly+gff3 -> the
# ONE exon/transcript source used for ALL four classes' target windows.
# Headers normalized to bare Ensembl transcript IDs (OsN22_xxG..._xx) so that
# guides, SAM refs, gene lists and windows all share ONE namespace.
# (This is what broke before: buckets were assembly k-mers, but the old target
#  windows/SAMs mixed NCBI-RefSeq rna-XM_066... ids -> guides never matched.)
# ----------------------------------------------------------------------------
set -euo pipefail
FRESH=~/krispmer/new_datasets/exps/n22/fresh
mkdir -p "$FRESH"; cd "$FRESH"

N22DIR=~/krispmer/new_datasets/oryza_sativa_n22
N22G="$N22DIR/Oryza_sativa_n22.OsN22RS2.dna.toplevel.fa"
N22GFF="$N22DIR/feb-23/Oryza_sativa_n22.OsN22RS2.57.gff3"

[ -s "$N22G" ]   || { echo "ERROR: missing genome $N22G"   >&2; exit 1; }
[ -s "$N22GFF" ] || { echo "ERROR: missing gff3 $N22GFF"   >&2; exit 1; }

echo "[*] gffread: extracting spliced transcripts ..."
gffread -w n22_exons_raw.fa -g "$N22G" "$N22GFF"
# strip the "transcript:" prefix -> headers become OsN22_01G000010_01
sed 's/^>transcript:/>/' n22_exons_raw.fa > n22_exons.fa
rm -f n22_exons_raw.fa

echo "[*] building bowtie index over the transcript source ..."
bowtie-build --quiet n22_exons.fa n22_exons_idx

echo "transcripts : $(grep -c '^>' n22_exons.fa)"
echo "sample hdrs :"; grep -m3 '^>' n22_exons.fa
echo "[done] n22_exons.fa + bowtie index (n22_exons_idx.*) ready in $FRESH"
