#!/usr/bin/env bash
# R01_build_japonica_exons.sh
# ----------------------------------------------------------------------------
# REFERENCE-anchored variant of step 1: build the exon/transcript source from
# the JAPONICA reference (GCF_034140825.1) instead of the N22 individual.
# Meeting to-do (2026-07): re-run the per-class head-to-head with targets built
# from the REFERENCE genome (japonica for rice, GRCh38 for human), mirroring the
# individual-anchored run. Output lives under fresh/ref_targets/ so it never
# clobbers the individual-anchored n22_exons.fa / targets/.
# ----------------------------------------------------------------------------
set -euo pipefail
FRESH=~/krispmer/new_datasets/exps/n22/fresh
mkdir -p "$FRESH"; cd "$FRESH"

JAPDIR=~/krispmer/new_datasets/rice_ref_annotated
JAPG="$JAPDIR/GCF_034140825.1_ASM3414082v1_genomic.fna"
JAPGFF="$JAPDIR/genomic_clean.gff"

[ -s "$JAPG" ]   || { echo "ERROR: missing japonica genome $JAPG" >&2; exit 1; }
[ -s "$JAPGFF" ] || { echo "ERROR: missing japonica gff $JAPGFF"  >&2; exit 1; }

# Keep ONLY the 12 nuclear chromosomes (NC_089035..NC_089046). The organellar
# contigs -- NC_011033.1 (mito), NC_001320.1 (plastid), NC_001751.1 (mito
# plasmid) -- have circular coordinates that crash gffread, and we don't want
# organellar guides in a nuclear off-target comparison anyway. Allowlist the
# nuclear accessions (keep GFF header '#' lines too).
echo "[*] filtering to the 12 nuclear chromosomes ..."
grep -E "^#|^NC_0890(3[5-9]|4[0-6])\.1\b" "$JAPGFF" > jap_nuclear.gff3
JAPGFF_USE=jap_nuclear.gff3

echo "[*] gffread: extracting spliced transcripts from japonica reference ..."
gffread -w jap_exons_raw.fa -g "$JAPG" "$JAPGFF_USE"
# normalize headers: strip a leading "rna-"/"transcript:" if present so the
# namespace is consistent within this reference-anchored run.
sed -e 's/^>rna-/>/' -e 's/^>transcript:/>/' jap_exons_raw.fa > jap_exons.fa
rm -f jap_exons_raw.fa

echo "[*] building bowtie index over the japonica transcript source ..."
bowtie-build --quiet jap_exons.fa jap_exons_idx

echo "transcripts : $(grep -c '^>' jap_exons.fa)"
echo "sample hdrs :"; grep -m3 '^>' jap_exons.fa
echo "[done] jap_exons.fa + bowtie index (jap_exons_idx.*) ready in $FRESH"
