#!/usr/bin/env bash
# 05b_run_krispmer_human_bothstrand.sh — kRISP-meR over the HUMAN per-class targets,
# BOTH STRANDS (adds -n / --consider_negative). Mirrors the N22 both-strand run:
# the original human scoring run had no -n, so ~half the class guides (minus-strand
# PAM) were "not-generated" as a strand-coverage artifact, not a rejection.
# Writes scores_bothstrand/ (does NOT clobber the +strand-only scores/).
#
# Reuses the human reads DB (jf/histo) + fitted EM via -ein. This script AUTO-
# DISCOVERS those inputs under H and PRINTS them first -- verify they match the
# files the original human run used before letting it proceed.
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate krispmer

H=~/human_assemblies_kmer_count/july-10
MAX_HD=2
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"

# --- auto-discover the reads DB / histo / EM (edit here if the guess is wrong) ---
JF=$(ls "$H"/mer_counts.jf "$H"/*.jf 2>/dev/null | head -1)
HISTO=$(ls "$H"/*histo* 2>/dev/null | head -1)
EM=$(ls "$H"/em_*.txt "$H"/*em*.txt 2>/dev/null | head -1)
READS=$(ls "$H"/*.fastq "$H"/*.fq 2>/dev/null | head -1)

echo "=== discovered inputs (verify before proceeding) ==="
echo "  H     = $H"
echo "  JF    = ${JF:-<none>}"
echo "  HISTO = ${HISTO:-<none>}"
echo "  EM    = ${EM:-<none>}"
echo "  READS = ${READS:-<none>  (ok if JF/HISTO/EM present)}"
echo "===================================================="
[ -s "${JF:-}" ]    || { echo "ERROR: no jf found under $H — set JF= manually"; exit 1; }
[ -s "${HISTO:-}" ] || { echo "ERROR: no histo found under $H — set HISTO= manually"; exit 1; }
[ -s "${EM:-}" ]    || { echo "ERROR: no EM found under $H — set EM= manually"; exit 1; }

cd "$H"
for CLS in $CLASSES; do
  mkdir -p "ind_scores_bothstrand/$CLS" "ind_logs_bothstrand/$CLS"
  for t in ind_targets/$CLS/target_*.fasta ind_targets/$CLS/target_*.fa; do
    [ -e "$t" ] || continue
    name=$(basename "$t" | sed 's/\.\(fasta\|fa\)$//')
    out="ind_scores_bothstrand/$CLS/${name}.csv"
    [ -s "$out" ] && { echo "[skip] $CLS/$name"; continue; }
    echo "[run ] $CLS/$name  -ein -n"
    krispmer "${READS:-$t}" "$t" "$out" "$MAX_HD" \
      -J "$JF" -H "$HISTO" -ein "$EM" -n \
      -l "ind_logs_bothstrand/$CLS/${name}.log" || echo "[FAIL] $CLS/$name"
  done
done
echo "[all done] human both-strand scores in $H/scores_bothstrand/"
