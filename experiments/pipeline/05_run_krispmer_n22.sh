#!/usr/bin/env bash
# 05_run_krispmer_n22.sh — kRISP-meR over the FRESH per-class targets.
# Reuses the real-N22 reads DB (jf/histo) + already-fitted EM (lambda=29.40),
# so every target hits the -ein branch (no re-fitting). max_hd=2.
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate krispmer

FRESH=~/krispmer/new_datasets/exps/n22/fresh
READS_DIR=~/krispmer/new_datasets/n22-real
READS=$(ls "$READS_DIR"/n22_reads.fastq "$READS_DIR"/*.fastq 2>/dev/null | head -1)
JF="$READS_DIR/mer_counts.jf"
HISTO="$READS_DIR/n22_histo"
EM="$READS_DIR/em_n22_real.txt"
MAX_HD=2
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"

cd "$FRESH"
[ -s "$JF" ]    || { echo "ERROR: missing $JF"; exit 1; }
[ -s "$HISTO" ] || { echo "ERROR: missing $HISTO"; exit 1; }
[ -s "$EM" ]    || { echo "ERROR: missing fitted EM $EM"; exit 1; }
# reads file optional here: -J/-H/-ein are all precomputed. warn only.
[ -n "${READS:-}" ] && [ -s "$READS" ] || echo "WARN: no reads fastq (ok: JF/HISTO/EM precomputed)"
echo "jf=$JF"; echo "histo=$HISTO"; echo "em=$EM (lambda=29.40, -ein)"

for CLS in $CLASSES; do
  mkdir -p "scores/$CLS" "logs/$CLS"
  for t in targets/$CLS/target_*.fasta; do
    [ -e "$t" ] || { echo "[warn] no targets in $CLS"; break; }
    name=$(basename "$t" .fasta)
    out="scores/$CLS/${name}.csv"
    [ -s "$out" ] && { echo "[skip] $CLS/$name"; continue; }
    echo "[run ] $CLS/$name  -ein"
    krispmer "${READS:-$t}" "$t" "$out" "$MAX_HD" \
      -J "$JF" -H "$HISTO" -ein "$EM" \
      -l "logs/$CLS/${name}.log" || echo "[FAIL] $CLS/$name"
  done
done
echo "[all done]"
