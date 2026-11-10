#!/usr/bin/env bash
# 05b_run_krispmer_n22_bothstrand.sh — kRISP-meR over the FRESH per-class targets,
# BOTH STRANDS (adds -n / --consider_negative). The original 05_* run had no -n,
# so it designed only + strand guides; ~54% of present_v_absent class guides were
# "not-generated" simply because their NGG PAM sits on the - strand (VERIFIED:
# not-generated set is 129/129 minus-strand-PAM guides; scratchpad/verify_strand_mechanism.py).
# This run realizes both strands so the head-to-head "not-generated" column
# reflects only genuine read-support drops (value1<=0), apples-to-apples with
# CRISPOR/GuideScan (which design both strands by default).
#
# Writes to scores_bothstrand/ (does NOT clobber the + strand-only scores/).
# Reuses the real-N22 reads DB (jf/histo) + fitted EM (lambda=29.40) via -ein.
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
[ -n "${READS:-}" ] && [ -s "$READS" ] || echo "WARN: no reads fastq (ok: JF/HISTO/EM precomputed)"
echo "jf=$JF"; echo "histo=$HISTO"; echo "em=$EM (lambda=29.40, -ein)  MODE=both-strand (-n)"

for CLS in $CLASSES; do
  mkdir -p "ref_scores_bothstrand/$CLS" "ref_logs_bothstrand/$CLS"
  for t in ref_targets/$CLS/target_*.fasta; do
    [ -e "$t" ] || { echo "[warn] no targets in $CLS"; break; }
    name=$(basename "$t" .fasta)
    out="ref_scores_bothstrand/$CLS/${name}.csv"
    [ -s "$out" ] && { echo "[skip] $CLS/$name"; continue; }
    echo "[run ] $CLS/$name  -ein -n"
    krispmer "${READS:-$t}" "$t" "$out" "$MAX_HD" \
      -J "$JF" -H "$HISTO" -ein "$EM" -n \
      -l "ref_logs_bothstrand/$CLS/${name}.log" || echo "[FAIL] $CLS/$name"
  done
done
echo "[all done] both-strand scores in scores_bothstrand/"
