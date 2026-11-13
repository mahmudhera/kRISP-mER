#!/bin/bash
# 06_run_crispor_n22.sh — CRISPOR on the FRESH per-class N22 targets (japonica reference)
# kRISP-meR reads the individual (N22); CRISPOR scores the SAME guides against the
# japonica reference genome -> individual-vs-reference head-to-head.
set -uo pipefail
CRISPOR_DIR=~/krispmer/new_datasets/feb-23-rice-n22/crispor/crisporWebsite
GENOME=rice_ref_annotated
FRESH=~/krispmer/new_datasets/exps/n22/fresh
TGT="$FRESH/targets"
OUT="$FRESH/crispor_class_scores"
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"
export MPLCONFIGDIR=/tmp/mpl_$$
cd "$CRISPOR_DIR"
for CLS in $CLASSES; do
  mkdir -p "$OUT/$CLS"
  for t in "$TGT/$CLS"/target_*.fasta; do
    [ -e "$t" ] || { echo "[warn] no targets in $CLS"; break; }
    name=$(basename "$t" .fasta)
    o="$OUT/$CLS/${name}_crispor.tsv"; off="$OUT/$CLS/${name}_off.tsv"
    [ -s "$o" ] && { echo "[skip] $CLS/$name"; continue; }
    echo "[run ] $CLS/$name"
    python crispor.py "$GENOME" "$t" "$o" -o "$off" --noEffScores \
      > "$OUT/$CLS/${name}.log" 2>&1 || echo "[FAIL] $CLS/$name"
  done
done
echo "[all done]"
