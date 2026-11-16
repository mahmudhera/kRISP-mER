#!/bin/bash
# 07_run_guidescan_n22.sh — GuideScan2 on the FRESH per-class N22 targets (japonica ref)
# Reuses the japonica-reference GS index built earlier (july-13/gs/rice_ref).
set -uo pipefail
GK=~/krispmer/new_datasets/rice_ref_annotated/generate_kmers.py
IDX=~/krispmer/new_datasets/exps/n22/july-13/gs/rice_ref      # genome index, reusable
FRESH=~/krispmer/new_datasets/exps/n22/fresh
TGT="$FRESH/ref_targets"
KMDIR="$FRESH/ref_guidescan/kmers"
OUT="$FRESH/ref_guidescan_class_scores/out"
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"
mkdir -p "$KMDIR" "$OUT"
for CLS in $CLASSES; do
  KM="$KMDIR/${CLS}.csv"; first=1
  for t in "$TGT/$CLS"/target_*.fasta; do
    [ -e "$t" ] || { echo "[warn] no targets in $CLS"; break; }
    base=$(basename "$t" .fasta)
    if [ $first = 1 ]; then python "$GK" --pam NGG --prefix "${base}-" "$t" > "$KM"; first=0
    else python "$GK" --pam NGG --prefix "${base}-" "$t" | tail -n +2 >> "$KM"; fi
  done
  echo "[$CLS] $(( $(wc -l < "$KM") - 1 )) guides -> enumerating"
  guidescan enumerate "$IDX" -f "$KM" -m 2 -a NGG --format csv --mode complete \
    -o "$OUT/${CLS}_ont.csv" 2> "$OUT/${CLS}.log"
  echo "[$CLS] done -> $(wc -l < "$OUT/${CLS}_ont.csv") rows"
done
echo "[all done]"
