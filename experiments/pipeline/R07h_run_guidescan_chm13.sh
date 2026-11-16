#!/usr/bin/env bash
# R07h_run_guidescan_chm13.sh — GuideScan2 on the INDIVIDUAL-anchored (CHM13)
# HUMAN targets. Paths from the ORIGINAL human run_guidescan_classes.sh:
#   FEB = ~/krispmer/experiments/krispmer-expreminets/manuscript-experiments/human/feb-20
#   GEN = $FEB/generate_kmers.py, IDX = $FEB/human  (GRCh38 GuideScan index).
# GuideScan still scores against the GRCh38 index (reference) -- only the target
# SOURCE is CHM13. Reads ind_targets/, writes ind_guidescan_class_scores/.
# Matches the original: -n 16 -m 2 -a NGG --mode complete, id prefix CLS::name::
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate krispmer

FEB=~/krispmer/experiments/krispmer-expreminets/manuscript-experiments/human/feb-20
GEN="$FEB/generate_kmers.py"
IDX="$FEB/human"
H=~/human_assemblies_kmer_count/july-10
TGT="$H/ind_targets"
OUT="$H/ind_guidescan_class_scores"
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"
mkdir -p "$OUT/kmers" "$OUT/out"
cd "$H"
for CLS in $CLASSES; do
  KM="$OUT/kmers/${CLS}.csv"
  echo "id,sequence,pam,chromosome,position,sense" > "$KM"
  for f in "$TGT/$CLS"/target_*.fasta "$TGT/$CLS"/target_*.fa; do
    [ -e "$f" ] || continue
    name=$(basename "$f" | sed 's/\.\(fasta\|fa\)$//')
    python "$GEN" --pam NGG --prefix "${CLS}::${name}::" "$f" 2>/dev/null | tail -n +2 >> "$KM"
  done
  echo "[enumerate] $CLS  ($(($(wc -l < "$KM")-1)) guides)"
  guidescan enumerate "$IDX" -f "$KM" -n 16 -m 2 -a NGG \
    --format csv --mode complete -o "$OUT/out/${CLS}.csv"
done
echo "[all done]"
