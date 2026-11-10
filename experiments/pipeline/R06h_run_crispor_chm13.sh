#!/usr/bin/env bash
# R06h_run_crispor_chm13.sh — CRISPOR on the INDIVIDUAL-anchored (CHM13) HUMAN
# targets. Paths taken from the ORIGINAL human run_crispor_classes.sh:
#   CRISPOR = feb-23-rice-n22/crispor/crisporWebsite/crispor.py, GENOME = hg38.
# CRISPOR still scores against hg38 (the reference) -- only the target SOURCE is
# CHM13-anchored now. Reads ind_targets/, writes ind_crispor_class_scores/.
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate krispmer
export MPLCONFIGDIR=/tmp/mpl-$USER; mkdir -p "$MPLCONFIGDIR"

CB=/home/atif/krispmer/new_datasets/feb-23-rice-n22/crispor/crisporWebsite/crispor.py
GENOME=hg38
H=~/human_assemblies_kmer_count/july-10
TGT="$H/ind_targets"
OUT="$H/ind_crispor_class_scores"
CLASSES="multiple_v_single present_v_absent ref_multiple_v_single ref_present_v_absent"
cd "$H"
for CLS in $CLASSES; do
  mkdir -p "$OUT/$CLS"
  for f in "$TGT/$CLS"/target_*.fasta "$TGT/$CLS"/target_*.fa; do
    [ -e "$f" ] || continue
    name=$(basename "$f" | sed 's/\.\(fasta\|fa\)$//')
    out="$OUT/$CLS/${name}_crispor.tsv"
    [ -s "$out" ] && { echo "[skip] $CLS/$name"; continue; }
    echo "[run ] $CLS/$name"
    python "$CB" "$GENOME" "$f" "$out" -o "$OUT/$CLS/${name}_off.tsv" --noEffScores 2>/dev/null \
      || echo "[FAIL] $CLS/$name"
  done
done
echo "[all done]"
