#!/usr/bin/env bash
# 00_build_and_verify.sh — clean N22 rebuild, steps 1-3 (STOPS at the gate).
# Run this first. Only if the GATE says PASS do you run 05_run_krispmer_n22.sh.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate krispmer
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "############ STEP 1: build n22_exons.fa + bowtie index ############"
bash "$HERE/01_build_n22_exons.sh"

echo; echo "############ STEP 2: gene lists + centered target windows ############"
for CLS in present_v_absent ref_present_v_absent multiple_v_single ref_multiple_v_single; do
  python3 "$HERE/02_build_targets_n22.py" "$CLS" 50
done

echo; echo "############ STEP 3: GATE — targets contain their own class guides? ############"
python3 "$HERE/03_verify_targets_n22.py"

echo
echo ">>> If the GATE says PASS, run:  bash $HERE/05_run_krispmer_n22.sh"
echo ">>> then:                        python3 $HERE/04_class_scoring_n22.py"
