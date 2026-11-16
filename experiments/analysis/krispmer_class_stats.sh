#!/usr/bin/env bash
# ============================================================================
# krispmer_class_stats.sh   --   prints to terminal (no file written)
# ----------------------------------------------------------------------------
# k-mer class-distribution statistics for a hawk.cpp comparison of two
# ASSEMBLIES (sample vs reference), whose canonical PAM 23-mers have been
# classified into four bucket files.
#
# NOTE: inputs are ASSEMBLY k-mer counts = genomic COPY NUMBERS (canonical,
# jellyfish -C, PAM-filtered ^CC || GG$). They are NOT raw reads, so there is
# no sequencing-error / coverage confound and count==1 means a genuine
# single-copy k-mer. (A reads-based, EM-normalized version would be a separate
# calculation.)
#
# Expects, inside DATA_DIR:
#   case_sorted_files.txt      -> points to CASE  = sample-assembly PAM dump(s)
#   control_sorted_files.txt   -> points to CONTROL = reference-assembly PAM dump(s)
#   present_v_absent.txt       sample>0, ref=0   (in sample, not in reference)
#   multiple_v_single.txt      sample>1, ref=1   (dup in sample, single in ref)
#   ref_present_v_absent.txt   sample=0, ref>0   (in reference, not in sample)
#   ref_multiple_v_single.txt  sample=1, ref>1   (single in sample, dup in ref)
# Each bucket file: SEQ <ws> caseCount <ws> controlCount
#
# Universe (denominator) = distinct canonical PAM 23-mers in sample OR ref =
# union of the two dumps named by the *_sorted_files.txt pointers. The four
# buckets are mutually disjoint; "none" = universe - sum(buckets).
#
# Usage:  krispmer_class_stats.sh [DATA_DIR] [DATASET_LABEL]
#   DATA_DIR       dir with the files above   (default: .)
#   DATASET_LABEL  name shown in the header    (default: basename of DATA_DIR)
#
# To save output, just redirect:  krispmer_class_stats.sh . n22 > july-9/n22_stats.txt
# ============================================================================
set -euo pipefail

DATA_DIR="${1:-.}"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
DATASET="${2:-$(basename "$DATA_DIR")}"

B1="$DATA_DIR/present_v_absent.txt"       # sample>0, ref=0
B2="$DATA_DIR/multiple_v_single.txt"      # sample>1, ref=1
B3="$DATA_DIR/ref_present_v_absent.txt"   # sample=0, ref>0
B4="$DATA_DIR/ref_multiple_v_single.txt"  # sample=1, ref>1

for f in "$DATA_DIR/case_sorted_files.txt" "$DATA_DIR/control_sorted_files.txt" \
         "$B1" "$B2" "$B3" "$B4"; do
  [ -f "$f" ] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

pct() { awk -v n="$1" -v d="$2" 'BEGIN{ if(d>0) printf "%.2f", 100*n/d; else printf "NA" }'; }

# ---- bucket sizes ----------------------------------------------------------
c1=$(wc -l < "$B1"); c2=$(wc -l < "$B2"); c3=$(wc -l < "$B3"); c4=$(wc -l < "$B4")
sum=$((c1 + c2 + c3 + c4))

# ---- universe = distinct canonical PAM 23-mers in sample U ref -------------
echo "[*] computing universe (sort -u over case + control dumps; slow step)..." >&2
tmpu="$(mktemp)"
(
  cd "$DATA_DIR"
  while read -r p; do [ -n "$p" ] && awk '{print $1}' "$p"; done < case_sorted_files.txt
  while read -r p; do [ -n "$p" ] && awk '{print $1}' "$p"; done < control_sorted_files.txt
) > "$tmpu"
union=$(LC_ALL=C sort -u "$tmpu" | wc -l)
rm -f "$tmpu"
none=$((union - sum))

# ---- print -----------------------------------------------------------------
printf '================================================================================\n'
printf 'kRISP-meR  -  k-mer Class Distribution   (dataset: %s)\n' "$DATASET"
printf 'assembly-vs-assembly, canonical PAM 23-mers   data dir: %s\n' "$DATA_DIR"
printf '================================================================================\n\n'

printf 'CLASS DISTRIBUTION   (universe = distinct PAM 23-mers in sample U ref)\n'
printf '%-22s %-15s %14s %8s\n' "Bucket" "(sample,ref)" "Count" "%"
printf '%-22s %-15s %14d %7s%%\n' "present_v_absent"       "sample>0,ref=0" "$c1"   "$(pct "$c1"   "$union")"
printf '%-22s %-15s %14d %7s%%\n' "multiple_v_single"      "sample>1,ref=1" "$c2"   "$(pct "$c2"   "$union")"
printf '%-22s %-15s %14d %7s%%\n' "ref_present_v_absent"   "sample=0,ref>0" "$c3"   "$(pct "$c3"   "$union")"
printf '%-22s %-15s %14d %7s%%\n' "ref_multiple_v_single"  "sample=1,ref>1" "$c4"   "$(pct "$c4"   "$union")"
printf '%-22s %-15s %14d %7s%%\n' "none of the four"       "concordant"     "$none" "$(pct "$none" "$union")"
printf '%-22s %-15s %14d %7s%%\n' "UNIVERSE"               ""               "$union" "100.00"
printf 'consistency: %d + %d = %d  [%s]\n' "$sum" "$none" "$((sum + none))" \
       "$( [ "$((sum + none))" -eq "$union" ] && echo OK || echo MISMATCH )"
