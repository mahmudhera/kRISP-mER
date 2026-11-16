#!/usr/bin/env python3
# circleseq_parse.py INPUT [OUTPUT.tsv]
# Parse a CIRCLE-seq off-target table (Nature supp MOESM5 'circleSeq' sheet, or
# the pipeline "identified" tsv) into a clean per-site table for Path-A mining.
# Emits: targetsite cell guide_seq offtarget_seq n_mm circle_reads chrom start end is_ontarget
import sys
import pandas as pd

IN  = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "circleseq_sites.tsv"

if IN.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(IN)               # first sheet (MOESM5 = 'circleSeq')
else:
    df = pd.read_csv(IN, sep=None, engine="python")

# accept either the supp-table names or the pipeline names
ALIASES = {
    "offtarget_seq": ["Off-target Sequence", "Site_SubstitutionsOnly.Sequence"],
    "circle_reads":  ["Read", "ReadCount"],
    "n_mm":          ["Distance", "Site_SubstitutionsOnly.NumSubstitutions"],
    "guide_seq":     ["TargetSequence"],
    "targetsite":    ["Targetsite"],
    "cell":          ["Cell"],
    "chrom":         ["Chromosome"],
    "start":         ["Start"],
    "end":           ["End"],
}

def pick(names):
    for n in names:
        if n in df.columns:
            return df[n]
    return None

cols = {k: pick(v) for k, v in ALIASES.items()}
for req in ("offtarget_seq", "circle_reads"):
    if cols[req] is None:
        sys.exit(f"ERROR: no column for {req}. present: {list(df.columns)}")

out = pd.DataFrame({k: (v if v is not None else "") for k, v in cols.items()})
out["offtarget_seq"] = out["offtarget_seq"].astype(str).str.upper().str.replace("-", "", regex=False)

# keep real 23-ish-mers; drop blanks / gap-only artifacts
out = out[out["offtarget_seq"].str.len() >= 20].copy()
out["is_ontarget"] = (pd.to_numeric(out["n_mm"], errors="coerce") == 0)

out.to_csv(OUT, sep="\t", index=False)

print(f"parsed {len(out)} rows "
      f"({(~out['is_ontarget']).sum()} off-target, {out['is_ontarget'].sum()} on-target)")
if cols["cell"] is not None:
    print("\ncell breakdown:\n", out["cell"].value_counts())
print("\nper guide:\n", out["targetsite"].value_counts())
print(f"\nwrote {OUT}")
