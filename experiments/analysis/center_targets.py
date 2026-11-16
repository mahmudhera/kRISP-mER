#!/usr/bin/env python3
# center_targets.py CLASS
# Slice each target in targets/<CLASS>/ down to a 150bp window centered on the
# class guide's binding position (taken from the runStageB alignment SAM),
# overwriting the target file in place. Matches the manuscript target unit
# (150bp centered on the guide) so runtime + results line up with feb-20.
import os, sys
from Bio import SeqIO

CLASS = sys.argv[1]
WIN   = 150
HALF  = WIN // 2

# which runStageB SAM holds this class's guide->transcript alignments
SAM = {
    "multiple_v_single":     "multiple_v_single_exons_max_0.sam",
    "present_v_absent":      "present_v_absent_exons_max_3.sam",   # CHM13-only -> <=3mm
    "ref_multiple_v_single": "ref_multiple_v_single_exons_max_0.sam",
    "ref_present_v_absent":  "ref_present_v_absent_exons_max_0.sam",
}[CLASS]

TGTDIR = f"targets/{CLASS}"

# 1) transcript ids we actually need (the 50 selected targets)
need = {}
for fn in os.listdir(TGTDIR):
    if fn.endswith(".fasta"):
        rec = next(SeqIO.parse(os.path.join(TGTDIR, fn), "fasta"))
        need[rec.id] = fn                      # transcript-id -> filename

# 2) one pass over the SAM; first alignment POS per needed transcript
pos = {}
with open(SAM) as f:
    for line in f:
        if line[0] == "@":
            continue
        c = line.split("\t")
        if len(c) < 4:
            continue
        rname = c[2]
        if rname in need and rname not in pos:
            pos[rname] = int(c[3])             # 1-based leftmost pos of the guide
            if len(pos) == len(need):
                break

# 3) slice + overwrite
done = skip = 0
for tx, fn in need.items():
    path = os.path.join(TGTDIR, fn)
    rec  = next(SeqIO.parse(path, "fasta"))
    seq  = str(rec.seq).upper()
    if tx not in pos:
        skip += 1
        continue
    center = (pos[tx] - 1) + 11                # ~middle of the 23bp guide
    start  = max(0, center - HALF)
    end    = min(len(seq), start + WIN)
    start  = max(0, end - WIN)                 # keep full 150 near the 3' end
    with open(path, "w") as out:
        out.write(f">{rec.id}\n{seq[start:end]}\n")
    done += 1

print(f"{CLASS}: centered {done}, skipped {skip} (no SAM hit)")
