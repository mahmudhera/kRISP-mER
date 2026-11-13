#!/bin/bash
# k-mer counting and filtering for PAM (NGG/CCN)

# Input FASTA genome file
GENOME="Oryza_sativa_n22.OsN22RS2.dna.toplevel.fa"  ## input-filename

# Output files
JF_OUT="mer_counts.jf"
DUMP_OUT="mer_counts_dumps.txt"
FILTERED_OUT="mer_counts_PAM.txt"

echo "Starting counting ... "
# Step 1: Count 23-mers with Jellyfish
jellyfish count -m 23 -s 500M -t 16 -C $GENOME -o $JF_OUT -L 0 -U 999999
echo "Starting to dump .... "
# Step 2: Dump counts into text format
jellyfish dump -c $JF_OUT > $DUMP_OUT
echo "Starting to filter ... "
# Step 3: Filter for PAM sequences (NGG or CCN)
awk '{if($1 ~ /^CC/ || $1 ~ /GG$/) print $0}' $DUMP_OUT > $FILTERED_OUT

echo "Done! Filtered k-mers with PAM saved in $FILTERED_OUT"

