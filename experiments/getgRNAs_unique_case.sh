awk '($1 ~ /^CC/ || $1 ~ /GG$/) && $2 == 1 {print $0}' lexo_mer_counts_PAM_sorted_n22.txt > unique_PAM_gRNAs_n22.txt
awk '{printf(">%d\n%s\n", NR, $1)}' unique_PAM_gRNAs_n22.txt > unique_PAM_gRNAs_n22.fasta
