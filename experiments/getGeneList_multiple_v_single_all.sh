prefix_1=multiple_v_single
ref=exons

grep 'XM:i:0' ${prefix_1}_${ref}_all.sam > ${prefix_1}_${ref}_all_max_0.sam
wc -l ${prefix_1}_${ref}_all_max_0.sam
fileName=${prefix_1}_${ref}_all_max_0.sam
outputFile=geneIDs_${prefix_1}_${ref}_all_max_0.txt
awk '{print $3}' $fileName | sed 's/^rna-//' > $outputFile

wc -l $outputFile

sort geneIDs_${prefix_1}_${ref}_all_max_0.txt | uniq > geneIDs_${prefix_1}_${ref}_all_max_0_sorted_uniq.txt

wc -l geneIDs_${prefix_1}_${ref}_all_max_0_sorted_uniq.txt


prefix_2=unique_PAM_gRNAs_n22

grep 'XM:i:0' ${prefix_2}_${ref}_all.sam > ${prefix_2}_${ref}_all_max_0.sam
fileName=${prefix_2}_${ref}_all_max_0.sam
outputFile=${prefix_2}_${ref}_all_max_0_genes.txt
awk '{print $3}' $fileName | sed 's/^rna-//' > $outputFile

