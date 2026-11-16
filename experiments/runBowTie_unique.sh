file=unique_PAM_gRNAs_n22
ref=exons
# /home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2 -p 16 -a -x ${ref}  -f ${file}.fasta --no-unal -S ${file}_${ref}_all.sam
/home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2 -p 16 -a -x ${ref} -f ${file}.fasta -S ${file}_${ref}_all.sam --no-unal