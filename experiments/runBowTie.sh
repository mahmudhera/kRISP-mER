file=multiple_v_single
ref=exons
/home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2 -p 16 -a -x ${ref}  -f ${file}.fasta --no-unal -S ${file}_${ref}_all.sam

file=present_v_absent
ref=exons
/home/atif/bowtie2-2.4.1-linux-x86_64/bowtie2 -p 16  -x ${ref}  -f ${file}.fasta --no-unal -S ${file}_${ref}_all.sam