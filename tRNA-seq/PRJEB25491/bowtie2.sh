#!/bin/bash
mkdir ../aligned_tRNA
for item in $(ls *.fastq)
do
        echo "bow_${item%.*}"
        bowtie2 -p 30 -N 1 -x /media/hp/disk2/DYY2/tRNA/Sc/sacCer3-tRNAs/index/tRNA_index -U ./${item%.*}.fastq -S ../aligned_tRNA/${item%.*}.sam 2>> ../aligned_tRNA/mapping_repo.txt
        
done



