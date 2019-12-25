#!/bin/bash
for item in $(ls *.fastq)
do
        echo "bow_${item%.*}"
        bowtie2 -p 30 -N 2 -x /media/hp/disk2/DYY2/tRNA/Sc/sacCer3-tRNAs/index_one/tRNA_index -U ./${item%.*}.fastq -S ../aligned_tRNA_one/${item%.*}.sam 2>> ../aligned_tRNA_one/mapping_repo.txt
        
done
