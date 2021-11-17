#!/bin/bash
for item in $(ls *.fastq)
do
        echo "bowtie2_${item%.*}"
        bowtie2 --sensitive -p 60 -N 1 -x ../sacCer3-tRNAs/index_one/tRNA_index -U ./${item%.*}.fastq -S ../aligned_tRNA_one/${item%.*}.sam 2>> ../aligned_tRNA_one/mapping_repo.txt
        
done
