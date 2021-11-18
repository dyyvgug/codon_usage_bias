#!/bin/bash
mkdir ./aligned
for item in $(ls *.sra)
do
        echo "bowtie2_${item%.*}"
        bowtie2 --sensitive -p 60 -N 1 -x ./sacCer3-tRNAs/index_one/tRNA_index -U ./clipper_fastq/${item%.*}_trimmed.fq -S ./aligned/${item%.*}.sam 2>> ./aligned/mapping_repo.txt

done
