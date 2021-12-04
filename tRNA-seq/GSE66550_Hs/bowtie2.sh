#!/bin/bash
mkdir ./aligned
for item in $(ls *.sra)
do
        echo "bowtie2_${item%.*}"
        bowtie2 --sensitive -p 30 -N 1 -x ./hg38-tRNA/index/tRNA -U ./clipper_fastq/${item%.*}.fq -S ./aligned/${item%.*}.sam 2>> ./aligned/mapping_repo.txt

done
