#!/bin/bash
mkdir ./aligned
for item in $(ls *.sra)
do
        echo "bowtie2_${item%.*}"
        bowtie2 --sensitive -p 60 -N 1 -x ./mm39-tRNA/index/tRNA -1 ./clipper_fastq/${item%.*}_1.fq -2 ./clipper_fastq/${item%.*}_2.fq -S ./aligned/${item%.*}.sam 2>> ./aligned/mapping_repo.txt

done

