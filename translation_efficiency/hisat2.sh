#!/bin/bash

for item in $(ls *.sra)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta --rna-strandness RF -x /media/hp/disk1/DYY/reference/index/hg19/genome_index -1 ~/Desktop/SRP136094/fq/${item%.*}_1.fastq.gz -2 ~/Desktop/SRP136094/fq/${item%.*}_2.fastq.gz -S ~/Desktop/SRP136094/${item%.*}.sam
done
