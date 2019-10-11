#!/bin/bash
mkdir ../aligned
for item in $(ls *.sra)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta -x /media/hp/disk1/DYY/reference/index/$1/genome_index -1 ./${item%.*}.sra_1.fastq -2 ./${item%.*}.sra_2.fastq -S ../aligned/${item%.*}.sam 2>> ../aligned/mapping_repo.txt
done
