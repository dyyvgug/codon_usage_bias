#!/bin/bash
mkdir aligned
for item in $(ls *.sra)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta --rna-strandness F -x /media/hp/disk1/DYY/reference/index/hg19/genome_index -U ./${item%.*}.sra.fastq -S ./aligned/${item%.*}.sam 2>> ./aligned/ribo_mapping_repo.txt
	
done
