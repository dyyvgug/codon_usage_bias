#!/bin/bash

for item in $(ls *.sra)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta --rna-strandness F -x /media/hp/disk1/DYY/reference/index/C_elegans_Ensl_WBcel235/genome_index -U ./fastq/${item%.*}.sra.fastq -S ./${item%.*}.sam 2>> ./ribo_mapping_repo.txt
	
done
