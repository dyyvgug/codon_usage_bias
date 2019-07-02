#!/bin/bash

for item in $(ls *.fq)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta --rna-strandness F -x /media/hp/disk1/DYY/reference/index/hg19/genome_index -U ./${item%.*}.fq -S ./${item%.*}.sam 2>>./ribo_mapping_repo.txt
	
done
