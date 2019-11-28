#!/bin/bash

cd /media/hp/disk2/DYY2/sra_fq_RNAseq/$1/experiment$2
for item in $(ls *.sra)
do
	echo "hi_${item%.*}"
	hisat2 -p 30 -t --dta --rna-strandness RF --max-intronlen 100000 -x /media/hp/disk1/DYY/reference/index/$1/genome_index -1 /media/hp/disk2/DYY2/sra_fq_RNAseq/$1/experiment$2/fastq/${item%.*}.sra_1.fastq.gz -2 /media/hp/disk2/DYY2/sra_fq_RNAseq/$1/experiment$2/fastq/${item%.*}.sra_2.fastq.gz -S /media/hp/Katniss/DYY/aligned/$1/experiment$2/${item%.*}.sam 2>>/media/hp/Katniss/DYY/aligned/$1/experiment$2/mapping_repo.txt
	
done
