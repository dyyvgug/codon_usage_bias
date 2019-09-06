#!/bin/bash
mkdir ../aligned_ri
for item in $(ls *.fastq)
do
        echo "hi_${item%.*}"
        hisat2 -p 30 -t --dta --rna-strandness F -x /media/hp/disk1/DYY/reference/index/Solanum_lycopersicum/genome_index -U ./${item%.*}.sra.fastq -S ../aligned_ri/${item%.*}.sam 2>> ../aligned_ri/ribo_mapping_repo.txt
        
done
