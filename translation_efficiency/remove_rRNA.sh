#!/bin/bash
for item in $(ls *.fastq)
do
	echo "{item%.*} remove rRNA"
	bowtie -v 3 --norc /media/hp/disk1/DYY/reference/rRNA/hg38_rRNA -q ${item%.*}.fastq --un ${item%.*}_rmrna.fastq
done

#paired-end 
# --un-conc ${i}_input_rmrRNA.fastq -1 ${i}_in.read1_Clean.fastq -2 ${i}_in.read2_Clean.fastq
