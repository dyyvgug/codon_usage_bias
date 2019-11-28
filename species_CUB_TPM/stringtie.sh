#!/bin/bash

conda activate bioinfo
for item in $(ls *.sam)
do
	echo "sam_${item%.*}"		
	
	samtools sort -@ 8 -o ${item%.*}.bam ${item%.*}.sam
	
	echo "str_${item%.*}"
	stringtie -p 30 --rf -e -G /media/hp/disk1/DYY/reference/annotation/$1/ref.gtf -o ./${item%.*}.gtf -l ${item%.*} ${item%.*}.bam -A ${item%.*}_abund.out

done


