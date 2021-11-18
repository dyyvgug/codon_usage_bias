#!/bin/bash


for item in $(ls *.sam)
do
	echo "sam_${item%.*}"		
	
	samtools sort -@ 8 -o ${item%.*}.bam ${item%.*}.sam
	
	samtools flagstat ${item%.*}.bam > ${item%.*}.bam.flag
	
	bamToBed -i ${item%.*}.bam > ${item%.*}.bed
done
