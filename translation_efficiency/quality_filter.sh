#!/bin/bash
mkdir clipper_quality
for item in $(ls *.sra)
do
	echo "quality_${item%.*}"
	fastq_quality_filter -q 20 -p 80 -Q 33 -i ./${item%.*}_tri.fastq -o ./clipper_quality/${item%.*}_qua.fastq
done

