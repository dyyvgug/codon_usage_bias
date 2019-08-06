#!/bin/bash
mkdir ../quality
for item in $(ls *.sra)
do
	echo "quality_${item%.*}"
	fastq_quality_filter -q 20 -p 80 -Q 33 -i ../clipper_fastq/${item%.*}_tri.fastq -o ../quality/${item%.*}_qua.fastq
done

