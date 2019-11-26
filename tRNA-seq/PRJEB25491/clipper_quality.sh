#!/bin/bash
cd ./fastq
mkdir ../clipper_fastq
mkdir ../quality
for item in $(ls *.sra)
do
	echo "clipper_${item%.*}"
	fastx_clipper -a TTCTGCTTGAAAAAAA -l 20 -d 0 -Q 33 -i ./${item%.*}.sra.fastq -o ../clipper_fastq/${item%.*}_tri.fastq

	echo "quality_${item%.*}"
	fastq_quality_filter -q 20 -p 80 -Q 33 -i ../clipper_fastq/${item%.*}_tri.fastq -o ../quality/${item%.*}_qua.fastq
done

cd ../quality
rename 's/_qua//' *
