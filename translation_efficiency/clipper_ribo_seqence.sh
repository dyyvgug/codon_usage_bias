#!/bin/bash
cd ./fq
mkdir ../clipper_fastq
for item in $(ls *.sra)
do
	echo "clipper_${item%.*}"
	fastx_clipper -a CTGTAGGCA -l 20 -d 0 -Q 33 -i ./${item%.*}.sra.fastq -o ../clipper_fastq/${item%.*}_tri.fastq
	


done

