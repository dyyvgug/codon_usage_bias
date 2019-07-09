#!/bin/bash
mkdir clipper_fastq2
for item in $(ls *.sra)
do
	echo "clipper_${item%.*}"
	fastx_clipper -a TGGAATT -l 20 -d 0 -Q 33 -i ./fastq/${item%.*}.sra.fastq -o ./clipper_fastq2/${item%.*}_tri.fastq
done
