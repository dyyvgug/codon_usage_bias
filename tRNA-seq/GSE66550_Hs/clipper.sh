#!/bin/bash

mkdir clipper_fastq

for item in $(ls *.sra)
do
	echo "clipper_${item%.*}"
	trim_galore -q 30 --phred33 --length 15 -e 0.1 --stringency 3 ./fastq/${item%.*}.sra.fastq -o ./clipper_fastq/${item%.*}.fastq
done
