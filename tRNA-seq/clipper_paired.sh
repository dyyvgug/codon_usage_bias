#!/bin/bash
mkdir clipper_fastq
cd fastq
gunzip *.gz

ls *_1.fastq >1
ls *_2.fastq >2

paste 1 2 > config
cat config |while read id
do
	echo "clipper_${id}"
	arr=(${id})
    fq1=${arr[0]}
    fq2=${arr[1]}
	trim_galore -q 33 --phred33 --length 15 -e 0.1 --stringency 3 --paired $fq1 $fq2 -o ../clipper_fastq/
done

rename 's/\.sr.*_1/_1/g' *
rename 's/\.sr.*_2/_2/g' *
