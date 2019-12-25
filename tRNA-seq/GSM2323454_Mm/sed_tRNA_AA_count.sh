#!/bin/bash
for i in $(ls *.txt)
do
	echo "$i"	
	sed -i 's/Mus.*-//g' $i
	sed -i 's/  /\t/g' $i
	sed -i 's/Val.*\t/.*\t/g' $i
	sed -i 's/Asp.*\t/.*\t/g' $i
	sed -i 's/Met.*\t/.*\t/g' $i
	sed -i 's/Gly.*\t/.*\t/g' $i
	sed -i 's/Gln.*\t/.*\t/g' $i
	sed -i 's/Phe.*\t/.*\t/g' $i
	sed -i 's/Lys.*\t/.*\t/g' $i
	sed -i 's/His.*\t/.*\t/g' $i
	sed -i 's/Glu.*\t/.*\t/g' $i
	sed -i 's/Leu.*\t/.*\t/g' $i
	sed -i 's/Asn.*\t/.*\t/g' $i
	sed -i 's/Arg.*\t/.*\t/g' $i
	sed -i 's/Ser.*\t/.*\t/g' $i
	sed -i 's/Ala.*\t/.*\t/g' $i
	sed -i 's/Ile.*\t/.*\t/g' $i
	sed -i 's/Thr.*\t/.*\t/g' $i
	sed -i 's/Trp.*\t/.*\t/g' $i
	sed -i 's/Pro.*\t/.*\t/g' $i
	sed -i 's/Tyr.*\t/.*\t/g' $i
	sed -i 's/Cys.*\t/.*\t/g' $i
	sed -i 's/SeC.*\t/.*\t/g' $i
	sed -i 's/Sup.*\t/.*\t/g' $i
	sed -i '/Undet/d' $i
	sed -i 's/total: //g' $i


done
