#!/bin/bash
for i in $(ls *.txt)
do
	echo "$i"	
	sed -i 's/Homo.*-//g' $i
	sed -i 's/Val.*\t/Val\t/g' $i
	sed -i 's/Asp.*\t/Asp\t/g' $i
	sed -i 's/Met.*\t/Met\t/g' $i
	sed -i 's/Gly.*\t/Gly\t/g' $i
	sed -i 's/Gln.*\t/Gln\t/g' $i
	sed -i 's/Phe.*\t/Phe\t/g' $i
	sed -i 's/Lys.*\t/Lys\t/g' $i
	sed -i 's/His.*\t/His\t/g' $i
	sed -i 's/Glu.*\t/Glu\t/g' $i
	sed -i 's/Leu.*\t/Leu\t/g' $i
	sed -i 's/Asn.*\t/Asn\t/g' $i
	sed -i 's/Arg.*\t/Arg\t/g' $i
	sed -i 's/Ser.*\t/Ser\t/g' $i
	sed -i 's/Ala.*\t/Ala\t/g' $i
	sed -i 's/Ile.*\t/Ile\t/g' $i
	sed -i 's/Thr.*\t/Thr\t/g' $i
	sed -i 's/Trp.*\t/Trp\t/g' $i
	sed -i 's/Pro.*\t/Pro\t/g' $i
	sed -i 's/Tyr.*\t/Tyr\t/g' $i
	sed -i 's/Cys.*\t/Cys\t/g' $i
	sed -i 's/SeC.*\t/STOP\t/g' $i
	sed -i 's/Sup.*\t/STOP\t/g' $i
	sed -i '/Undet/d' $i
	sed -i '/-/d' $i
	sed -i '/Ecoli/d' $i
	sed -i '/Yeast/d' $i
	sed -i '1d' $i

done
