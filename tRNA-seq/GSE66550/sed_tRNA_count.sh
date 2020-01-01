#!/bin/bash
for i in $(ls *.txt)
do
	echo "$i"	
	sed -i 's/  /\t/g' $i
	sed -i 's/.*Val//g' $i
	sed -i 's/.*Asp//g' $i
	sed -i 's/.*Met//g' $i
	sed -i 's/.*Gly//g' $i
	sed -i 's/.*Gln//g' $i
	sed -i 's/.*Phe//g' $i
	sed -i 's/.*Lys//g' $i
	sed -i 's/.*His//g' $i
	sed -i 's/.*Glu//g' $i
	sed -i 's/.*Leu//g' $i
	sed -i 's/.*Asn//g' $i
	sed -i 's/.*Arg//g' $i
	sed -i 's/.*Ser//g' $i
	sed -i 's/.*Ala//g' $i
	sed -i 's/.*Ile//g' $i
	sed -i 's/.*Thr//g' $i
	sed -i 's/.*Trp//g' $i
	sed -i 's/.*Pro//g' $i
	sed -i 's/.*Tyr//g' $i
	sed -i 's/.*Cys//g' $i
	sed -i 's/.*SeC//g' $i
	sed -i 's/.*Sup//g' $i
	sed -i '/Undet/d' $i
	sed -i 's/total: //g' $i
	sed -i '2d' $i



done
