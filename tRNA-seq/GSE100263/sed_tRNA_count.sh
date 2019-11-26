#!/bin/bash
for i in $(ls *.csv)
do
	echo "$i"	
	sed -i 's/tRNA-//g' $i
	sed -i 's/-.*,/,/g' $i
	sed -i '1d' $i
	
done
