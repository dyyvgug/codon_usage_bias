#!/bin/bash
cp Ce_anticodon_count.txt Ce_copy_num.txt
sed -i 's/Ala.* A/A/g' Ce_copy_num.txt
sed -i 's/Gly.* A/A/g' Ce_copy_num.txt
sed -i 's/Pro.* A/A/g' Ce_copy_num.txt
sed -i 's/Thr.* A/A/g' Ce_copy_num.txt
sed -i 's/Val.* A/A/g' Ce_copy_num.txt
sed -i 's/Ser.* A/A/g' Ce_copy_num.txt
sed -i 's/Arg.* A/A/g' Ce_copy_num.txt
sed -i 's/Leu.* A/A/g' Ce_copy_num.txt
sed -i 's/Phe.* A/A/g' Ce_copy_num.txt
sed -i 's/Asn.* A/A/g' Ce_copy_num.txt
sed -i 's/Lys.* C/C/g' Ce_copy_num.txt
sed -i 's/Glu.* C/C/g' Ce_copy_num.txt
sed -i 's/Gln.* C/C/g' Ce_copy_num.txt
sed -i 's/Met.* C/C/g' Ce_copy_num.txt
sed -i 's/Trp.* C/C/g' Ce_copy_num.txt
sed -i 's/Asp.* A/A/g' Ce_copy_num.txt
sed -i 's/His.* A/A/g' Ce_copy_num.txt
sed -i 's/Ile.* A/A/g' Ce_copy_num.txt
sed -i 's/Tyr.* A/A/g' Ce_copy_num.txt
sed -i 's/Supres.* C/C/g' Ce_copy_num.txt
sed -i 's/SelCys.* T/T/g' Ce_copy_num.txt
sed -i '/Isotype/d' Ce_copy_num.txt
sed -i 's/         /\t/g' Ce_copy_num.txt
sed -i 's/      /\n/g' Ce_copy_num.txt
sed -i 's/\t/\n/g' Ce_copy_num.txt
sed -i 's/ //g' Ce_copy_num.txt
sed -i '/Supres/d' Ce_copy_num.txt
sed -i '/:/!d' Ce_copy_num.txt


