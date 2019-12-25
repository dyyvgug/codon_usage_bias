#!/bin/bash
cp copy_num.txt copy_num_backup.txt
sed -i 's/Ala.* A/A/g' copy_num.txt
sed -i 's/Gly.* A/A/g' copy_num.txt
sed -i 's/Pro.* A/A/g' copy_num.txt
sed -i 's/Thr.* A/A/g' copy_num.txt
sed -i 's/Val.* A/A/g' copy_num.txt
sed -i 's/Ser.* A/A/g' copy_num.txt
sed -i 's/Arg.* A/A/g' copy_num.txt
sed -i 's/Leu.* A/A/g' copy_num.txt
sed -i 's/Phe.* A/A/g' copy_num.txt
sed -i 's/Asn.* A/A/g' copy_num.txt
sed -i 's/Lys.* C/C/g' copy_num.txt
sed -i 's/Asp.* A/A/g' copy_num.txt
sed -i 's/Glu.* C/C/g' copy_num.txt
sed -i 's/His.* A/A/g' copy_num.txt
sed -i 's/Gln.* C/C/g' copy_num.txt
sed -i 's/Ile.* A/A/g' copy_num.txt
sed -i 's/Met.* C/C/g' copy_num.txt
sed -i 's/Tyr.* A/A/g' copy_num.txt
sed -i 's/Supres.* C/C/g' copy_num.txt
sed -i 's/Cys.* A/A/g' copy_num.txt
sed -i 's/Trp.* C/C/g' copy_num.txt
sed -i 's/SelCys.* T/T/g' copy_num.txt
sed -i '/Isotype/d' copy_num.txt
sed -i 's/         /\t/g' copy_num.txt
sed -i 's/      /\n/g' copy_num.txt
sed -i 's/\t/\n/g' copy_num.txt
sed -i 's/ //g' copy_num.txt
sed -i '/:/!d' copy_num.txt
sed -i 's/:/\t/g' copy_num.txt

Rscript comp_copy_num.R
