#!/bin/bash

cp exportstAiGenes.csv exportstAiGenesBackup.csv
sed -i 's/ATG.*//g' exportstAiGenesBackup.csv
sed -i /^$/d exportstAiGenesBackup2.csv
R
sed -i '/,$/N; s/\n//' tAI.txt
sed -i 's/,,\"/\t/g' tAI.txt
sed -i 's/"//g' tAI.txt

