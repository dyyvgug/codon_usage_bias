sed 's/title/transcription_id/g' CBI_CAI_bycodonW.txt > temp
mv temp CBI_CAI_bycodonW.txt
sed 's/           	/\t/g' CBI_CAI_bycodonW.txt > temp
mv temp CBI_CAI_bycodonW.txt
sed 's/\t\t/\t/g' CBI_CAI_bycodonW.txt > temp
mv temp CBI_CAI_bycodonW.txt
sed -i 's/         //g' CBI_CAI_bycodonW.txt
sed -i 's/      //g' CBI_CAI_bycodonW.txt
sed -i 's/\t\t/\t/g' CBI_CAI_bycodonW.txt

