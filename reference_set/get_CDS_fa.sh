
sed 's/gene.*//g' transcripts.fa > transcripts_only.fa
sed 's/[a-z]/\u&/g' transcripts_only.fa > CDS.fa
