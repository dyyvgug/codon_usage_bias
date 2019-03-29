bedtools closest -a GSM1480325_K562_GROseq_minus.bedGraph -b ref.bed > GSM1480325_minus_closer.txt
bedtools closest -a GSM1480325_K562_GROseq_plus.bedGraph -b ref.bed > GSM1480325_plus_closer.txt
#注释文件来自于gtf，gtf to bed 的命令gtf2bed < foo.gtf > foo.bed命令不可用，bash提示格式出错，所以使用awak， awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' input.gtf | gtf2bed - > output.bed
bedtools closest -a 0325out.dREG.peak.score.bed -b ref_sort.bed > 0325peak_closer_sort.txt
#注释文件来自于gff，gff to bed 的命令为 sortBed -i myfile.gff | gff2bed > my_sorted_file.bed.
bedtools closest -a 0325out.dREG.peak.score.bed -b ref.bed > 0325peak_closer.txt	
