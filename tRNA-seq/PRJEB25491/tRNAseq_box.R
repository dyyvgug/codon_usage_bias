setwd("/media/hp/disk2/DYY2/tRNA/Sc/seqbackup/aligned_tRNA_one/")
bed = read.table(file.choose(), sep = "\t")
b1 = as.data.frame(table(bed$V1))
b1$id = gsub(".+tRNA-", "", b1$Var1)
b1$id = gsub("-\\d-\\d+$", "", b1$id)
b1$codon = sub(".+-", "", b1$id)

library(dplyr)
require(Biostrings)
b2 = b1 %>% group_by(codon) %>% summarise(Frq = sum(Freq))
b2$Codon = as.character(reverseComplement( DNAStringSet(b2$codon) ))
names(b2) = c("anticodon","anti_fre","codon")

rscu = read.table("/media/hp/disk1/DYY/reference/annotation/Saccharomyces_cerevisiae/ribo_codon_fre_RSCU.txt",
                  header = T, stringsAsFactors = F)
bed_rscu = merge(b2,rscu,by = "codon", all = T)
bed_rscu = bed_rscu[-grep("N", bed_rscu$anticodon),]
bed_rscu[is.na(bed_rscu)] = 0
bed_rscu$anticodon = as.character(reverseComplement( DNAStringSet(bed_rscu$codon) ))

anti_rscu_r = cor(bed_rscu$anti_fre,bed_rscu$hits)
anti_rscu_p = cor.test(log2(bed_rscu$anti_fre),bed_rscu$hits)
anti_rscu_p
plot(log2(bed_rscu$anti_fre),bed_rscu$hits)

whole = read.table("/media/hp/disk1/DYY/reference/annotation/Saccharomyces_cerevisiae/ref/codon_frequency.txt")
names(whole) = c("codon","aa","whole_fre")
bed_whole = merge(b2,whole,by = "codon")
anti_whole_p = cor.test(log2(bed_whole$anti_fre),bed_whole$whole_fre)
anti_whole_p

#==================================================================================================
# Eight major degenerate codon pairs
#==================================================================================================
bed_rscu$deg = bed_rscu$anticodon
bed_rscu$deg = sub("GGC","AGC",bed_rscu$deg)
bed_rscu$deg = sub("GGG","AGG",bed_rscu$deg)
bed_rscu$deg = sub("GGT","AGT",bed_rscu$deg)
bed_rscu$deg = sub("GAC","AAC",bed_rscu$deg)
bed_rscu$deg = sub("GGA","AGA",bed_rscu$deg)
bed_rscu$deg = sub("GCG","ACG",bed_rscu$deg)
bed_rscu$deg = sub("GAG","AAG",bed_rscu$deg)
bed_rscu$deg = sub("GAT","AAT",bed_rscu$deg)
df <- data.frame(bed_rscu$AA,bed_rscu$hits,bed_rscu$anti_fre,bed_rscu$deg)
names(df) = c("AA","ribo_hits","anti_fre","deg")
df2 <- df %>% group_by(deg) %>% summarise_at(vars(ribo_hits,anti_fre), list(name = sum)) 
df2$codon = as.character(reverseComplement(DNAStringSet(df2$deg)))
df3 <- data.frame(bed_rscu$codon,bed_rscu$AA)
names(df3) = c("codon","AA")
df4 = merge(df2,df3,by = "codon",all.x = TRUE)
names(df4) = c("codon","deg","ribo_hits","anti_fre","AA")
