#==========================================================================================
# Yingying Dong.2019-12-18.tRNA seq anticodon & the reference set codon
#==========================================================================================
library(dplyr)
require(Biostrings)
library(ggplot2)
spA = "Sc"
setwd(paste0("G:\\学习专用\\tRNA\\Sc\\seqbackup\\aligned_tRNA_one"))
bed = read.table("ERR2382482.bed", sep = "\t")
b1 = as.data.frame(table(bed$V1))
b1$id = gsub(".+tRNA-", "", b1$Var1)
b1$id = gsub("-\\d-\\d+$", "", b1$id)
b1$codon = sub(".+-", "", b1$id)
b2 = b1 %>% group_by(codon) %>% summarise(Frq = sum(Freq))
b2$Codon = as.character(reverseComplement( DNAStringSet(b2$codon) ))
names(b2) = c("anticodon","anti_fre","codon")
View(b2)
rscu = read.table("ribo_codon_fre_RSCU.txt",header = T, stringsAsFactors = F)
bed_rscu = merge(b2,rscu,by = "codon", all = T)
View(bed_rscu)
bed_rscu = bed_rscu[-grep("NNN", bed_rscu$anticodon),]
bed_rscu = bed_rscu[-grep("STOP",bed_rscu$AA),]
bed_rscu = bed_rscu[-grep("M",bed_rscu$AA),]
bed_rscu = bed_rscu[-grep("W",bed_rscu$AA),]
View(bed_rscu)
#bed_rscu[is.na(bed_rscu)] = 0
bed_rscu = bed_rscu[complete.cases(bed_rscu),]
bed_rscu$anticodon = as.character(reverseComplement( DNAStringSet(bed_rscu$codon) ))
write.csv(bed_rscu,"tRNAseq_anti_ref_codon_full.csv",quote = F,row.names = F)
anti_rscu_p = cor.test(log2(bed_rscu$anti_fre+1),bed_rscu$hits)
anti_rscu_p
p1 <- ggplot(bed_rscu,aes(x = log2(bed_rscu$anti_fre),y = bed_rscu$hits,color = bed_rscu$AA))+
  geom_point(shape = 16,size = 2.5)+
  labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon count)")+
  ylab("Reference set codon count")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p1
library(MASS)
p1 <- ggplot(bed_rscu,aes(x = log2(bed_rscu$anti_fre),y = bed_rscu$hits,color = bed_rscu$AA))+
  geom_point(shape = 16,size = 2.5)+
  labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon count)")+
  ylab("Reference set codon count")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
#==========================================================================================
# Yingying Dong.2019-12-18.tRNA seq anticodon & the reference set codon
#==========================================================================================
library(dplyr)
require(Biostrings)
library(ggplot2)
library(MASS)
?trans_breaks
library(scales)
p1 <- ggplot(bed_rscu,aes(x = log2(bed_rscu$anti_fre),y = bed_rscu$hits,color = bed_rscu$AA))+
  geom_point(shape = 16,size = 2.5)+
  labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon count)")+
  ylab("Reference set codon count")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p1
p1 <- ggplot(bed_rscu,aes(x = log2(bed_rscu$anti_fre),y = bed_rscu$hits,color = bed_rscu$AA))+
  geom_point(shape = 16,size = 2.5)+
  labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon counts)")+
  ylab("Cytosolic RP codon counts")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p1
ggsave("Sc_RP_anticodon_codon.pdf")
write.csv(bed_rscu,"tRNAseq_RP_anti_ref_codon.csv",quote = F,row.names = F)
whole = read.table("codon_frequency.txt")
names(whole) = c("codon","aa","whole_fre")
bed_whole = merge(b2,whole,by = "codon")
bed_whole = bed_whole[-grep("M",bed_whole$aa),]
bed_whole = bed_whole[-grep("W",bed_whole$aa),]
anti_whole_p = cor.test(log2(bed_whole$anti_fre),bed_whole$whole_fre)
anti_whole_p
p2 <- ggplot(bed_whole,aes(x = log2(bed_whole$anti_fre),y = bed_whole$whole_fre,color = bed_whole$aa))+
  geom_point(shape = 16,size = 2.5)+
  labs(title = paste0(spA," cor_anticodon_genome_codon    ","r=",round(anti_whole_p$estimate,5),"  p=",round(anti_whole_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon countS)")+
  ylab("Whole genome codon countS")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))
p2
ggsave("Sc_genome_codon.pdf")
write.csv(bed_rscu,"tRNAseq_anti_genome_codon.txt",quote = F,row.names = F)
