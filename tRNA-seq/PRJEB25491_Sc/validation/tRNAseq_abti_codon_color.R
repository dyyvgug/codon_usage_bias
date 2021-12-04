#==========================================================================================
# Yingying Dong.2019-12-18.modify:2021.tRNA seq anticodon & the reference set codon
#==========================================================================================
library(dplyr)
if (!require("Biostrings"))
  BiocManager::install("Biostrings")
require(Biostrings)
library(ggplot2)
library(scales)
spA = "Sc"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(paste0("G:\\Ñ§Ï°×¨ÓÃ\\tRNA\\Sc\\seqbackup\\aligned_tRNA_one"))
bed = read.table("ERR2382482.bed", sep = "\t")
b1 = as.data.frame(table(bed$V1))
b1$id = gsub(".+tRNA-", "", b1$Var1)
b1$id = gsub("-\\d-\\d+$", "", b1$id)
b1$codon = sub(".+-", "", b1$id)

b2 = b1 %>% group_by(codon) %>% summarise(Frq = sum(Freq))
b2$Codon = as.character(reverseComplement( DNAStringSet(b2$codon) ))
names(b2) = c("anticodon","anti_fre","codon")

rscu_array<-c("Cytosolic RP","Mitochondria RP","HTT")
for (i in rscu_array) {
  if(F){
    i = "Mitochondria RP"
    # color = rgb(200,175,212,maxColorValue = 255)
    i = "HTT"
    # rgb(252,215,24)
  }
  rscu = read.table(paste0(i,"_codon_fre_RSCU.txt"),header = T, stringsAsFactors = F)
  bed_rscu = merge(b2,rscu,by = "codon", all = T)
  bed_rscu = bed_rscu[-grep("NNN", bed_rscu$anticodon),]
  bed_rscu = bed_rscu[-grep("STOP",bed_rscu$AA),]
  bed_rscu = bed_rscu[-grep("M",bed_rscu$AA),]
  bed_rscu = bed_rscu[-grep("W",bed_rscu$AA),]
  
  #bed_rscu[is.na(bed_rscu)] = 0
  bed_rscu = bed_rscu[complete.cases(bed_rscu),]
  bed_rscu$anticodon2 = as.character(reverseComplement( DNAStringSet(bed_rscu$codon) ))
  #write.csv(bed_rscu,"tRNAseq_anti_ref_codon_full.csv",quote = F,row.names = F)
  
  # cytosolic ribosomal protein genes
  anti_rscu_p = cor.test(log2(bed_rscu$anti_fre+1),bed_rscu$hits)
  anti_rscu_p
  
  p1 <- ggplot(bed_rscu,aes(x = log2(bed_rscu$anti_fre),y = bed_rscu$hits))+
    geom_point(shape = 16,size = 2.5,color = rgb(252,215,24,maxColorValue = 255))+
    labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("Log2 (anticodon counts)")+
    ylab(paste0(i," codon counts"))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
          axis.text = element_text(size=14))
  p1
  ggsave(paste0(i,"_anticodon_cor.pdf"),width = 15, height = 10, units = "cm")
  write.csv(bed_rscu,paste0(i,"_anti_ref_codon.csv"),quote = F,row.names = F)
}


# whole genome genes
whole = read.table("codon_frequency.txt")
names(whole) = c("codon","aa","whole_fre")
bed_whole = merge(b2,whole,by = "codon")
bed_whole = bed_whole[-grep("M",bed_whole$aa),]
bed_whole = bed_whole[-grep("W",bed_whole$aa),]

anti_whole_p = cor.test(log2(bed_whole$anti_fre),bed_whole$whole_fre)
anti_whole_p
p2 <- ggplot(bed_whole,aes(x = log2(bed_whole$anti_fre),y = bed_whole$whole_fre))+
  geom_point(shape = 16,size = 2.5,color = rgb(101,200,208,maxColorValue = 255))+
  labs(title = paste0(spA," cor_anticodon_genome_codon    ","r=",round(anti_whole_p$estimate,5),"  p=",round(anti_whole_p$p.value,5)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Log2 (anticodon counts)")+
  ylab("Whole genome codon counts")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16),
        axis.text = element_text(size = 14))
p2
ggsave("Sc_genome_codon_forPPT.pdf",width = 15,height = 10,units = "cm")
write.csv(bed_rscu,"tRNAseq_anti_genome_codon.txt",quote = F,row.names = F)

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
cor_df4 = cor.test(df4$ribo_hits,log2(df4$anti_fre))
cor_df4
plot(df4$ribo_hits,df4$anti_fre,log = "xy")
