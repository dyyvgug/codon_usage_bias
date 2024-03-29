#==========================================================================================
# Yingying Dong.2019-12-18.tRNA seq anticodon & the reference set codon
#==========================================================================================
library(dplyr)
require(Biostrings)
library(ggplot2)
spA = "Hs"
spe = "hg38"
setwd("G:\\tRNA\\Hs\\seq(HEK293)\\ribo")
pathw = paste0("G:\\tRNA\\Hs\\seq(HEK293)\\ribo")
tRNA_array = list.files(getwd(),pattern = "txt$") 
tRNA_array
for(i in tRNA_array){
  if(FALSE){
    tRNA = read.table('GSM1624822.txt',sep = '\t',header = T)
    name = 'GSM1624822.txt'
  }
  tRNA = read.table(i, sep = "\t",header = T)
  name = sub("^([^.]*).*", "\\1",i) 
  
  names(tRNA) = c('anticodon','anti_fre')
  tRNA$anticodon = sub('(e)TCA',"TCA",tRNA$anticodon,fixed = TRUE)  # Ignore regular expressions.
  
  t = tRNA %>% group_by(anticodon) %>% summarise(Frq = sum(anti_fre))
  t$codon = as.character(reverseComplement(DNAStringSet(t$anticodon)))
  
  rscu = read.table(paste0("ribo_codon_fre_RSCU.txt"),
                    header = T, stringsAsFactors = F)
  
  tRNA_rscu = merge(t,rscu,by = "codon", all = T)
  tRNA_rscu = tRNA_rscu[-grep("STOP",tRNA_rscu$AA),]
  tRNA_rscu = tRNA_rscu[-grep("M",tRNA_rscu$AA),]
  tRNA_rscu = tRNA_rscu[-grep("W",tRNA_rscu$AA),]
  #tRNA_rscu = tRNA_rscu[-grep('-',tRNA_rscu$codon),]
  tRNA_rscu$anticodon = as.character(reverseComplement( DNAStringSet(tRNA_rscu$codon) ))
  #tRNA_rscu[is.na(tRNA_rscu)] = 0
  #write.csv(paste0(pathw,tRNA_rscu),"tRNAseq_anti_hE_codon_full.csv",quote = F,row.names = F)
  tRNA_rscu[tRNA_rscu == 0] <- NA
  tRNA_rscu = tRNA_rscu[complete.cases(tRNA_rscu),]
  
  anti_rscu_p = cor.test(tRNA_rscu$Frq,tRNA_rscu$hits)
  anti_rscu_p
  
  p1 <- ggplot(tRNA_rscu,aes(x = tRNA_rscu$Frq,y = tRNA_rscu$hits,color = tRNA_rscu$AA))+
    geom_point(shape = 16,size = 2.5)+
    labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("anticodon count")+
    ylab("high expression level genes codon count")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p1+geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  ggsave(paste0(pathw,name,"_anticodon_codon.pdf"))
  write.table(paste0(pathw,tRNA_rscu),"tRNAseq_anti_hE_codon.txt",sep = '\t',quote = F,row.names = F)
  #------------------------remove large discrete value------------------------------------------
  tRNA_rscu = tRNA_rscu[-grep("ATC",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("TAT",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("AGT",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("AAT",tRNA_rscu$codon),]
  
  anti_rscu_p = cor.test(log2(tRNA_rscu$Frq),tRNA_rscu$hits)
  anti_rscu_p
  
  p3 <- ggplot(tRNA_rscu,aes(x = log2(tRNA_rscu$Frq),y = tRNA_rscu$hits,color = tRNA_rscu$AA))+
    geom_point(shape = 16,size = 2.5)+
    labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("Log2 (anticodon count)")+
    ylab("Cytosolic RP codon count")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p3
  ggsave(p3,filename = paste0(pathw,'ribo_codon_re_dis_anticodon.pdf'))
  write.table(tRNA_rscu,paste0(pathw,"ribo_codon_re_dis_anticodon.txt"),sep = '\t',quote = F,row.names = F)
  
  
  whole = read.table("G:\\学习专用\\annotation\\hg19\\ref\\codon_frequency.txt")
  names(whole) = c("codon","aa","whole_fre")
  
  tRNA_whole = merge(t,whole,by = "codon")
  tRNA_whole = tRNA_whole[-grep("*",tRNA_whole$aa,fixed = TRUE),]
  tRNA_whole = tRNA_whole[-grep("M",tRNA_whole$aa),]
  tRNA_whole = tRNA_whole[-grep("W",tRNA_whole$aa),]
  tRNA_whole[tRNA_whole == 0] <- NA
  tRNA_whole = tRNA_whole[complete.cases(tRNA_whole),]
  
  anti_whole_p = cor.test(tRNA_whole$Frq,tRNA_whole$whole_fre)
  anti_whole_p
  p2 <- ggplot(tRNA_whole,aes(x = tRNA_whole$Frq,y = tRNA_whole$whole_fre,color = tRNA_whole$aa))+
    geom_point(shape = 16,size = 2.5)+
    labs(title = paste0(spA," cor_anticodon_genome_codon    ","r=",round(anti_whole_p$estimate,5),"  p=",round(anti_whole_p$p.value,5)))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("anticodon count")+
    ylab("Whole genome codon count")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))+
    geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  p2+geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  ggsave(paste0(pathw,spA,"_anticodon_genome_codon.pdf"))
  write.table(paste0(pathw,tRNA_whole),"tRNAseq_anti_genome_codon.txt",sep = '\t',quote = F,row.names = F)
  
  tRNA_whole = tRNA_whole[-grep("ATC",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("TAT",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("AGT",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("AAT",tRNA_whole$codon),]
  
  
  anti_rscu_p = cor.test(log2(tRNA_whole$Frq),tRNA_whole$whole_fre)
  anti_rscu_p
  
  p4 <- ggplot(tRNA_whole,aes(x = log2(tRNA_whole$Frq),y = tRNA_whole$whole_fre,color = tRNA_whole$aa))+
    geom_point(shape = 16,size = 2.5)+
    labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("log2 (anticodon count)")+
    ylab("Whole genome codon count")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p4
  ggsave(p4,filename = paste0(pathw,'genome_codon_re_dis_anticodon.pdf'))
  write.table(tRNA_rscu,paste0(pathw,"genome_codon_re_dis_anticodon.txt"),sep = '\t',quote = F,row.names = F)
  
  #==================================================================================================
  # Eight major degenerate codon pairs
  #==================================================================================================
  tRNA_rscu$deg = tRNA_rscu$anticodon
  tRNA_rscu$deg = sub("GGC","AGC",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GGG","AGG",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GGT","AGT",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GAC","AAC",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GGA","AGA",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GCG","ACG",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GAG","AAG",tRNA_rscu$deg)
  tRNA_rscu$deg = sub("GAT","AAT",tRNA_rscu$deg)
  df <- data.frame(tRNA_rscu$AA,tRNA_rscu$hits,tRNA_rscu$anti_fre,tRNA_rscu$deg)
  names(df) = c("AA","ribo_hits","anti_fre","deg")
  df2 <- df %>% group_by(deg) %>% summarise_at(vars(ribo_hits,anti_fre), list(name = sum)) 
  df2$codon = as.character(reverseComplement(DNAStringSet(df2$deg)))
  df3 <- data.frame(tRNA_rscu$codon,tRNA_rscu$AA)
  names(df3) = c("codon","AA")
  df4 = merge(df2,df3,by = "codon",all.x = TRUE)
  names(df4) = c("codon","deg","ribo_hits","anti_fre","AA")
  cor_df4 = cor.test(df4$ribo_hits,log2(df4$anti_fre))
  cor_df4
  plot(df4$ribo_hits,df4$anti_fre,log = "xy")
}
