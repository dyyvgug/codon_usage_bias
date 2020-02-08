#==========================================================================================
# Yingying Dong.2019-12-18.tRNA seq anticodon & the reference set codon
#==========================================================================================
library(dplyr)
require(Biostrings)
library(ggplot2)
library(scales)
spA = "Mm"
spe = "mouse_mm10"
setwd(paste0("G:\\学习专用\\tRNA\\Mm\\seq\\anti_codon"))
tRNA_array = list.files(getwd(),pattern = ".txt$")
tRNA_array
for (i in tRNA_array) {
  if(FALSE){
    tRNA = read.table('GSM2323453_wt1Count.txt',header = F,sep = '\t')
  }
  tRNA = read.table(i,header = F,sep = '\t')
  
  names(tRNA) = c('anticodon','anti_fre')
  tRNA$anticodon = sub('(e)TCA',"TCA",tRNA$anticodon,fixed = TRUE)  # Ignore regular expressions.
  
  t = tRNA %>% group_by(anticodon) %>% summarise(Frq = sum(anti_fre))
  t$codon = as.character(reverseComplement(DNAStringSet(t$anticodon)))
  
  rscu = read.table("G:\\学习专用\\annotation\\mouse_mm10\\ribo_codon_fre_RSCU.txt",
                    header = T, stringsAsFactors = F)
  
  tRNA_rscu = merge(t,rscu,by = "codon", all = T)
  tRNA_rscu = tRNA_rscu[-grep("STOP",tRNA_rscu$AA),]
  tRNA_rscu = tRNA_rscu[-grep("M",tRNA_rscu$AA),]
  tRNA_rscu = tRNA_rscu[-grep("W",tRNA_rscu$AA),]
  #tRNA_rscu[is.na(tRNA_rscu)] = 0
  tRNA_rscu = tRNA_rscu[complete.cases(tRNA_rscu),]
  
  anti_rscu_p = cor.test(log2(tRNA_rscu$Frq+1),tRNA_rscu$hits)
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
  ylab("Reference set codon count")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p1
  p1+geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  ggsave(paste0(spA,"_anticodon_codon.pdf"))
  write.table(tRNA_rscu,"tRNAseq_anti_ref_codon.txt",sep = '\t',quote = F,row.names = F)
  
  #------------------------remove large discrete value------------------------------------------
  tRNA_rscu = tRNA_rscu[-grep("ATC",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("TCC",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("CAT",tRNA_rscu$codon),]
  tRNA_rscu = tRNA_rscu[-grep("GTC",tRNA_rscu$codon),]
  
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
    xlab("log2 (anticodon count)")+
    ylab("Cytosolic RP codon count)")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p3
  ggsave(p3,filename = paste0('ribo_codon_re_dis_anticodon.pdf'))
  write.table(tRNA_rscu,paste0(pathw,"hE_codon_re_dis_anticodon.txt"),sep = '\t',quote = F,row.names = F)
  
  whole = read.table("G:\\学习专用\\annotation\\mouse_mm10\\ref\\codon_frequency.txt")
  names(whole) = c("codon","aa","whole_fre")
  
  tRNA_whole = merge(t,whole,by = "codon")
  tRNA_whole = tRNA_whole[-grep("*",tRNA_whole$aa,fixed = TRUE),]
  tRNA_whole = tRNA_whole[-grep("M",tRNA_whole$aa),]
  tRNA_whole = tRNA_whole[-grep("W",tRNA_whole$aa),]
  
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
  axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p2+geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  ggsave(paste0(spA,"_anticodon_genome_codon.pdf"))
  write.table(tRNA_whole,"tRNAseq_anti_genome_codon.txt",sep = '\t',quote = F,row.names = F)
  
  tRNA_whole = tRNA_whole[-grep("ATC",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("TCC",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("CAT",tRNA_whole$codon),]
  tRNA_whole = tRNA_whole[-grep("GTC",tRNA_whole$codon),]
  
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
    xlab("Log2 (anticodon count)")+
    ylab("Whole genome codon count)")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))
  p4
  ggsave(p4,filename = paste0(pathw,'genome_codon_re_dis_anticodon.pdf'))
  write.table(tRNA_rscu,paste0(pathw,"genome_codon_re_dis_anticodon.txt"),sep = '\t',quote = F,row.names = F)
}
