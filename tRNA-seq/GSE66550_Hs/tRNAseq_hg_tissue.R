#==========================================================================================
# Yingying Dong.2019-12-18.Modification data: 2022-4.
# hg38 HEK293 tRNA-seq anticodon & every tissue
#==========================================================================================
library(dplyr)
require(Biostrings)
library(ggplot2)
library(scales)
spA = "Hs"
spe = "hg38"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# GSM1624820 demethylases 82.06%	mapped rate
dir.create("tis_anti_cor")

tRNA = read.table('GSM1624820.txt',sep = '\t',header = T)
name = 'GSM1624820'
names(tRNA) = c('anticodon','anti_fre')
tRNA$anticodon = sub('(e)TCA',"TCA",tRNA$anticodon,fixed = TRUE)  # Ignore regular expressions.

t = tRNA %>% group_by(anticodon) %>% summarise(Frq = sum(anti_fre))
t$codon = as.character(reverseComplement(DNAStringSet(t$anticodon)))

#getwd()
tis_path = "G:/研究生阶段/tRNA/Hs/seq(HEK293)/cellLine_RSCU/"
tissue_array = list.files("G:/研究生阶段/tRNA/Hs/seq(HEK293)/cellLine_RSCU",pattern = "txt$") 
tissue_array
wri_path = "G:/研究生阶段/tRNA/Hs/seq(HEK293)/tis_anti_cor/"
tis <- vector()
corr <- vector()
ls <- list()
for(i in tissue_array){
  if(FALSE){
    tissue = read.table(paste0(tis_path,'adipose tissue_RSCU.txt'),sep = '\t',header = T)
    name = 'adipose tissue'
  }
  tissue = read.table(paste0(tis_path,i), sep = "\t",header = T)
  name = sub("_RSCU.txt", "",i) 
  
  
  tRNA_rscu = merge(t,tissue,by = "codon", all = T)
  #tRNA_rscu = tRNA_rscu[-grep("STOP",tRNA_rscu$AA),]
  #tRNA_rscu = tRNA_rscu[-grep("M",tRNA_rscu$AA),]
  #tRNA_rscu = tRNA_rscu[-grep("W",tRNA_rscu$AA),]

  tRNA_rscu$anticodon = as.character(reverseComplement( DNAStringSet(tRNA_rscu$codon) ))
  tRNA_rscu[is.na(tRNA_rscu)] = 0
  #write.csv(paste0(pathw,tRNA_rscu),"tRNAseq_anti_hE_codon_full.csv",quote = F,row.names = F)
  tRNA_rscu[tRNA_rscu == 0] <- NA
  tRNA_rscu = tRNA_rscu[complete.cases(tRNA_rscu),]
  
  anti_rscu_p = cor.test(tRNA_rscu$Frq+1,tRNA_rscu$hits)
  anti_rscu_p
  
  p1 <- ggplot(tRNA_rscu,aes(x = Frq,y = hits))+
    geom_point(shape = 16,size = 2.5,color = "DarkTurquoise")+
    labs(title = paste0(spA," cor_anticodon_codon    ","r=",round(anti_rscu_p$estimate,5),"  p=",round(anti_rscu_p$p.value,5)))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
    theme_bw()+
    xlab("Log2 (anticodon counts)")+
    ylab(paste0(i," HT codon count"))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title.x =element_text(size=16), axis.title.y=element_text(size=16),
          axis.text = element_text(size = 14))
  #p1+geom_text(aes(label=codon), size=4,vjust = 0, nudge_y = 0.1,check_overlap = TRUE)
  p1
  ggsave(paste0(wri_path,name,".pdf"),width = 15, height = 10, units = "cm")
  
  ls[name] <- round(anti_rscu_p$estimate,5)
  
  tis <- append(tis,name)
  corr <- append(corr,round(anti_rscu_p$estimate,5))
}
print(ls)
print(tis)
print(corr)
df <- as.data.frame(ls) # convert list to data.frame
# method2
df2 <- cbind(tis,corr)
class(df2)
df2 <- as.data.frame(df2)
write.table(df2, file = paste0("tissue_correlation.csv"),sep = ",",quote = FALSE,
            row.names = F, col.names = T)
le <- length(df2$tis)
df2_p <- cbind(df2,1:le)
p <- ggplot(df2,aes(x=1:le,y=corr))+
  geom_boxplot()
p

p1 <- ggplot(df2,aes(x = tis,y = corr))+
  geom_point(shape = 16,size = 1,color = "blue")+
  labs(title = "Correlation of HEK293 tRNA-seq anticodon counts and different tissue codon counts")+
  annotation_logticks(sides="bl")+
  stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "black",size = 0.75)+
  theme_bw()+
  xlab("Tissue name")+
  ylab("Correlation of anticodon counts and codon counts")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text = element_text(size=14))
p1
ggsave("DIFtissueCodon_tRNAanticodon.pdf",width = 15, height = 10, units = "cm")
