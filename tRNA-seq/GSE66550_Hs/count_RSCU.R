#=================================================================================================
# 2019-11-25.Author:Yingying Dong.
#=================================================================================================
library(ggplot2)
spA = 'Hs'
species = 'hg19'
rscu_path = paste0('/media/hp/disk1/DYY/reference/annotation/',species,'/')

setwd(paste0("/media/hp/disk2/DYY2/tRNA/",spA,"/seq/"))
count_array = list.files(getwd(),pattern = ".txt$") 
count_array
for (i in count_array) {
  if(FALSE){
    tRNA_count = read.table('GSM1624818.txt',header = F)
    name = 'GSM1624818'
  }
  tRNA_count = read.table(i,header = T, quote = "")
  name = sub("^([^.]*).*", "\\1",i) 
  
  names(tRNA_count) = c("aa","count")
  tRNA_aa = aggregate(tRNA_count$count, by=list(Category=tRNA_count$aa), FUN=sum)
  names(tRNA_aa) = c("aa","acc_count")
  
  aa_AA = read.table(paste0("/media/hp/disk2/DYY2/tRNA/aa_AA.txt"),header = T)
  
  RSCU = read.table(paste0(rscu_path,'ribo_codon_fre_RSCU.txt'),header = T)
  hits_acc = aggregate(RSCU$hits,by=list(Category=RSCU$AA),FUN=sum)
  RSCU_acc = aggregate(RSCU$RSCU,by=list(Category=RSCU$AA),FUN=sum)
  names(hits_acc) = c("AA","acc_hits")
  
  tRNA_AA = merge(tRNA_aa,aa_AA,by = "aa",all = T)
  tRNA_ribo = merge(tRNA_AA,hits_acc,by = 'AA',all = T)
  tRNA_ribo_cor = cor(tRNA_ribo$acc_count,tRNA_ribo$acc_hits)
  tRNA_ribo_cor
  plot(tRNA_ribo,df$FPKM,log = "y",main =paste0(name,"_all_Nc_FPKM  ",all_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  
  p <- ggplot(tRNA_ribo,aes(x = tRNA_ribo$acc_count ,y = tRNA_ribo$acc_hits))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_tRNA_ribo    ","r=",tRNA_ribo_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("tRNA AA count")+
    ylab("ribosomal protein AA count")
  p
  
}
