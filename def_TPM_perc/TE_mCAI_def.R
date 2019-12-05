#========================================================================================================
# 2019-12-5.Author:Yingying Dong.Correlation analysis of modified CAI and TE(ribosome number) in all 
#  samples.The weight from ribosomal protein genes.
#========================================================================================================
library(getopt)
library(ggplot2)
library(MASS)
library(scales)

species = "C_elegans_Ensl_WBcel235"
exp = "2"
speA_exp = "Ce"
RNAseq = "/SRR7160566_abund.out"

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/", species, "/ref/",speA_exp,"_mCAI_ribo.txt"), 
                  sep = "\t", header = T,quote = "",fill = T)
cai$gene_id = gsub(">gene-", "",cai$gene_id)
RNA = read.table(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment3",RNAseq),sep = "\t",header = T,quote = "")
RNA = RNA[,-c(3,4,5,6,7)]
dir.create(paste0("picture_TE_bymCAI",exp))
dir.create(paste0("correlation_TE_bymCAI",exp))
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment",exp,"/aligned/"))
ribo_array = list.files(getwd(),pattern = ".out$")
ribo_array
for (i in ribo_array){
  if(FALSE) # examination
  { 
    ribo = read.table("SRR1804340_abund.out",sep = "\t",header = T,quote = "")
    name = "SRR1804340_abund.out"
    name = sub("^([^.]*).*", "\\1",name)
    name = gsub("_abund","",name)
  }
  ribo = read.table(i,sep = "\t",header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i)
  name = sub("_abund","",name)
  
  ribo = ribo[,-c(2,3,4,5,6,7)]
  names(ribo) = c("Gene.ID","ribo_FPKM","ribo_TPM")
  RNA_ribo = merge(RNA,ribo,by = "Gene.ID",all = T)
  RNA_ribo[RNA_ribo == 0] <-NA
  RNA_ribo = RNA_ribo[complete.cases(RNA_ribo),]
  #RNA_ribo_num = cbind(RNA_ribo,round(RNA_ribo$TPM / RNA_ribo$ribo_TPM,3))
  ribo_num = round(RNA_ribo$ribo_TPM / RNA_ribo$TPM,3)
  RNA_ribo_num = cbind(RNA_ribo,ribo_num)
  RNA_ribo_num = RNA_ribo_num[order(RNA_ribo_num$ribo_num,decreasing = T),]
  #RNA_ribo_num$Gene.ID = gsub("_.*", "", RNA_ribo_num[,1])
  write.table(RNA_ribo_num,file = paste0("./ribo_num/",name,"_riboNum.txt"),
              sep = "\t",quote = F,row.names = F)
  #=================================================================================================================
  # Extract genes encoding only proteins 
  #=================================================================================================================
  only_protein = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CBI_CAI.txt"),header = T)
  only_protein = only_protein[,-c(3,4,5)]
  only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "transcription_id",all = T)
  only_protein_num = only_protein_num[complete.cases(only_protein_num),]
  only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]

  q = quantile(only_protein_num$ribo_num,probs = seq(0,1,0.01))
  q
  hiTE = only_protein_num[only_protein_num$ribo_num >= q[99],]
  
  cai_tpm = merge(cai,only_protein_num, by.x = "gene_id", by.y = "Gene.ID", all = T)
  cai_tpm[cai_tpm == 0] <- NA
  cai_tpm = cai_tpm[complete.cases(cai_tpm),] # Or cai_tpm = na.omit(cai_tpm)
  df <- cai_tpm
  #df$ribo_num[df$ribo_num<1] <- NA
  #df = df[complete.cases(df),]
  
  q = quantile(df$ribo_num[df$ribo_num > 0], prob = seq(0,1,0.01))
  q
  df2 = df[df$ribo_num >= q[100],] # top 1%
  df3 = df[df$ribo_num >= q[90],] # top 10%
  df_l = df[df$ribo_num <= q[50],] # low 50%
  df4 = df[df$ribo_num >= q[100],] # top 1%
  df5 = df[df$ribo_num >= q[48],] # top 10%
  df_def = subset(df5, !df5$ribo_num %in%c(df4$ribo_num)) # Remove the a data frame from the b data frame.
  df6 = df[df$ribo_num >= q[50],] # top 50%
  df_50 = subset(df6, !df6$ribo_num %in%c(df4$ribo_num))
  #===========================================================================
  # For all gene correlation mCAI and ribo_num 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/",name,"_CAI_ribo_num.svg"))
  CAI_cor = cor(df$mCAI_value,df$ribo_num)
  CAI_cor
  CAI_cor_p = cor.test(df$mCAI_value,df$ribo_num)
  CAI_cor_p
  p <- ggplot(df,aes(x = df$mCAI_value ,y = df$ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_mCAI_TPM    ","r=",round(CAI_cor,3),"  p=",round(CAI_cor_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p)
  dev.off()
  #===========================================================================
  # For Top2%Gene correlation mCAI and ribo_num 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","TOP1_",name,"_CAI_ribo_num.svg"))
  top_cai_cor = cor(df2$mCAI_value,df2$ribo_num)
  top_cai_p = cor.test(df2$mCAI_value,df2$ribo_num)
  p2 <- ggplot(df2,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_TOP1_mCAI_ribo_num    ","r=",round(top_cai_cor,3),"  p=",round(top_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p2)
  dev.off()
  #===========================================================================
  # For Top10%Gene correlation mCAI and ribo_num 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","h_",name,"_CAI_ribo_num.svg"))
  high_cai_cor = cor(df3$mCAI_value,df3$ribo_num)
  high_cai_p = cor.test(df3$mCAI_value,df3$ribo_num)
  p3 <- ggplot(df3,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_h_mCAI_ribo_num    ","r=",round(high_cai_cor,3),"  p=",round(high_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p3)
  dev.off()
  #==========================================================================
  # For 0-48%Gene correlation mCAI and ribo_num 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","l50_",name,"_CAI_ribo_num.svg"))
  low_cai_cor = cor(df_l$mCAI_value,df_l$ribo_num)
  low_cai_p = cor.test(df_l$mCAI_value,df_l$ribo_num)
  p4 <- ggplot(df_l,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_l50_mCAI_ribo_num    ","r=",round(low_cai_cor,3),"  p=",round(low_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p4)
  dev.off()
  #===========================================================================
  # For defGene correlation mCAI and ribo_num .def: 48-99%
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","def_",name,"_CAI_ribo_num.svg"))
  def_cai_cor = cor(df_def$mCAI_value,df_def$ribo_num)
  def_cai_p = cor.test(df_def$mCAI_value,df_def$ribo_num)
  def_cai_p
  p5 <- ggplot(df_def,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_def_mCAI_ribo_num    ","r=",round(def_cai_cor,3),"  p=",round(def_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p5
  print(p5)
  dev.off()
  #===============================================================================
  # For 50-100%Gene correlation mCAI and ribo_num 
  #===============================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","half_",name,"_CAI_ribo_num.svg"))
  half_cai_cor = cor(df6$mCAI_value,df6$ribo_num)
  half_cai_p = cor.test(df6$mCAI_value,df6$ribo_num)
  p6 <- ggplot(df6,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_half_mCAI_ribo_num    ","r=",round(half_cai_cor,3),"  p=",round(half_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p6)
  dev.off()
  #===============================================================================
  # For 50-99%Gene correlation mCAI and ribo_num 
  #===============================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_TE_bymCAI",exp,"/","test50_",name,"_CAI_ribo_num.svg"))
  test_cai_cor = cor(df_50$mCAI_value,df_50$ribo_num)
  test_cai_p = cor.test(df_50$mCAI_value,df_50$ribo_num)
  p7 <- ggplot(df_50,aes(x = mCAI_value ,y = ribo_num))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_50-99_mCAI_ribo_num    ","r=",round(test_cai_cor,3),"  p=",round(test_cai_p$p.value,5)))+
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("TE value")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p7)
  dev.off()
  #===============================================================================
  # Write out data,convenient to calculate the average
  #===============================================================================
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,test_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_TE_bymCAI",exp,"/",name,"correlation_mCAI"),
              quote = FALSE,row.names = "correlation", 
              col.names = "all_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor\ttest_cor")
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,test_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_TE_bymCAI",exp,"/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_TE_bymCAI",exp,"/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
names(a) = c("all_CAI_cor","top_CAI_cor","high_CAI_cor","low_CAI_cor","def_CAI_cor","half_CAI_cor","test_cor")
write.table(paste(mean(a$all_CAI_cor),mean(a$top_CAI_cor),mean(a$high_CAI_cor),
                  mean(a$low_CAI_cor),mean(a$def_CAI_cor),mean(a$half_CAI_cor),mean(a$test_cor),sep = '\t'),file = "mean_cor.txt",
            quote = FALSE,row.names = "mean_mCAI_correlation", 
            col.names = "\tall_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor\ttest_cor")
