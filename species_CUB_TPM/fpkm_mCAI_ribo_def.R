#========================================================================================================
# 2019-11-14.Author:Yingying Dong.Correlation analysis of modified CAI and 
#  global gene expression in all samples.The gene expression data quantified by other laboratories.
#========================================================================================================
library(getopt)
library(ggplot2)
library(MASS)
library(scales)

species = "Komagataella_pastoris"
speA_exp = "Kp2"

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/", species, "/ref/",speA_exp,"_mCAI.txt"), 
                  sep = "\t", header = T,quote = "",fill = T)
cai$gene_id = gsub(">gene-", "",cai$gene_id)
dir.create(paste0("picture_bymCAI_ribo"))
dir.create(paste0("correlation_bymCAI_ribo"))
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/"))
gtf_array = list.files(getwd(),pattern = ".+csv$") 
gtf_array

#graphics.off()
for (i in gtf_array) {
  if (FALSE){
    fpkm = read.csv("W_D_1.csv", header=T,quote = "")
    name = "W_D_1"
  }
  fpkm = read.csv(i,header = T,fill = T)
  name = sub("^([^.]*).*", "\\1",i) 
  
  cai_FPKM = merge(cai, fpkm, by.x = "gene_id", by.y = "Locus.ID", all = T)
  cai_FPKM[cai_FPKM == 0] <- NA
  cai_FPKM = cai_FPKM[complete.cases(cai_FPKM),]
  
  df <- cai_FPKM
  q = quantile(df$FPKM[df$FPKM > 0], prob = seq(0,1,0.01))
  q
  df2 = df[df$FPKM >= q[99],] # top 2%
  df3 = df[df$`FPKM` >= q[90],] # top 10%
  df_l = df[df$FPKM <= q[48],] # low 48%
  df4 = df[df$FPKM >= q[100],] # top 1%
  df5 = df[df$FPKM >= q[48],] # top 10%
  df_def = subset(df5, !df5$FPKM %in%c(df4$FPKM)) # Remove the a data frame from the b data frame.
  df6 = df[df$`FPKM` >= q[50],] # top 50%
  #===========================================================================
  # For all gene correlation mCAI and FPKM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/",name,"_CAI_FPKM.svg"))
  CAI_cor = cor(cai_FPKM$mCAI_value,cai_FPKM$FPKM)
  CAI_cor
  p <- ggplot(cai_FPKM,aes(x = cai_FPKM$mCAI_value ,y = cai_FPKM$FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_mCAI_FPKM    ","r=",CAI_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p)
  dev.off()
  #===========================================================================
  # For Top2%Gene correlation mCAI and FPKM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/","TOP1_",name,"_CAI_FPKM.svg"))
  top_cai_cor = cor(df2$mCAI_value,df2$FPKM)
  p2 <- ggplot(df2,aes(x = mCAI_value ,y = FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_TOP1_mCAI_FPKM    ","r=",top_cai_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p2)
  dev.off()
  #===========================================================================
  # For Top10%Gene correlation mCAI and FPKM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/","h_",name,"_CAI_FPKM.svg"))
  high_cai_cor = cor(df3$mCAI_value,df3$FPKM)
  p3 <- ggplot(df3,aes(x = mCAI_value ,y = FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_h_mCAI_FPKM    ","r=",high_cai_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p3)
  dev.off()
  #==========================================================================
  # For 0-48%Gene correlation mCAI and FPKM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/","l48_",name,"_CAI_FPKM.svg"))
  low_cai_cor = cor(df_l$mCAI_value,df_l$FPKM)
  p4 <- ggplot(df_l,aes(x = mCAI_value ,y = FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_l48_mCAI_FPKM    ","r=",low_cai_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p4)
  dev.off()
  #===========================================================================
  # For defGene correlation mCAI and FPKM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/","def_",name,"_CAI_FPKM.svg"))
  def_cai_cor = cor(df_def$mCAI_value,df_def$FPKM)
  p5 <- ggplot(df_def,aes(x = mCAI_value ,y = FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_def_mCAI_FPKM    ","r=",def_cai_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p5)
  dev.off()
  #===============================================================================
  # For 50-100%Gene correlation mCAI and FPKM 
  #===============================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo/","half_",name,"_CAI_FPKM.svg"))
  half_cai_cor = cor(df6$mCAI_value,df6$FPKM)
  p6 <- ggplot(df6,aes(x = mCAI_value ,y = FPKM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_half_mCAI_FPKM    ","r=",half_cai_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq FPKM")
  print(p6)
  dev.off()
  #===============================================================================
  # Write out data,convenient to calculate the average
  #===============================================================================
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bymCAI_ribo/",name,"correlation_mCAI"),
              quote = FALSE,row.names = "correlation", 
              col.names = "all_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor")
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bymCAI_ribo/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bymCAI_ribo/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
names(a) = c("all_CAI_cor","top_CAI_cor","high_CAI_cor","low_CAI_cor","def_CAI_cor","half_CAI_cor")
write.table(paste(mean(a$all_CAI_cor),mean(a$top_CAI_cor),mean(a$high_CAI_cor),
                  mean(a$low_CAI_cor),mean(a$def_CAI_cor),mean(a$half_CAI_cor),sep = '\t'),file = "mean_cor.txt",
            quote = FALSE,row.names = "mean_mCAI_correlation", 
            col.names = "\tall_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor")
