#========================================================================================================
# 2019-11-11.Modified date:2019-12-4.Author:Yingying Dong.Correlation analysis of modified CAI and 
#  global gene expression in all samples.The weight from hE_hTsome genes.
#========================================================================================================
library(getopt)
library(ggplot2)
library(MASS)
library(scales)
command=matrix(c("species","s",1,"character",
                 "experiment","e",1,"numeric",
                 "species_abb","A",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$species) || is.null(args$experiment) || is.null(args$species_abb)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

species = args$species
exp = args$experiment
speA_exp = args$species_abb

if(FALSE){
  species = "C_elegans_Ensl_WBcel235"
  exp = "1"
  speA_exp = "Ce"
}

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/", species, "/ref/",speA_exp,"_mCAI_ribo.txt"), 
                  sep = "\t", header = T,quote = "",fill = T)
cai$gene_id = gsub(">gene-", "",cai$gene_id)
dir.create(paste0("picture_bymCAI_ribo",exp))
dir.create(paste0("correlation_bymCAI_ribo",exp))
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment",exp,"/"))
gtf_array = list.files(getwd(),pattern = "[SE]RR\\d.+out$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  if (FALSE){
    gtf = read.table("SRR6815557_abund.out", sep = "\t", header=T,quote = "")
    name = "SRR6815557"
  }
  
  gtf = read.table(i, sep = "\t", header = T,quote = "",fill = T)
  name = sub("^([^.]*).*", "\\1",i) 
  name = sub("_abund","",name)
  
  gtf = gtf[,-c(2,3,4,5,6,7)]
  cai_tpm = merge(cai, gtf, by.x = "gene_id", by.y = "Gene.ID", all = T)
  cai_tpm[cai_tpm == 0] <- NA
  cai_tpm = cai_tpm[complete.cases(cai_tpm),] # Or cai_tpm = na.omit(cai_tpm)
  df <- cai_tpm
  q = quantile(df$TPM[df$TPM > 0], prob = seq(0,1,0.01))
  q
  df2 = df[df$TPM >= q[100],] # top 1%
  df3 = df[df$`TPM` >= q[90],] # top 10%
  df_l = df[df$TPM <= q[48],] # low 48%
  df4 = df[df$TPM >= q[100],] # top 1%
  df5 = df[df$TPM >= q[48],] # top 10%
  df_def = subset(df5, !df5$TPM %in%c(df4$TPM)) # Remove the a data frame from the b data frame.
  df6 = df[df$`TPM` >= q[50],] # top 50%
  df_50 = subset(df6, !df6$TPM %in%c(df4$TPM))
  #===========================================================================
  # For all gene correlation mCAI and TPM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/",name,"_CAI_TPM.svg"))
  CAI_cor = cor(cai_tpm$mCAI_value,cai_tpm$TPM)
  CAI_cor
  CAI_cor_p = cor.test(cai_tpm$mCAI_value,cai_tpm$TPM)

  p <- ggplot(cai_tpm,aes(x = cai_tpm$mCAI_value ,y = cai_tpm$TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_mCAI_TPM    ","r=",round(CAI_cor,3),"  p=",round(CAI_cor_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p)
  
  dev.off()
  #===========================================================================
  # For Top2%Gene correlation mCAI and TPM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","TOP1_",name,"_CAI_TPM.svg"))
  top_cai_cor = cor(df2$mCAI_value,df2$TPM)
  top_cai_p = cor.test(df2$mCAI_value,df2$TPM)
  p2 <- ggplot(df2,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_TOP1_mCAI_TPM    ","r=",round(top_cai_cor,3),"  p=",round(top_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p2)
  dev.off()
  #===========================================================================
  # For Top10%Gene correlation mCAI and TPM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","h_",name,"_CAI_TPM.svg"))
  high_cai_cor = cor(df3$mCAI_value,df3$TPM)
  high_cai_p = cor.test(df3$mCAI_value,df3$TPM)
  p3 <- ggplot(df3,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_h_mCAI_TPM    ","r=",round(high_cai_cor,3),"  p=",round(high_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p3)
  dev.off()
  #==========================================================================
  # For 0-48%Gene correlation mCAI and TPM 
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","l48_",name,"_CAI_TPM.svg"))
  low_cai_cor = cor(df_l$mCAI_value,df_l$TPM)
  low_cai_p = cor.test(df_l$mCAI_value,df_l$TPM)
  p4 <- ggplot(df_l,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_l48_mCAI_TPM    ","r=",round(low_cai_cor,3),"  p=",round(low_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p4)
  dev.off()
  #===========================================================================
  # For defGene correlation mCAI and TPM .def: 48-99%
  #===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","def_",name,"_CAI_TPM.svg"))
  def_cai_cor = cor(df_def$mCAI_value,df_def$TPM)
  def_cai_p = cor.test(df_def$mCAI_value,df_def$TPM)
  def_cai_p
  p5 <- ggplot(df_def,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_def_mCAI_TPM    ","r=",round(def_cai_cor,3),"  p=",round(def_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p5
  print(p5)
  dev.off()
  #===============================================================================
  # For 50-100%Gene correlation mCAI and TPM 
  #===============================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","half_",name,"_CAI_TPM.svg"))
  half_cai_cor = cor(df6$mCAI_value,df6$TPM)
  half_cai_p = cor.test(df6$mCAI_value,df6$TPM)
  p6 <- ggplot(df6,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_half_mCAI_TPM    ","r=",round(half_cai_cor,3),"  p=",round(half_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p6)
  dev.off()
  #===============================================================================
  # For 50-99%Gene correlation mCAI and TPM 
  #===============================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_ribo",exp,"/","test50_",name,"_CAI_TPM.svg"))
  test_cai_cor = cor(df_50$mCAI_value,df_50$TPM)
  test_cai_p = cor.test(df_50$mCAI_value,df_50$TPM)
  p7 <- ggplot(df_50,aes(x = mCAI_value ,y = TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_50-99_mCAI_TPM    ","r=",round(test_cai_cor,3),"  p=",round(test_cai_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p7)
  dev.off()
  #===============================================================================
  # Write out data,convenient to calculate the average
  #===============================================================================
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,test_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bymCAI_ribo",exp,"/",name,"correlation_mCAI"),
              quote = FALSE,row.names = "correlation", 
              col.names = "all_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor\ttest_cor")
  write.table(paste(CAI_cor,top_cai_cor,high_cai_cor,low_cai_cor,def_cai_cor,half_cai_cor,test_cai_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bymCAI_ribo",exp,"/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bymCAI_ribo",exp,"/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
names(a) = c("all_CAI_cor","top_CAI_cor","high_CAI_cor","low_CAI_cor","def_CAI_cor","half_CAI_cor","test_cor")
write.table(paste(mean(a$all_CAI_cor),mean(a$top_CAI_cor),mean(a$high_CAI_cor),
                  mean(a$low_CAI_cor),mean(a$def_CAI_cor),mean(a$half_CAI_cor),mean(a$test_cor),sep = '\t'),file = "mean_cor.txt",
            quote = FALSE,row.names = "mean_mCAI_correlation", 
            col.names = "\tall_CAI_cor\ttop_CAI_cor\thigh_CAI_cor\tlow_CAI_cor\tdef_CAI_cor\thalf_CAI_cor\ttest_cor")
