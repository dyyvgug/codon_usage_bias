#2019-4-12.Correlation analysis of CAI and global gene expression in all samples.DYY.
species = "Apis_mellifera"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_bycodonW.txt", sep = ""), 
                      sep = "\t", header = T,quote = "")
dir.create("picture_Nc1")
dir.create("correlation_Nc1")
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  gtf = read.table("SRR6833960.gtf.T", sep = "\t", header=F,quote = "") # Ignore quotes
  gtf = read.table(i, sep = "\t", header=F,quote = "")
  #name = sub(".gtf.T","", gtf)
  name = sub("^([^.]*).*", "\\1",i) 
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df = df [-16]
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
  df[df == '*****'] <-NA
  df = df[complete.cases(df),]
  #sapply(df, class)
  #sapply(df$Nc, is.factor)
  #mean(as.numeric(as.character(df$Nc)))
  df$Nc = as.numeric(as.character(df$Nc))
  q = quantile(df$TPM[df$TPM > 0], prob = seq(0,1,0.01))
  q
  df2 = df[df$TPM >= q[99],] # top 2%
  df3 = df[df$`TPM` >= q[90],] # top 10%
  df_l = df[df$TPM <= q[48],] # low 48%
  df4 = df[df$TPM >= q[100],] # top 1%
  df5 = df[df$TPM >= q[48],] # top 10%
  df_def = subset(df5, !df5$TPM %in%c(df4$TPM)) # Remove the a data frame from the b data frame.
#===========================================================================
# For allGene correlation Nc and TPM 
#===========================================================================
  all_Nc_cor = cor(df$Nc,df$TPM)
  all_Nc_cor
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/",name,"all_Nc_TPM.jpg"))
  plot(df$Nc,df$TPM,log = "y",main =paste0(name,"_all_Nc_TPM  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For Top2%Gene correlation Nc and TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","TOP1_",name,"_Nc_TPM.jpg"))
  top_Nc_cor = cor(df2$Nc,df2$TPM)
  top_Nc_cor
  plot(df2$Nc,df2$TPM,log = "y",main =paste0(name,"_Nc_TPM_1%  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For Top10%Gene correlation Nc and TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","h_",name,"_Nc_TPM.jpg"))
  high_Nc_cor = cor(df3$Nc,df3[,"TPM"])
  high_Nc_cor
  plot(df3$Nc,df3[,"TPM"],log = "y",main =paste0(name,"_Nc_TPM_h  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#==========================================================================
# For 0-48%Gene correlation Nc and TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","l48_",name,"_Nc_TPM.jpg"))
  low_Nc_cor = cor(df_l$Nc,df_l[,"TPM"])
  low_Nc_cor
  plot(df_l$Nc,df_l[,"TPM"],log = "y",main =paste0(name,"_Nc_TPM_l48  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For defGene correlation Nc and TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","def_",name,"_Nc_TPM.jpg"))
  def_Nc_cor = cor(df_def$Nc,df_def$TPM)
  def_Nc_cor
  plot(df_def$Nc,df_def$TPM,log = "y",main =paste0(name,"_Nc_TPM_def  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================   
  write.table(paste(all_Nc_cor,top_Nc_cor,high_Nc_cor,low_Nc_cor,def_Nc_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_Nc1/",name,"correlation_bycodonW"),
            quote = FALSE,row.names = "correlation", 
            col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM\tNc_TPM")
  write.table(paste(all_Nc_cor,top_Nc_cor,high_Nc_cor,low_Nc_cor,def_Nc_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_Nc1/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  }

