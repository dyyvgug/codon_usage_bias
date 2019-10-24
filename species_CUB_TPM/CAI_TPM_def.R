#2019-4-23.Correlation analysis of CAI and high gene expression in all samples.DYY.
species = "Saccharomyces_cerevisiae"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_bycodonW.txt", sep = ""), sep = "\t", header = T)
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
for (i in gtf_array) {
  gtf = read.table(i, sep = "\t", header=F)
  #gtf = read.table("SRR2182451.gtf.T", sep = "\t", header=F)
  name = sub("^([^.]*).*", "\\1",i) 
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df = df [-16]
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
q = quantile(df$`TPM`[df$`TPM`> 0], prob = seq(0,1,0.01))
q
df2 = df[df$TPM >= q[100],] # top 1%
df3 = df[df$TPM >= q[48],] # top 10%
df_def = subset(df3, !df3$TPM %in%c(df2$TPM)) # Remove the a data frame from the b data frame.

#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","def_",name,"_CAI_TPM.jpg"))
  CAI_cor = cor(df_def$CAI,df_def$TPM)
  CAI_cor
  plot(df_def$CAI,df_def$TPM,log = "y",main = paste0(name,"_CAI_TPM_def  ",CAI_cor),
     xlab="CAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CBI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","def_",name,"_CBI_TPM.jpg"))
  CBI_cor = cor(df_def[,"CBI"],df_def[,"TPM"])
  CBI_cor
  plot(df_def[,"CBI"],df_def[,"TPM"],log = "y",main =paste0(name,"_CBI_TPM_def  ",CBI_cor),
     xlab="CBI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_Nc_TPM 
#===========================================================================  
  df_def[df_def == '*****'] <-NA
  df_def = df_def[complete.cases(df_def),]
  #sapply(df, class)
  #sapply(df$Nc, is.factor)
  #mean(as.numeric(as.character(df$Nc)))
  df_def$Nc = as.numeric(as.character(df_def$Nc)) #It is convenient to calculate the NC column from factor type to number type.
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","def_",name,"_Nc_TPM.jpg"))
  Nc_cor = cor(df_def$Nc,df_def[,"TPM"])
  Nc_cor
  plot(df_def$Nc,df_def[,"TPM"],log = "y",main =paste0(name,"_Nc_TPM_def  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================
  write.table(paste(CAI_cor,CBI_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bycodonW/",name,"_def_correlation"),
            quote = FALSE,row.names = "correlation", col.names = "\tCAI_TPM\tCBI_TPM")
  write.table(paste(CAI_cor,CBI_cor,Nc_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_bycodonW/","def_allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  
}

