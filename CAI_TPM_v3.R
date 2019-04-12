#2019-4-12.Correlation analysis of CAI and global gene expression in all samples.DYY.
species = "Danio_rerio"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI.txt", sep = ""), sep = "\t", header = T)
dir.create("picture")
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment1/"))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
for (i in gtf_array) {
  gtf = read.table(i, sep = "\t", header=F)
  name = sub(".gtf.T", "", i)
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
#===========================================================================
# Correlation_GC_TPM 
#===========================================================================
  GC_cor = cor(df[,"GC"],df[,"TPM"])
  GC_cor
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture/",name,"_GC_TPM.jpg"))
  plot(df[,"GC"],df[,"TPM"],log = "y",main = paste0(name,"_GC_TPM"," GC_cor =",GC_cor),
     xlab="GC",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture/",name,"_CAI_TPM.jpg"))
  CAI_cor = cor(df$avg_CAI,df$TPM)
  CAI_cor
  plot(df$avg_CAI,df[,"TPM"],log = "y",main = paste0(name,"_CAI_TPM"," CAI_cor =",CAI_cor),
     xlab="CAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CBI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture/",name,"_CBI_TPM.jpg"))
  CBI_cor = cor(df[,"CBI"],df[,"TPM"])
  CBI_cor
  plot(df[,"CBI"],df[,"TPM"],log = "y",main =paste0(name,"_CBI_TPM"," CBI_cor =",CBI_cor),
     xlab="CBI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()

  write.table(paste(GC_cor,CAI_cor,CBI_cor,sep = '\t'),file = paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/",name,"correlation"),
            quote = FALSE,row.names = "correlation", col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM")
}

