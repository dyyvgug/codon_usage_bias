#2019-4-12.Correlation analysis of CAI and global gene expression in all samples.DYY.
species = "Escherichia_coli"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_bycodonW.txt", sep = ""), sep = "\t", header = T)
dir.create("picture_bycodonW")
dir.create("correlation_bycodonW")
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  #myfile <- "path1/path2/myoutput.txt"
  #gtf = read.table("SRR7866023.gtf.T", sep = "\t", header=F)
  #gtf <- "SRR7866023.gtf.T"
  gtf = read.table(i, sep = "\t", header=F)
  #name = sub(".gtf.T","", gtf)
  name = sub("^([^.]*).*", "\\1",i) 
  #name
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df = df [-16]
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
#===========================================================================
# Correlation_GC_TPM 
#===========================================================================
  GC_cor = cor(df[,"GC"],df[,"TPM"])
  GC_cor
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/",name,"_GC_TPM.jpg"))
  
  plot(df$GC,df$TPM,log = "y",main = paste0(name,"_GC_TPM  ",GC_cor),
     xlab="GC",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/",name,"_CAI_TPM.jpg"))
  CAI_cor = cor(df$CAI,df$TPM)
  CAI_cor
  plot(df$CAI,df[,"TPM"],log = "y",main = paste0(name,"_CAI_TPM  ",CAI_cor),
     xlab="CAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CBI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/",name,"_CBI_TPM.jpg"))
  CBI_cor = cor(df[,"CBI"],df[,"TPM"])
  CBI_cor
  plot(df[,"CBI"],df[,"TPM"],log = "y",main =paste0(name,"_CBI_TPM  ",CBI_cor),
     xlab="CBI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#==========================================================================
# Correlation_Nc_TPM
#===========================================================================
  df[df == '*****'] <-NA
  df = df[complete.cases(df),]
  #sapply(df, class)
  #sapply(df$Nc, is.factor)
  #mean(as.numeric(as.character(df$Nc)))
  df$Nc = as.numeric(as.character(df$Nc))
  
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/",name,"_Nc_TPM.jpg"))
  Nc_cor = cor(df$Nc,df[,"TPM"])
  Nc_cor
  plot(df$Nc,df[,"TPM"],log = "y",main =paste0(name,"_Nc_TPM  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================   
  write.table(paste(GC_cor,CAI_cor,CBI_cor,Nc_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation/",name,"correlation_bycodonW"),
            quote = FALSE,row.names = "correlation", 
            col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM\tNc_TPM")
  write.table(paste(GC_cor,CAI_cor,CBI_cor,Nc_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_bycodonW/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  
  }

