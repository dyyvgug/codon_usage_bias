#2019-4-12.Correlation analysis of CAI and high gene expression in all samples.DYY.
species = "Escherichia_coli"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_bycodonW.txt", sep = ""), sep = "\t", header = T)
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
for (i in gtf_array) {
  gtf = read.table(i, sep = "\t", header=F)
  #gtf = read.table("SRR7866023.gtf.T", sep = "\t", header=F)
  #gtf <- "SRR7866023.gtf.T"
  #name = sub(".gtf.T", "", i)
  #name = "SRR7866023"
  name = sub("^([^.]*).*", "\\1",i) 
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df = df [-16]
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
q = quantile(df$`TPM`[df$`TPM`> 0], prob = seq(0,1,0.01))
q
df2 = df[df$`TPM` >= q[99],] # top 1%
#===========================================================================
# Correlation_GC_TPM 
#===========================================================================
  GC_cor = cor(df2[,"GC"],df2[,"TPM"])
  GC_cor
  jpeg(file = paste("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","TOP1_",name,"_GC_TPM.jpg",sep = ""))
  plot(df2[,"GC"],df2[,"TPM"],log = "y",main = paste0(name,"_GC_TPM_1%  ",GC_cor),
     xlab="GC",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
#library(ggplot2)
#DF <- data.frame(df2[,"GC"], df2[,"TPM"])
#p <- ggplot(df,aes(x=degree,y=complex) + geom_point(shape=19)) +
#  xlab("Degree") + ylab("Number of complexes") + 
#  geom_smooth(method = lm)
#ggplot(DF, aes(x = df2[,"GC"], y = df2[,"TPM"])) + geom_point() +
#  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE) +
#  stat_smooth(method = 'lm', formula = y ~ poly(x,2), aes(colour = 'polynomial'), se= FALSE) +
#  stat_smooth(method = 'nls', formula = y ~ a * log(x) +b, aes(colour = 'logarithmic'), se = FALSE, start = list(a=1,b=1)) +
#  stat_smooth(method = 'nls', formula = y ~ a*exp(b *x), aes(colour = 'Exponential'), se = FALSE, start = list(a=1,b=1)) +
#  theme_bw() +
#  scale_colour_brewer(name = 'Trendline', palette = 'Set2')''
  dev.off()
#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","TOP1_",name,"_CAI_TPM.jpg"))
  CAI_cor = cor(df2$CAI,df2$TPM)
  CAI_cor
  plot(df2$CAI,df2[,"TPM"],log = "y",main = paste0(name,"_CAI_TPM_1%  ",CAI_cor),
     xlab="CAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_CBI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","TOP1_",name,"_CBI_TPM.jpg"))
  CBI_cor = cor(df2[,"CBI"],df2[,"TPM"])
  CBI_cor
  plot(df2[,"CBI"],df2[,"TPM"],log = "y",main =paste0(name,"_CBI_TPM_1%  ",CBI_cor),
     xlab="CBI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# Correlation_Nc_TPM 
#===========================================================================  
  df2[df2 == '*****'] <-NA
  df2 = df2[complete.cases(df2),]
  #sapply(df, class)
  #sapply(df$Nc, is.factor)
  #mean(as.numeric(as.character(df$Nc)))
  df2$Nc = as.numeric(as.character(df2$Nc))
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bycodonW/","TOP1_",name,"_Nc_TPM.jpg"))
  Nc_cor = cor(df2$Nc,df2[,"TPM"])
  Nc_cor
  plot(df2$Nc,df2[,"TPM"],log = "y",main =paste0(name,"_Nc_TPM_1%  ",Nc_cor),
       xlab="Nc",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================
  write.table(paste(GC_cor,CAI_cor,CBI_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_bycodonW/",name,"_1%_correlation"),
            quote = FALSE,row.names = "correlation", col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM")
  write.table(paste(GC_cor,CAI_cor,CBI_cor,Nc_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_bycodonW/","1%_allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  
}

