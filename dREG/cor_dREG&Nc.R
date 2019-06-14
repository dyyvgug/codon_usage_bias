setwd("/media/hp/disk2/DYY/dREG/0325out/dREG_peak_analysis/")
dREG_ana = read.table("peak_analysis.txt",sep = "\t",header = T)
Nc_file = read.table("/media/hp/disk1/DYY/reference/annotation/hg38/ref/CBI_Nc_bycodonW.txt",sep = '\t',header = T,quote = "")

df = merge(Nc_file,dREG_ana,by.x = "transcription_id",by.y = "gene",all = T)
df = df [-16]
df[df == 0] <- NA
df = df[complete.cases(df),] # Or df = na.omit(df)
df[df == '*****'] <-NA
df = df[complete.cases(df),]
#sapply(df, class)
#sapply(df$Nc, is.factor)
#mean(as.numeric(as.character(df$Nc)))
df$Nc = as.numeric(as.character(df$Nc))

jpeg(filename = "dREG_score_Nc_cor.jpg")
score_Nc_cor = cor(df$dRGE_score,df$Nc)
score_Nc_cor
plot(df$dRGE_score,df$Nc,log = "x",main = paste0("cor-dRGE_score-Nc  ",score_Nc_cor),
     xlab="dREG_score",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()
jpeg(filename = "dREG_signal_density_Nc_cor.jpg")
dREGdensity_Nc_cor = cor(df$dREG_signal_density,df$Nc)
dREGdensity_Nc_cor
plot(df$dREG_signal_density,df$Nc,log = "x",main = paste0("dREG_signal_density-Nc  ",dREGdensity_Nc_cor),
     xlab="dREG_signal_density",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()
jpeg(filename = "dREG_GEOdensity_Nc_cor.jpg")
GROdensity_cor = cor(df$Gro_density,df$Nc)
GROdensity_cor
plot(df$Gro_density,df$Nc,log = "x",main = paste0("GRO_density-Nc  ",GROdensity_cor),
     xlab="GRO_density",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()

