setwd("/media/hp/disk2/DYY2/dREG/0325out/dREG_peak_analysis/")
dREG_ana = read.table("peak_analysis.txt",sep = "\t",header = T)
Nc_file = read.table("/media/hp/disk1/DYY/reference/annotation/hg38/ref/CBI_CAI_bycodonW.txt",sep = '\t',header = T,quote = "")

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

pdf("dREG_score_Nc_cor.pdf")
score_Nc_cor = cor(df$dRGE_score,df$Nc)
score_Nc_cor
score_Nc_p = cor.test(df$dRGE_score,df$Nc)

plot(df$dRGE_score,df$Nc,log = "x",main = paste0("cor-dRGE_score-Nc  ",round(score_Nc_cor,5),"  p =",round(score_Nc_p$p.value,5)),
     xlab="dREG_score",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()

pdf("dREG_signal_density_Nc_cor.jpg")
dREGdensity_Nc_cor = cor(df$dREG_signal_density,df$Nc)
dREGdensity_Nc_cor
dREGdensity_Nc_p = cor.test(df$dREG_signal_density,df$Nc)
plot(df$dREG_signal_density,df$Nc,log = "x",main = paste0("dREG_signal_density-Nc  ",round(dREGdensity_Nc_cor,5),"  p=",round(dREGdensity_Nc_p$p.value,5)),
     xlab="dREG_signal_density",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()

pdf("dREG_GEOdensity_Nc_cor.jpg")
GROdensity_cor = cor(df$Gro_density,df$Nc)
GROdensity_cor
GROdensity_cor_p = cor.test(df$Gro_density,df$Nc)
plot(df$Gro_density,df$Nc,log = "x",main = paste0("GRO_density-Nc  ",round(GROdensity_cor,5),round(GROdensity_cor_p$p.value,5)),
     xlab="GRO_density",ylab="Nc",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()

