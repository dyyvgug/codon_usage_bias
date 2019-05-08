setwd("/media/hp/disk2/DYY/dREG/0325out/dREG_peak_analysis/")
dREG_ana = read.table("peak_analysis.txt",sep = "\t",header = T)
CAI_file = read.table("/media/hp/disk1/DYY/reference/annotation/hg38/ref/CAI_byJcat.txt")
names(CAI_file) = c("gene_name","CAI")
df = merge(dREG_ana,CAI_file,by.x = "gene",by.y = "gene_name",all = T)
df[df == 0] <- NA
df = df[complete.cases(df),] # Or df = na.omit(df)
jpeg(filename = "dREG_score_cor.jpg")
score_cor = cor(df$dRGE_score,df$CAI)
score_cor
plot(df$dRGE_score,df$CAI,log = "x",main = paste0("cor-dRGE_score-CAI  ",score_cor),
     xlab="dREG_score",ylab="CAI",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()
jpeg(filename = "dREG_signal_density_cor.jpg")
dREGdensity_cor = cor(df$dREG_signal_density,df$CAI)
dREGdensity_cor
plot(df$dREG_signal_density,df$CAI,log = "x",main = paste0("dREG_signal_density-CAI  ",dREGdensity_cor),
     xlab="dREG_signal_density",ylab="CAI",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()
GROdensity_cor = cor(df$Gro_density,df$CAI)
GROdensity_cor
plot(df$Gro_density,df$CAI,log = "x",main = paste0("GRO_density-CAI  ",GROdensity_cor),
     xlab="GRO_density",ylab="CAI",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
dev.off()

