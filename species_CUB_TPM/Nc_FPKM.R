#2019-6-6.Correlation analysis of Nc and various gene expression levels in all samples.DYY.
species = "Komagataella_pastoris"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_bycodonW.txt", sep = ""), 
                      sep = "\t", header = T,quote = "")
dir.create("picture_Nc1")
dir.create("correlation_Nc1")
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
gtf_array = list.files(getwd(),pattern = ".+csv$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  fpkm = read.csv(i,header = T,fill = T)
  name = sub("^([^.]*).*", "\\1",i) 
  df = merge(cbi_cai, fpkm, by.x = "transcription_id", by.y = "Locus.ID", all = T)
  df = df [-16]
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
  df[df == '*****'] <-NA
  df = df[complete.cases(df),]
  #sapply(df, class)
  #sapply(df$Nc, is.factor)
  #mean(as.numeric(as.character(df$Nc)))
  df$Nc = as.numeric(as.character(df$Nc))
  q = quantile(df$FPKM[df$FPKM > 0], prob = seq(0,1,0.01))
  q
  df2 = df[df$FPKM >= q[99],] # top 2%
  df3 = df[df$`FPKM` >= q[90],] # top 10%
  df_l = df[df$FPKM <= q[48],] # low 48%
  df4 = df[df$FPKM >= q[100],] # top 1%
  df5 = df[df$FPKM >= q[48],] # top 10%
  df_def = subset(df5, !df5$FPKM %in%c(df4$FPKM)) # Remove the a data frame from the b data frame.
#===========================================================================
# For allGene correlation Nc and FPKM 
#===========================================================================
  all_Nc_cor = cor(df$Nc,df$FPKM)
  all_Nc_cor
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/",name,"all_Nc_FPKM.jpg"))
  plot(df$Nc,df$FPKM,log = "y",main =paste0(name,"_all_Nc_FPKM  ",all_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For Top2%Gene correlation Nc and FPKM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","TOP1_",name,"_Nc_FPKM.jpg"))
  top_Nc_cor = cor(df2$Nc,df2$FPKM)
  top_Nc_cor
  plot(df2$Nc,df2$FPKM,log = "y",main =paste0(name,"_Nc_FPKM_1%  ",top_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For Top10%Gene correlation Nc and FPKM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","h_",name,"_Nc_FPKM.jpg"))
  high_Nc_cor = cor(df3$Nc,df3[,"FPKM"])
  high_Nc_cor
  plot(df3$Nc,df3[,"FPKM"],log = "y",main =paste0(name,"_Nc_FPKM_h  ",high_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#==========================================================================
# For 0-48%Gene correlation Nc and FPKM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","l48_",name,"_Nc_FPKM.jpg"))
  low_Nc_cor = cor(df_l$Nc,df_l[,"FPKM"])
  low_Nc_cor
  plot(df_l$Nc,df_l[,"FPKM"],log = "y",main =paste0(name,"_Nc_FPKM_l48  ",low_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===========================================================================
# For defGene correlation Nc and FPKM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_Nc1/","def_",name,"_Nc_FPKM.jpg"))
  def_Nc_cor = cor(df_def$Nc,df_def$FPKM)
  def_Nc_cor
  plot(df_def$Nc,df_def$FPKM,log = "y",main =paste0(name,"_Nc_FPKM_def  ",def_Nc_cor),
       xlab="Nc",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================   
  write.table(paste(all_Nc_cor,top_Nc_cor,high_Nc_cor,low_Nc_cor,def_Nc_cor,sep = '\t'),file = paste0
              ("/media/hp/disk1/DYY/reference/annotation/", 
                species,"/correlation_Nc1/",name,"correlation_bycodonW"),
            quote = FALSE,row.names = "correlation", 
            col.names = "\tGC_FPKM\tCAI_FPKM\tCBI_FPKM\tNc_FPKM")
  write.table(paste(all_Nc_cor,top_Nc_cor,high_Nc_cor,low_Nc_cor,def_Nc_cor,sep = '\t'),file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_Nc1/","allcor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  }
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_Nc1/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
names(a) = c("all_Nc_cor","top_Nc_cor","high_Nc_cor","low_Nc_cor","def_Nc_cor")
write.table(paste(mean(a$all_Nc_cor),mean(a$top_Nc_cor),mean(a$high_Nc_cor),
                  mean(a$low_Nc_cor),mean(a$def_Nc_cor),sep = '\t'),file = "mean_cor.txt",
            quote = FALSE,row.names = "mean_Nc_correlation", 
            col.names = "all_Nc_cor\ttop_Nc_cor\thigh_Nc_cor\tlow_Nc_cor\tdef_Nc_cor")
