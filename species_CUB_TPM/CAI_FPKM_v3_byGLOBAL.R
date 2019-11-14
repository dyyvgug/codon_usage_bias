#========================================================================================================
# 2019-11-5.Author:Yingying Dong.Correlation analysis of CAI by no reference set (weight from global genome) 
#  and global gene expression in all samples.
#========================================================================================================

species = "Komagataella_pastoris"
speA_exp = "Kp2"

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                        species, "/ref/CBI_CAI2.txt", sep = ""), sep = "\t", header = T,quote = '')
cai = data.frame(cai$transcription_id,cai$CAI)
names(cai) = c("gene_id","CAI")
dir.create(paste0("picture_byGLOBAL"))
dir.create(paste0("correlation_byGLOBAL"))
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/"))
gtf_array = list.files(getwd(),pattern = ".+csv$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  if(FALSE) # examination
  { 
    fpkm = read.csv("W_D_1.csv", header=T,quote = "")
    name = "W_D_1"
  }
  #write.table(gtf$gene_id,file = "id_wait.txt",quote = F,row.names = F, col.names = F)
  fpkm = read.csv(i,header = T,fill = T)
  name = sub("^([^.]*).*", "\\1",i) 
  
  cai_FPKM = merge(cai, fpkm, by.x = "gene_id", by.y = "Locus.ID", all = T)
  cai_FPKM[cai_FPKM == 0] <- NA
  cai_FPKM = cai_FPKM[complete.cases(cai_FPKM),]
  df <- cai_FPKM
  #===========================================================================
  # Correlation_CAI_FPKM 
  #===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,
                   "/picture_byGLOBAL/",name,"_CAI_FPKM.jpg"))
  CAI_cor = cor(df$CAI,df$FPKM)
  CAI_cor
  plot(df$CAI,df$FPKM,log = "y",main = paste0(name,"_CAI_FPKM  ",CAI_cor),
       xlab="CAI",ylab="FPKM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
  #===============================================================================
  # Write out data,convenient to calculate the average
  #===============================================================================   
  write.table(CAI_cor,file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_byGLOBAL/","CAI_cor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_byGLOBAL/" ))
a = read.table("CAI_cor.txt",sep = '\t',header = F)

aveCAI = mean(a$V1)
aveCAI

write.table(aveCAI,file = "mean_CAI.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_FPKM")