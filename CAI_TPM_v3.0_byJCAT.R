#2019-4-12.Correlation analysis of CAI and global gene expression in all samples.DYY.
species = "C_elegans_Ensl_WBcel235"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/CBI_CAI_byJCAT.txt", sep = ""), sep = "\t", header = F)
dir.create("picture_byJCAT3")
dir.create("correlation_byJCAT3")
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment3/"))
gtf_array = list.files(getwd(),pattern = "SRR\\d.+T$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  gtf = read.table("SRR7160571.gtf.T", sep = "\t", header=F)
  #write.table(gtf$gene_id,file = "id_wait.txt",quote = F,row.names = F, col.names = F)
  gtf = read.table(i, sep = "\t", header = F)
  name = sub("^([^.]*).*", "\\1",i) 
  name = "SRR7160571"
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  names(cbi_cai) = c("transcription_id","CAI")
  cbi_cai$transcription_id = sub(" ","",cbi_cai$transcription_id)
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_byJCAT3/",name,"_CAI_TPM.jpg"))
  CAI_cor = cor(df$CAI,df$TPM)
  CAI_cor
  plot(df$CAI,df$TPM,log = "y",main = paste0(name,"_CAI_TPM  ",CAI_cor),
     xlab="CAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================   
  write.table(CAI_cor,file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_byJCAT3/","CAI_cor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  }

