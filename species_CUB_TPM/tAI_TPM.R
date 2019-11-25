#========================================================================================================
# 2019-11-22.Author:Yingying Dong. Correlation analysis of tAI by stAIcalc and 
# global gene expression in all samples.
#========================================================================================================
library(getopt)
command=matrix(c("species","s",1,"character",
                 "experiment","e",1,"numeric",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$species) || is.null(args$experiment)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}
species = args$species
exp = args$experiment

species = "Escherichia_coli"
exp = "2"

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/tAI/"))
tai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                        species, "/tAI/exportstAiGenesBackup.csv", sep = ""), sep = ",", fill = T,header = T,quote = '')
tai$X.Gene.Name. = sub("\"","",tai$X.Gene.Name.)
write.table(tai,file = "tAI.txt",sep = ',',quote = F,row.names = F,col.names = F)

tAI = read.table("tAI.txt",header = F,sep = '\t')
dir.create(paste0("picture_bytAI",exp))
dir.create(paste0("correlation_bytAI",exp))
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment",exp,"/"))
gtf_array = list.files(getwd(),pattern = "[SE]RR\\d.+out$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  if(FALSE) # examination
  { 
    gtf = read.table("SRR6111094_abund.out", sep = "\t", header=T,quote = "",fill = T)
    name = "SRR6111094"
  }
  #write.table(gtf$gene_id,file = "id_wait.txt",quote = F,row.names = F, col.names = F)
  gtf = read.table(i, sep = "\t", header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i) 
  name = sub("_abund","",name)
  
  gtf = gtf[,-c(2,3,4,5,6,7)]
  names(tAI) = c("gene_id","tAI")
  df = merge(tAI, gtf, by.x = "gene_id",by.y = "Gene.ID", all = T)
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)
  #===========================================================================
  # Correlation_tAI_TPM 
  #===========================================================================
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,
                   "/picture_bytAI",exp,"/",name,"_tAI_TPM.jpg"))
  tAI_cor = cor(df$tAI,df$TPM)
  tAI_cor
  tAI_p = cor.test(df$tAI,df$TPM)
  tAI_p
  plot(df$tAI,df$TPM,log = "y",main = paste0(name,"_tAI_TPM  ",tAI_cor),
       xlab="tAI",ylab="TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=255))
  dev.off()
  #===============================================================================
  # Write out data,convenient to calculate the average
  #===============================================================================   
  write.table(tAI_cor,file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_bytAI",exp,"/","tAI_cor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bytAI",exp,"/" ))
a = read.table("tAI_cor.txt",sep = '\t',header = F)

avetAI = mean(a$V1)
avetAI

write.table(avetAI,file = "mean_tAI.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\ttAI_TPM")