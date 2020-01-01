#========================================================================================================
#2019-12-31.Author:Yingying Dong.Taking high expression level genes.
#========================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
dir.create("high_exp")
gtf_array = list.files(getwd(),pattern = "out$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  if(FALSE) # examination
  { 
    gtf = read.table("SRR3589956_abund.out", sep = "\t", header=T,quote = "",fill = T)
    name = "SRR3589956"
  }
  #write.table(gtf$gene_id,file = "id_wait.txt",quote = F,row.names = F, col.names = F)
  gtf = read.table(i, sep = "\t", header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i) 
  name = sub("_abund","",name)
  
  gtf = gtf[,-c(2,3,4,5,6,7)]
  df <- gtf
  df[df == 0] <- NA
  threshhold <- 1
  df = subset(df, df[,3] > threshhold) 
  df = df[complete.cases(df),] # Or df = na.omit(df)
  qRNA = quantile(df$TPM,probs = seq(0,1,0.01))
  qRNA
  hE = df[df$TPM > qRNA[100],]   #Top 2%
  write.table(hE$Gene.ID,file = paste0(name,"_high_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F, col.names = F)
}
