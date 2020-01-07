#========================================================================================================
#2019-12-31.Author:Yingying Dong.Taking high and low expression level genes.
#========================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
gtf_array = list.files(getwd(),pattern = "out$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
  if(FALSE) # examination
  { 
    gtf = read.table("SRR7757133_abund.out", sep = "\t", header=T,quote = "",fill = T)
    name = "SRR7757133"
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
  hE = df[df$TPM >= qRNA[100],]   #Top 2%
  write.table(hE$Gene.ID,file = paste0(name,"_high_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F, col.names = F)
  lE = df[df$TPM <= qRNA[50],]    #Low 50%
  write.table(lE$Gene.ID,file = paste0(name,"_low_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  #-----------------------divided into four----------------------------------------
  q2 = quantile(df$TPM,probs = seq(0,1,0.25))
  q2[5]
  df1 = df[df$TPM < q2[2],]      # 0%-25%(not included)
  df2 = df[df$TPM < q2[3],]      # 0%-50%
  df3 = df[df$TPM < q2[4],]      # 0%-75%
  df4 = df[df$TPM < q2[5],]      # 0%-100%
  
  group1 <- df1 
  group2 <- subset(df2, !df2$TPM %in%c(df1$TPM))  # 25%-50%
  group3 <- subset(df3, !df3$TPM %in%c(df2$TPM))  # 50%-75%
  group4 <- subset(df4, !df4$TPM %in%c(df3$TPM))  # 75%-100%
  write.table(group1,file = paste0(name,"_0-25_exp_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  write.table(group2,file = paste0(name,"_25-50_exp_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  write.table(group3,file = paste0(name,"_50-75_exp_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  write.table(group4,file = paste0(name,"_75-100_exp_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  
  #-------------------low expression level explore in detail----------------------------
  lE10 = df[df$TPM <= qRNA[10],]    #Low 10%
  write.table(lE10$Gene.ID,file = paste0(name,"_low10_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  lE20 = df[df$TPM <= qRNA[20],]    #Low 20%
  write.table(lE20$Gene.ID,file = paste0(name,"_low20_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
  lE30 = df[df$TPM <= qRNA[30],]    #Low 10%
  write.table(lE30$Gene.ID,file = paste0(name,"_low30_exp_only_name.txt"),
              sep = '\n',quote = F,row.names = F,col.names = F)
}

