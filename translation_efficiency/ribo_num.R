#2019-6-26.DYY.A cursory look at several ribosomes on each gene.
RNA = read.table("~/Desktop/SRP136094/SRR6869737_abund.out",sep = "\t",header = T,quote = "")
RNA = RNA[,-c(3,4,5,6,7)]
setwd("~/Desktop/SRP136094/fq")
ribo_array = list.files(getwd(),pattern = ".out$")
ribo_array
for (i in ribo_array){
  #ribo = read.table("~/Desktop/SRP136094/fq/SRR6869751_25_abund.out",sep = "\t",header = T,quote = "")
  #name = "SRR6869751_25_abund.out"
  #name = sub("^([^.]*).*", "\\1",name)
  #name = gsub("_25_abund","",name)
  ribo = read.table(i,sep = "\t",header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i)
  name = gsub("_25_abund","",i)
  ribo = ribo[,-c(2,3,4,5,6,7)]
  names(ribo) = c("Gene.ID","ribo_FPKM","ribo_TPM")
  RNA_ribo = merge(RNA,ribo,by = "Gene.ID",all = T)
  RNA_ribo[RNA_ribo == 0] <-NA
  RNA_ribo = RNA_ribo[complete.cases(RNA_ribo),]
  #RNA_ribo_num = cbind(RNA_ribo,round(RNA_ribo$TPM / RNA_ribo$ribo_TPM,3))
  ribo_num = round(RNA_ribo$ribo_TPM / RNA_ribo$TPM,3)
  RNA_ribo_num = cbind(RNA_ribo,ribo_num)
  RNA_ribo_num = RNA_ribo_num[order(RNA_ribo_num$ribo_num,decreasing = T),]
  RNA_ribo_num$Gene.ID = gsub("_.*", "", RNA_ribo_num[,1])
  write.table(RNA_ribo_num,file = paste0("~/Desktop/SRP136094/fq/",name,"_riboNum.txt"),
              sep = "\t",quote = F,row.names = F)
#**************************Extract genes encoding only proteins********************************  
  only_protein = read.table("/media/hp/disk1/DYY/reference/annotation/hg19/ref/only_protein.txt",header = T)
  only_protein = only_protein[,-c(3,4,5)]
  only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "transcription_id",all = T)
  only_protein_num = only_protein_num[complete.cases(only_protein_num),]
  only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]
  write.table(only_protein_num,file = paste0("~/Desktop/SRP136094/fq/",name,"_ProtRiboNum.txt"),
              sep = '\t',quote = F,row.names = F)
}

