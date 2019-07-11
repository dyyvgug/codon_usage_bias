#====================================================================================== 
# 2019-7-11.Author:Dong Yingying.Roughly  observe the translation efficiency.
# Translation efficiency is subtracted from the TPM value quantified by a sample of 
# RNAseq transcripts and the corresponding TPM value of the RIBO-seq transcript(same lab)
#======================================================================================
  species = "hg19" 
  SRRnum = "SRP173217"
  SRAnum = "SRR8309258_abund.out"
  setwd(paste0("~/Desktop/ribosome_profiling/",SRRnum,"/aligned/"))
  dir.create("ribo_num")
  RNA = read.table(SRAnum,sep = "\t",header = T,quote = "")
  RNA = RNA[,-c(3,4,5,6,7)]
  RIBOnum = "SRR8309248_abund.out"
  ribo = read.table(paste0("../aligned_ri/",RIBOnum),sep = "\t",header = T,quote = "")
  name = RIBOnum
  name = sub("^([^.]*).*", "\\1",name)
  name = gsub("_abund","",name)
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
  write.table(RNA_ribo_num,file = paste0("./ribo_num/",name,"_riboNum.txt"),
              sep = "\t",quote = F,row.names = F)
  #**************************Extract genes encoding only proteins********************************  
  only_protein = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CBI_CAI.txt"),header = T)
  only_protein = only_protein[,-c(3,4,5)]
  only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "transcription_id",all = T)
  only_protein_num = only_protein_num[complete.cases(only_protein_num),]
  only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]
  write.table(only_protein_num,file = paste0("./ribo_num/",name,"_ProtRiboNum.txt"),
              sep = '\t',quote = F,row.names = F)
  
