#2019-6-26.DYY.A cursory look at several ribosomes on each gene.Same species, different laboratory RNAseq and ribo-seq.
species = "C_elegans_Ensl_WBcel235" 
RNAseq_path = "/RNAseq1/experiment2/SRR1056314_abund.out"
RNA = read.table(paste0("/home/hp/Desktop/other_riboseq/",species,RNAseq_path),sep = "\t",header = T,quote = "")
RNA = RNA[,-c(3,4,5,6,7)]
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned"))
dir.create("ribo_num")
ribo_array = list.files(getwd(),pattern = ".out$")
ribo_array
for (i in ribo_array){
  #ribo = read.table("SRR1804340_abund.out",sep = "\t",header = T,quote = "")
  #name = "SRR1804340_abund.out"
  #name = sub("^([^.]*).*", "\\1",name)
  #name = gsub("_abund","",name)
  ribo = read.table(i,sep = "\t",header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i)
  name = sub("_abund","",name)
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
  #=================================================================================================================
  # Extract genes encoding only proteins 
  #=================================================================================================================
  only_protein = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CBI_CAI.txt"),header = T)
  only_protein = only_protein[,-c(3,4,5)]
  only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "transcription_id",all = T)
  only_protein_num = only_protein_num[complete.cases(only_protein_num),]
  only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]
  write.table(only_protein_num,file = paste0("./ribo_num/",name,"_ProtRiboNum.txt"),
              sep = '\t',quote = F,row.names = F)
  #==================================================================================================================
  # Observe the correlation between RNAseq and RIBOseq
  #==================================================================================================================
  co_RNA_ri = cor(only_protein_num$TPM,only_protein_num$ribo_TPM)
  co_RNA_ri
  jpeg(paste0(name,"cor_RNA_ri.jpg"))
  plot(only_protein_num$TPM,only_protein_num$ribo_TPM,main = paste0(name,"cor_RNA_ri  ",co_RNA_ri),
       xlab="RNA_TPM",ylab="ri_TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=205))
  dev.off()
  write.table(co_RNA_ri,file = "cor_RNA_ri.txt",sep = '\t',append = T,quote = FALSE,
              row.names = F, col.names = F)
}
