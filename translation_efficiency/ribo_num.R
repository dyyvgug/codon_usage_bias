#====================================================================================== 
# 2019-7-11.Author:Dong Yingying.Roughly  observe the translation efficiency.
# Translation efficiency is subtracted from the TPM value quantified by a sample of 
# RNAseq transcripts and the corresponding TPM value of the RIBO-seq transcript(same lab)
#======================================================================================
species = "Drosophila_melanogaster" 
SRAnum = "SRR5667266_abund.out"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment1/aligned/"))
dir.create("ribo_num")
RNA = read.table(SRAnum,sep = "\t",header = T,quote = "")
RNA = RNA[,-c(3,4,5,6,7)]
RIBOnum = "SRR5667267_abund.out"
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
#RNA_ribo_num$Gene.ID = gsub("_.*", "", RNA_ribo_num[,1])
write.table(RNA_ribo_num,file = paste0("./ribo_num/",name,"_riboNum.txt"),
              sep = "\t",quote = F,row.names = F)
#=======================================================================================================
#  Extract genes encoding only proteins
#=======================================================================================================  
only_protein = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CBI_CAI.txt"),sep = ",",header = T, quote = "")
only_protein = only_protein[,-c(3,4,5)]
only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "Gene_Symbol",all = T)
only_protein_num = only_protein_num[complete.cases(only_protein_num),]
only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]
write.table(only_protein_num,file = paste0("./ribo_num/",name,"_ProtRiboNum.txt"),
              sep = '\t',quote = F,row.names = F)
#=======================================================================================================
# Observe the correlation between RNAseq and RIBOseq
#=======================================================================================================
co_RNA_ri = cor(only_protein_num$TPM,only_protein_num$ribo_TPM)
co_RNA_ri
jpeg(paste0(name,"cor_RNA_ri.jpg"))
plot(only_protein_num$TPM,only_protein_num$ribo_TPM,log ="xy",main = paste0(name,"cor_RNA_ri  ",co_RNA_ri),
       xlab="RNA_TPM",ylab="ri_TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=205))
dev.off()
write.table(co_RNA_ri,file = "cor_RNA_ri.txt",sep = '\t',append = T,quote = FALSE,
              row.names = F, col.names = F)

  
