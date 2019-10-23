#2019-9-15.DYY.The ribosome-encoding gene in the gene with high expression high translation 
# and the ribosome encoding gene in other genes except this part of the gene were 
# observed, and the difference of TPM value and TPM value of gene expression between 
# the two groups of ribosomal genes was observed.
species = "mouse_mm10"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment1/aligned/ribo_num"))
other_ribo_gene = read.table("other_ribo_filter.txt")
hE_hT_ribo_gene = read.table("hE_ht__ribo_filter.txt")

names(other_ribo_gene) = c("gene_ID","ribo_name","RNA_TPM","ribo_TPM")
names(hE_hT_ribo_gene) = c("gene_ID","ribo_name","RNA_TPM","ribo_TPM")

png("compare_ribo_rnaTPM.png")
pr1 <- boxplot(other_ribo_gene$RNA_TPM,hE_hT_ribo_gene$RNA_TPM,log ="y",
             names=c('other_ribo','hE_hT_ribo'),col=c("green","red"),ylab = 'RNAseq_TPM')
dev.off()
png("compare_ribo_riboTPM.png")
pr2 <- boxplot(other_ribo_gene$ribo_TPM,hE_hT_ribo_gene$ribo_TPM,log ="y",
             names=c('other_ribo','hE_hT_ribo'),col=c("green","red"),ylab = 'Ribo_seq_TPM')
dev.off()

