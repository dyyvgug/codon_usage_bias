species = "C_elegans_Ensl_WBcel235"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned/ribo_num"))
other_ribo_gene = read.table("40_other_ribo_filter.txt")
hE_hT_ribo_gene = read.table("40_hE_ht__ribo_filter.txt")
other_ribo_gene = cbind(other_ribo_gene,1)
hE_hT_ribo_gene = cbind(hE_hT_ribo_gene,2)
names(other_ribo_gene) = c("ribo_name","RNA_TPM","ribo_TPM","tag")
names(hE_hT_ribo_gene) = c("ribo_name","RNA_TPM","ribo_TPM","tag")
ribo_gene = merge(other_ribo_gene,hE_hT_ribo_gene,all = T)
png("40_compare_ribo_rnaTPM.png")
p <- boxplot(log(other_ribo_gene$RNA_TPM),log(hE_hT_ribo_gene$RNA_TPM),
             names=c('OTHER','hiEXPhiTRA'),col=c("green","red"),ylab = 'log(TPM)')
dev.off()
png("40_compare_ribo_riboTPM.png")
p <- boxplot(log(other_ribo_gene$ribo_TPM),log(hE_hT_ribo_gene$ribo_TPM),
             names=c('OTHER','hiEXPhiTRA'),col=c("green","red"),ylab = 'log(TPM)')
dev.off()
