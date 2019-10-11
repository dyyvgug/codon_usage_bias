#================================================================================================
# 2019-10-07.DYY.Extract GO three classification most significant set,then get their intersection,
# except the intersection between MF & BP.
# The TPM value of the intersection gene is compared with the TPM value of the gene with high 
# expression & high translation.
#================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
s = read.table("KEGG_only_symbol.txt",header = F,quote = '',sep = '\t')
TPM =  read.table('SRR1804340_riboNum.txt',header = T,quote = '',sep = '\t')
hE_hT = read.table('SRR1804340_hE_ht_gene.txt',header = T,quote = '',sep = '\t')
names(s) = "Gene.Name"
s_TPM = merge(s,TPM,by = "Gene.Name",all.x = T)
s_TPM = s_TPM[complete.cases(s_TPM),]

png('oxi_pho_vs_hE_hT_TPMbox.png')
p <- boxplot(s_TPM$TPM,hE_hT$TPM,TPM$TPM,log = 'y',
             names=c('oxi_pho','hE_hT','global_gene'),col=c("red","yellow","blue"),ylab = 'TPM')
dev.off()
png('oxi_pho_vs_hE_hT_riboTPMbox.png')
p <- boxplot(s_TPM$ribo_TPM,hE_hT$ribo_TPM,TPM$ribo_TPM,log = 'y',
             names=c('oxi_pho','hE_hT','global_gene'),col=c("red","yellow","blue"),ylab = 'ribo_TPM')
dev.off()
