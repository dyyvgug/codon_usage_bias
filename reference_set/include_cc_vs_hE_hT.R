#================================================================================================
# 2019-10-07.DYY.Extract GO three classification most significant set,then get their intersection,
# except the intersection between MF & BP.
# The TPM value of the intersection gene is compared with the TPM value of the gene with high 
# expression & high translation.
#================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
s = read.table("include_cc.txt",header = F,quote = '',sep = '\t')
TPM =  read.table('SRR1804340_riboNum.txt',header = T,quote = '',sep = '\t')
hE_hT = read.table('SRR1804340_hE_ht_gene.txt',header = T,quote = '',sep = '\t')
names(s) = "Gene.Name"
s_TPM = merge(s,TPM,by = "Gene.Name",all.x = T)
s_TPM = s_TPM[complete.cases(s_TPM),]

png('include_cc_vs_hE_hT_TPMbox.png')
p <- boxplot(log(s_TPM$TPM),log(hE_hT$TPM),
             names=c('include_cc','hE_hT'),col=c("green","red"),ylab = 'log(TPM)')
dev.off()
png('include_cc_vs_hE_hT_riboTPMbox.png')
p <- boxplot(log(s_TPM$ribo_TPM),log(hE_hT$ribo_TPM),
             names=c('include_cc','hE_hT'),col=c("green","red"),ylab = 'log(ribo_TPM)')
dev.off()
