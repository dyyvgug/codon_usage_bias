#================================================================================================
# 2019-10-07.DYY.Extract GO three classification most significant set,then get their intersection.
# The TPM value of the intersection gene is compared with the TPM value of the gene with high 
# expression & high translation.
#================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
ss = read.table("inCC_inter.txt",header = F,quote = '',sep = '\t')
TPM =  read.table('SRR1804340_riboNum.txt',header = T,quote = '',sep = '\t')
h_e_s = read.table('hE_hT_except_CCinter.txt',header = F,quote = '',sep = '\t')
s_e_s = read.table('in_CC_except_inter.txt',header = F,quote = '',sep = '\t')
names(ss) = "Gene.Name"
ss_TPM = merge(ss,TPM,by = "Gene.Name",all.x = T)
ss_TPM = ss_TPM[complete.cases(ss_TPM),]
names(h_e_s) = "Gene.Name"
h_e_s_TPM = merge(h_e_s,TPM,by = "Gene.Name",all.x = T)
h_e_s_TPM = h_e_s_TPM[complete.cases(h_e_s_TPM),]
names(s_e_s) = "Gene.Name"
s_e_s_TPM = merge(s_e_s,TPM,by = "Gene.Name",all.x = T)
s_e_s_TPM = s_e_s_TPM[complete.cases(s_e_s_TPM),]

png('inter_except_i.png')
p <- boxplot(log(ss_TPM$TPM),log(h_e_s_TPM$TPM),log(s_e_s_TPM$TPM),
             names=c('inter_inter','hEhT_except_inter','inter_except_inter'),col=c("red","yellow",'blue'),ylab = 'log(TPM)')
dev.off()


