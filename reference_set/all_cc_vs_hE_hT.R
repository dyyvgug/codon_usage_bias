#================================================================================================
# 2019-10-07.Modified date:2020-1-14.DYY.Extract GO three classification most significant set,
# then get their intersection,except the intersection between MF & BP.
# The TPM value of the intersection gene is compared with the TPM value of the gene with high 
# expression & high translation.
#================================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
s = read.table("3_completion_set.txt",header = F,quote = '',sep = '\t')
TPM =  read.table('SRR1804340_riboNum.txt',header = T,quote = '',sep = '\t')
hE_hT = read.table('SRR1804340_hE_ht_gene.txt',header = T,quote = '',sep = '\t')
names(s) = "Gene.Name"
s_TPM = merge(s,TPM,by = "Gene.Name",all.x = T)
s_TPM = s_TPM[complete.cases(s_TPM),]
s_TPM = unique(s_TPM)

pdf('40_three_set_vs_hE_hT_TPMbox.pdf')
p1 <- boxplot(s_TPM$TPM,hE_hT$TPM,TPM$TPM,log = 'y',boxwex = .5,cex = .5,cex.lab = 1.4,pch = 18,
             names=c('three intersection','hEhT','whole genome'),col=c("red","yellow","lightblue"),
             ylab = 'RNA-seq (TPM)',main = '(n = 84)\t\t(n=192)\t\t(n=12444)')

p2 <- boxplot(s_TPM$ribo_TPM,hE_hT$ribo_TPM,TPM$TPM,log = 'y',boxwex = .5,cex = .5,cex.lab = 1.4,
             names=c('three intersection','hEhT','whole genome'),pch = 18,
             col=c("red","yellow","lightblue"),ylab = 'Ribo-seq (TPM)',main = '(n = 84)\t\t(n=192)\t\t(n=12444)')
dev.off()
#=======================================================================================================
# Extract GO three classification most significant set,then get their intersection.
# The TPM value of the intersection gene is compared with the TPM value of the gene with high 
# expression & high translation.
#=======================================================================================================
ss = read.table("inter_inter.txt",header = F,quote = '',sep = '\t')
h_e_s = read.table('hE_hT_except_inter.txt',header = F,quote = '',sep = '\t')
s_e_s = read.table('inter_except_inter.txt',header = F,quote = '',sep = '\t')
names(ss) = "Gene.Name"
ss_TPM = merge(ss,TPM,by = "Gene.Name",all.x = T)
ss_TPM = ss_TPM[complete.cases(ss_TPM),]
names(h_e_s) = "Gene.Name"
h_e_s_TPM = merge(h_e_s,TPM,by = "Gene.Name",all.x = T)
h_e_s_TPM = h_e_s_TPM[complete.cases(h_e_s_TPM),]
names(s_e_s) = "Gene.Name"
s_e_s_TPM = merge(s_e_s,TPM,by = "Gene.Name",all.x = T)
s_e_s_TPM = s_e_s_TPM[complete.cases(s_e_s_TPM),]

pdf('inter_except_inter_TPMbox.pdf')
p3 <- boxplot(ss_TPM$TPM,h_e_s_TPM$TPM,s_e_s_TPM$TPM,TPM$TPM,log = 'y',
             names=c('group I','group II','group III','whole genome'),
             col=c("orange","yellow",'red','lightblue'),ylab = 'RNA-seq (TPM)',
             boxwex = .5,cex = .5,cex.lab = 1.4,pch = 18,
             main = '(n = 69) (n = 123) (n = 15) (n = 12444)')

p4 <- boxplot(ss_TPM$ribo_TPM,h_e_s_TPM$ribo_TPM,s_e_s_TPM$ribo_TPM,TPM$ribo_TPM,log = 'y',
             names=c('group I','group II','group III','whole genome'),
             col=c("orange","yellow",'red','lightblue'),ylab = 'Ribo-seq (TPM)',
             boxwex = .5,cex = .5,cex.lab = 1.4,pch = 18,
             main = '(n = 69) (n = 123) (n = 15) (n = 12444)')
dev.off()


