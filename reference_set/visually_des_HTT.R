#===========================================================================================================
# 2020-10-27.Author:Dong Yingying.Visually describe HTT.
#===========================================================================================================
library(pheatmap)

species = 'Saccharomyces_cerevisiae'
setwd('G:\\研究生阶段\\RIBOseq_RNAseq\\Saccharomyces_cerevisiae\\experiment2')
ribo = read.table('SRR4175354_abund.out',header = T,sep = '\t',quote = "")
rna = read.table('SRR4175342_abund.out',header = T,sep = '\t',quote = "")

ribo2 = data.frame(ribo$Gene.Name,ribo$TPM)
rna2 = data.frame(rna$Gene.Name,rna$TPM)
names(ribo2) = c("name","riboTPM")
names(rna2) = c("name","rnaTPM")

df = merge(ribo2,rna2)
df[df == 0] <-NA
df = df[complete.cases(df),]
rownames(df) = df[,1]
df = df[,-1]
threshhold <- 1
df = subset(df, df[,1] > threshhold) 
df = subset(df, df[,2] > threshhold)
df = df[order(-df[,1],-df[,2]),]

pheatmap(log10(df),cluster_cols = FALSE,cluster_rows = FALSE,
         color = colorRampPalette(c("red", "black", "green"))(50),
         show_rownames=F,filename = "vis_HTT.pdf")

df2 = df[order(-df[,2],-df[,1]),]
pheatmap(log10(df2),cluster_cols = FALSE,cluster_rows = FALSE,
         color = colorRampPalette(c("red", "black", "green"))(50),
         show_rownames=F,filename = "vis_HTT_rna_order.pdf")
