#===========================================================================================================
# 2019-9-19.Differences in codon frequency between DNA sequence codon frequencies 
# of genes that express high translation and DNA sequences of genes with low expression and low translation.
#===========================================================================================================
species = "C_elegans_Ensl_WBcel235"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned/"))
hE_hT_gene_fre = read.table("40_hE_hT_codon_fre.txt",header = T,sep = '\t',quote = "")
lE_lT_gene_fre = read.table("40_lE_lT_codon_fre.txt",header = T,sep = '\t',quote = "")
jpeg("40_hE_hT_gene_fre.jpg",width = 1200, height = 600,quality = 100)
plot(hE_hT_gene_fre$AA,hE_hT_gene_fre$hits,col = "yellow",xlab = "amino acid",ylab = "count")
dev.off()
jpeg("40_lE_lT_gene_fre.jpg",width = 1200, height = 600,quality = 100)
plot(lE_lT_gene_fre$AA,lE_lT_gene_fre$hits,col = "yellow",xlab = "amino acid",ylab = "count")
dev.off()
library(pheatmap)

codon_fre_diff = read.table("40_codon_fre_dif.csv",header = T,sep = ',',quote = "")
codon_fre_diff = codon_fre_diff[,-c(1,3,5)]
rownames(codon_fre_diff) = codon_fre_diff[,1]
codon_fre_diff = codon_fre_diff[-1]
colnames(codon_fre_diff) = c("h_codon_fre","l_codon_fre")
dev.new()
jpeg("40_codon_fre_dif.jpg",width = 600, height = 1200, quality = 100)
pheatmap(codon_fre_diff,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off() 
