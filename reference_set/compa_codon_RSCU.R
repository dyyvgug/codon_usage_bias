#===========================================================================================================
# 2019-9-19.Modified date:2019-10-28.Author:Dong Yingying.Differences in codon frequency and RSCU value 
#  between DNA sequence codon frequencies of genes that high express high translation and DNA sequences 
#  of genes with low expression and low translation.
#===========================================================================================================
library(pheatmap)
species = "C_elegans_Ensl_WBcel235"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned/ribo_num"))
hE_hT_gene_fre = read.table("hE_hT_codon_fre_RSCU.txt",header = T,sep = '\t',quote = "")
df <- hE_hT_gene_fre
df$Weights <- ave(df$RSCU,df$AA,FUN=function(x) x/max(x)) # calculate weight
df = df[,-c(1,3,4,5)]
df = df[-c(30,61,62,63,64),]
write.table(df,file = "hE_hT_RSCU_weight.txt",sep = '\t',quote = F,row.names = F,col.names = F) # for calculate CAI

lE_lT_gene_fre = read.table("lE_lT_codon_fre_RSCU.txt",header = T,sep = '\t',quote = "")
global_gene = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CDS_codon_fre_RSCU.txt")
                         ,header = T,sep = '\t',quote = "")
hE_hT_rm_ribo = read.table("hE_hT_rm_ribo_codon_fre_RSCU.txt",header = T,sep = '\t',quote = "")
hE_hT_RIBO_only = read.table("hE_hT_ribo_only_codon_fre_RSCU.txt",header = T,sep = '\t',quote = "")
if(FALSE){
  jpeg("40_hE_hT_gene_fre.jpg",width = 1200, height = 600,quality = 100)
  plot(hE_hT_gene_fre$AA,hE_hT_gene_fre$hits,col = "yellow",xlab = "amino acid",ylab = "count")
  dev.off()
  jpeg("40_lE_lT_gene_fre.jpg",width = 1200, height = 600,quality = 100)
  plot(lE_lT_gene_fre$AA,lE_lT_gene_fre$hits,col = "yellow",xlab = "amino acid",ylab = "count")
  dev.off()
}
hE_hT_RSCU = data.frame(hE_hT_gene_fre$AA,hE_hT_gene_fre$codon,hE_hT_gene_fre$RSCU)
names(hE_hT_RSCU) = c("AA","codon","hE_hT_RSCU")
lE_lT_RSCU = data.frame(lE_lT_gene_fre$codon,lE_lT_gene_fre$RSCU)
names(lE_lT_RSCU) = c("codon","lE_lT_RSCU")
global_RSCU = data.frame(global_gene$codon,global_gene$RSCU)
names(global_RSCU) = c("codon","global_RSCU")
hEhT_rm_ribo_RSCU = data.frame(hE_hT_rm_ribo$codon,hE_hT_rm_ribo$RSCU)
names(hEhT_rm_ribo_RSCU) = c("codon","hEhT_rm_ribo_RSCU")
hEhT_RIBO_only_RSCU = data.frame(hE_hT_RIBO_only$codon,hE_hT_RIBO_only$RSCU)
names(hEhT_RIBO_only_RSCU) = c("codon","hEhT_RIBO_only_RSCU")

codon_RSCU = merge(hE_hT_RSCU,lE_lT_RSCU,by = "codon")
codon_RSCU = merge(codon_RSCU,global_RSCU,by = "codon")
codon_RSCU = merge(codon_RSCU,hEhT_rm_ribo_RSCU,by = "codon")
codon_RSCU = merge(codon_RSCU,hEhT_RIBO_only_RSCU,by = "codon")
codon_RSCU = codon_RSCU[order(codon_RSCU$AA,decreasing = F),] # order by AA.
codon_RSCU$codon = gsub(" ", "",codon_RSCU$codon)
codon_RSCU$AA_codon = paste(codon_RSCU$AA,codon_RSCU$codon,sep = "_") # AA_codon
rownames(codon_RSCU) = codon_RSCU[,8] # AA_codon as row name
codon_RSCU = codon_RSCU[,-c(1,2,8)]

pdf(file = "40_codon_RSCU.pdf")
pheatmap(codon_RSCU,cluster_row = FALSE,cluster_cols = FALSE,
         color = colorRampPalette(c("red", "black", "green"))(50))
dev.off() 
