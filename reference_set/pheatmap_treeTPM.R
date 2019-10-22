#============================================================================================
# 2019-10-22.Yingying Dong.According to the sorting of species phylogenetic trees, draw 
#  the correlation heat map of TPM and Nc, and the correlation heat map of TPM and CAI (JCAT), 
#  so as to compare the effect of improved CAI algorithm.
#============================================================================================
library(pheatmap)
setwd("G:\\学习专用\\species_tree")
tree_TPM = read.table("species_TPM.csv",sep = ",",quote = "",header = T)
TPM_diff <- data.frame(tree_TPM$NC_TPM,tree_TPM$CAI_TPM)
names(TPM_diff) = c("Nc_TPM_cor","CAI(JCat)_TPM_cor")
jpeg(filename = "pheatmap_TPM.jpg")
pheatmap(TPM_diff,cluster_row = FALSE,color = colorRampPalette(rainbow(7))(50),
         display_numbers = matrix(ifelse(TPM_diff > 0.2 |TPM_diff < -0.2 , "*", ""), nrow(TPM_diff)),
         show_rownames=F,number_color = "black")
dev.off()
dev.new()         
