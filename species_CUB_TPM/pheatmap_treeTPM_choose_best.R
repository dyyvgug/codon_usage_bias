#============================================================================================
# 2019-12-2.Yingying Dong.According to the sorting of species phylogenetic trees, draw 
#  the correlation heat map of TPM and Nc, and the correlation heat map of TPM and CAI (JCAT), 
#  and the correlation heat map of TPM and modified CAI,so as to compare the effect of 
#  improved CAI algorithm.
#============================================================================================
library(pheatmap)
setwd("G:\\Ñ§Ï°×¨ÓÃ\\species_tree")
tree_TPM = read.csv("species_TPM_only_best.csv",quote = "",header = T)
TPM_diff <- data.frame(tree_TPM$all_Nc_cor,tree_TPM$jcatCAI_TPM)
names(TPM_diff) = c("Nc_TPM_cor","CAI(JCat)_TPM_cor")

pheatmap(TPM_diff,cluster_row = FALSE,color = colorRampPalette(c("green", "black", "red"))(50),
         display_numbers = TRUE,show_rownames=F,number_color = "white",
         cellwidth = 27, cellheight = 12, fontsize = 7,
         filename = "heatmap_species_Nc_Jcat.pdf")

        
