#===============================================================================================
# Author:Yingying Dong.2019-12-3.Choose best Nc&TPM,CAI(JCAT)&TPM,CAI(no reference set)&TPM,
#   modified CAI&TPM.Drawing heatmap.
#===============================================================================================
library(dplyr)
library(pheatmap)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_cor = read.csv('Nc_JCat_mCAI_choose.csv',header = T,quote = "")
names(all_cor) = c("Species","Nc","JCat","no_reference_set","modified")
all_cor[,2:5] <- lapply(all_cor[,2:5],as.numeric)

df <- all_cor %>%
  select(Nc,JCat,no_reference_set,modified) %>%
  group_by('Species') %>%
  summarise(Nc = mean(Nc),JCat = mean(JCat),no_reference_set = mean(no_reference_set),
            modified = mean(modified))

write.table(df,file = "allcor_Nc_JCat_no_mCAI_choose.csv",sep = ',',row.names = F,quote = F)

options(digits=3)
rownames(all_cor) = all_cor[,1]
all_cor = all_cor[-1]

pheatmap(all_cor,cluster_cols = FALSE,color = colorRampPalette(c("green", "black", "red"))(50),
         display_numbers = T,show_rownames=T,number_color = "white",clustering_method = "average",
         cellwidth = 27, cellheight = 12, fontsize = 7, cluster_rows = FALSE,
         filename = "heatmap_Nc_JCat_mCAI_choose.pdf"
         )
