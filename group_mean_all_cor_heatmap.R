#===============================================================================================
# Yingying Dong.2019-12-2.Group mean Nc&TPM,CAI(JCAT)&TPM,CAI(no reference set)&TPM,
#   modified CAI&TPM.
#===============================================================================================
library(dplyr)
library(pheatmap)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_cor = read.csv('Nc_JCat_mCAI.csv',header = T,quote = "")
all_cor = all_cor[-c(44,45,46,47),-1]
names(all_cor) = c("Species","Nc","JCat","no_reference_set","modified")
all_cor[1:43,2:5] <- lapply(all_cor[1:43,2:5],as.numeric)

df <- all_cor %>%
  select(Nc,JCat,no_reference_set,modified) %>%
  group_by('Species') %>%
  summarise(Nc = mean(Nc),JCat = mean(JCat),no_reference_set = mean(no_reference_set),
            modified = mean(modified))

df2 <- aggregate(all_cor[,2:5],list(all_cor$Species),mean)

write.table(df,file = "allcor_Nc_JCat_no_mCAI.csv",sep = ',',row.names = F,quote = F)
write.table(df2,file = "group_Nc_JCat_no_mCAI.txt",sep = '\t',row.names = F,quote = F)

options(digits=3)
rownames(df2) = df2[,1]
df2 = df2[-1]

pheatmap(df2,cluster_cols = FALSE,color = colorRampPalette(c("green", "black", "red"))(50),
         display_numbers = T,show_rownames=T,number_color = "white",clustering_method = "average",
         cellwidth = 27, cellheight = 12, fontsize = 7, filename = "heatmap_Nc_JCat_mCAI.pdf"
         )
