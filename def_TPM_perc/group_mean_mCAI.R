#===============================================================================================
# Yingying Dong.2019-11-12.Modified date:2019-12-2.Group mean def mCAI and TPM and pheatmap.
#===============================================================================================
library(ggplot2)
library(dplyr)
library(pheatmap)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
def_cor = read.csv('def_mCAI_TPM.csv',header = T)
names(def_cor) = c('ID','Species','all','top2','top10','low48','def48_99')
def_cor[1:41,3:7] <- lapply(def_cor[1:41,3:7],as.numeric) 

df <- def_cor %>%
  select(all,top2,top10,low48,def48_99) %>%
  group_by('Species') %>%
  summarise(all = mean(all), top2 = mean(top2), top10 = mean(top10),
            low48 = mean(low48),def48_99= mean(def48_99))
df2 <- aggregate(def_cor[,3:7], list(def_cor$Species), mean)
df2$Group.1 <- as.factor(df2$Group.1)
names(df2) <- c('species','all','top2','top10','low48','def48_99')
df3 <- as.data.frame(t(df2))   # converting rows into columns
names(df3) <- as.matrix(df3[1, ]) # take the first row as the column names
df3 <- df3[-1, ]
df3[] <- lapply(df3, function(x) type.convert(as.character(x))) # converting the columns to their appropriate types
write.table(df,file = "all_ave_mCAItpm_cor.txt",sep = '\t',row.names = F,quote = F)
write.table(df2,file = "ave_mCAItpm_cor.txt",sep = '\t',row.names = F,quote = F)
write.table(df3,file = 'ave_mCAItpm_cor_invert.txt',sep = '\t',row.names = F,quote = F)

options(digits=3)
rownames(df2) = df2[,1]
df2 = df2[-1]
pheatmap(df2,cluster_cols = FALSE,color = colorRampPalette(c("white", "blue"))(50),
         display_numbers = T,show_rownames=T,number_color = "black",clustering_method = "average",
         cellwidth = 27, cellheight = 12, fontsize = 7,border=FALSE,
         filename = "heatmap_def_mCAI.pdf"
)
