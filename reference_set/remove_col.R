species = "C_elegans_Ensl_WBcel235"
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned/ribo_num"))
hE_ht_ribo = read.table("40_hE_ht__ribo_filter.txt",header = F,sep = '\t',quote = "")
hE_ht_remove_ribo = read.table("40_hE_ht-ribo.txt",header = T,sep = '\t',quote = "")
hE_ht_ribo_onlyName = hE_ht_ribo[-1]
hE_ht_remove_ribo_onlyName = hE_ht_remove_ribo[-1]
write.table(hE_ht_ribo_onlyName,file = "hE_ht_ribo_onlyName.txt",sep = "\t",quote = F,row.names = F)
write.table(hE_ht_remove_ribo_onlyName,file = "hE_ht_remove_ribo_onlyName.txt",sep = "\t",quote = F,row.names = F)

