#=======================================================================
# Yingying Dong. Calculating weight from RSCU.
#===================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
directory_path <- './weight/'

if (dir.exists(directory_path)) {
  cat("Directory exists.\n")
} else {
  dir.create(directory_path)
}

list.files(getwd())
all_files = list.files(getwd(),pattern = "*")
rscu_array = all_files[!grepl("\\.R$", all_files, ignore.case = TRUE)]
for (i in rscu_array){
  print(i)
  gene_fre = read.table(i,header = T,sep = '\t',quote = "")
  df <- gene_fre
  df$Weights <- ave(df$RSCU,df$AA,FUN=function(x) x/max(x)) # calculate weight
  df = df[,-c(1,3,4,5)]
  df = df[-c(30,61,62,63,64),]
  write.table(df,file = paste0(directory_path,i),sep = '\t',quote = F,row.names = F,col.names = F)
}

