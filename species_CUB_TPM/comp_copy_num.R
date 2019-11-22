copy_num = read.table("copy_num.txt",header = F,sep = '\t',fill = T)
names(copy_num) = c("AntiCodonsList","tGCN")
copy_num[is.na(copy_num)] <- 0
write.table(copy_num,file = "copy_num.txt",sep = '\t',quote = F,row.names = F)
