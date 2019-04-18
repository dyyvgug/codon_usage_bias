species = "Escherichia_coli"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bycodonW" ))
a = read.table("allcor.txt",sep = '\t',header = F)
b = read.table("1%_allcor.txt",sep = '\t',header = F)
c = read.table("10%_allcor.txt",sep = '\t',header = F)
aveGC = mean(a$V1)
aveGC
aveCAI = mean(a$V2)
aveCAI
aveCBI = mean(a$V3)
aveCBI
aveNc = mean(a$V4)
aveNc
write.table(paste(aveGC,aveCAI,aveCBI,aveNc,sep = '\t'),file = "mean_index.txt",
          quote = FALSE,row.names = "mean_correlation", 
          col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(b$V1),mean(b$V2),mean(b$V3),mean(b$V4),sep = '\t'),file = "mean_1%_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(c$V1),mean(c$V2),mean(c$V3),mean(b$V4),sep = '\t'),file = "mean_10%_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM\tNc_TPM")
