species = "Escherichia_coli"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bycodonW2/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
b = read.table("1%_allcor.txt",sep = '\t',header = F)
c = read.table("10%_allcor.txt",sep = '\t',header = F)
d = read.table("def_allcor.txt",sep = '\t',header = F)
e = read.table("l48_allcor.txt",sep = '\t',header = F)

aveGC = mean(a$V1)
aveGC
aveCAI = mean(a$V2)
aveCAI
aveCBI = mean(a$V3)
aveCBI
aveNc = mean(a$V4)
aveNc
write.table(paste(aveCAI,aveCBI,aveNc,sep = '\t'),file = "mean_index.txt",
          quote = FALSE,row.names = "mean_correlation", 
          col.names = "\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(b$V2),mean(b$V3),mean(b$V4),sep = '\t'),file = "mean_1%_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(c$V2),mean(c$V3),mean(c$V4),sep = '\t'),file = "mean_10%_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(d$V1),mean(d$V2),mean(d$V3),sep = '\t'),file = "mean_def_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM\tCBI_TPM\tNc_TPM")
write.table(paste(mean(e$V1),mean(e$V2),mean(e$V3),sep = '\t'),file = "mean_l48_index.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM\tCBI_TPM\tNc_TPM")
