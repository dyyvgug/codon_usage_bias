library(DOSE)
library(org.Sc.sgd.db)
library(topGO)
library(clusterProfiler)
library(pathview)

keytypes(org.Sc.sgd.db)
species = "Saccharomyces_cerevisiae"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bycodonW2/" ))

data = read.table("gene_id.txt",header=FALSE)      #单列基因名文件
data$V1 = as.character(data$V1)                    #需要character格式，然后进行ID转化
                                                   #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
#test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Sc.sgd.db")
test1 = bitr(data$V1, fromType="GENENAME", toType= "ENTREZID", OrgDb="org.Sc.sgd.db")
head(test1,2)
write.table(test1,file = "genename_id.txt",sep = '\t',quote = FALSE,
            row.names = FALSE)
