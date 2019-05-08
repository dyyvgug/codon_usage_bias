library(DOSE)
library(org.Sc.sgd.db)
library(topGO)
library(clusterProfiler)
library(pathview)

keytypes(org.Sc.sgd.db)
species = "Saccharomyces_cerevisiae"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bycodonW2/" ))

data <- read.table("gene_id.txt",header=FALSE)      #单列基因名文件
data$V1 <- as.character(data$V1)                    #需要character格式，然后进行ID转化
                                                    #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
#test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Sc.sgd.db")
test1 = bitr(data$V1, fromType="GENENAME", toType= "ENTREZID", OrgDb="org.Sc.sgd.db")
head(test1,2)
write.table(test1,file = "genename_id.txt",sep = '\t',quote = FALSE,
            row.names = FALSE)

data(geneList, package="DOSE")                      #富集分析的背景基因集
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "GENENAME"), OrgDb = org.Sc.sgd.db)
head(gene.df,2)
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Sc.sgd.db, ont = "CC",level = 3,readable = TRUE)
ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                    universe = names(geneList),     #背景基因集
                    OrgDb = org.Sc.sgd.db,           #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENTREZID',
                    ont = "ALL",                    #也可以是 CC  BP  MF中的一种
                    #ont： 是BP（Biological Process）, CC（Cellular Component）, MF（Molecular Function）
                    pAdjustMethod = "BH",           #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”,
                    #“BH”, “BY”, “fdr”, “none”中的一种
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,               #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 0.2,
                    readable = TRUE)                #Gene ID 转成gene Symbol ，易读
head(ego_ALL,2)
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)
# setReadable
ego_MF <- enrichGO(gene = test1$ENTREZID,ont = "MF",
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, qvalueCutoff = 0.2,
                   minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
ego_MF1 <- setReadable(ego_MF,OrgDb = org.Sc.sgd.db)

dotplot(ego_MF,title="EnrichmentGO_MF_dot")
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")
plotGOgraph(ego_MF)
#=====================================================================================
#KEGG
#=====================================================================================
gene_df <- read.table("genename_id.txt",sep = '\t',header = T,quote = "")
kk <- enrichKEGG(gene = gene_df$GENENAME,
                 organism = 'sce',                 #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
gene_df2 = bitr_kegg(gene_df$ENTREZID,fromType = "ENTREZID",toType = 'PATH',organism='sce')
ego <- enrichKEGG(gene = gene_df$ENTREZID, keyType = "ENTREZID", 
                  organism = 'sce', pvalueCutoff = 0.05, 
                  pAdjustMethod = "BH", qvalueCutoff = 0.05 )


