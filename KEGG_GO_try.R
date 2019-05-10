# GO and KEGG analysis.Dong Yingying.2019-5-9.
#library(DOSE)
library(org.EcK12.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

keytypes(org.EcK12.eg.db)
species = "Escherichia_coli"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bycodonW3/" ))
#==============================================================================================
# Read data and conversion id
#==============================================================================================
symbol_id = read.table("gene_id.txt",header=FALSE)      
symbol_id = as.character(symbol_id$V1)                    
entrez_id = bitr(symbol_id, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.EcK12.eg.db")
head(entrez_id,2)
write.table(entrez_id,file = "symbol_entrez_id.txt",sep = '\t',quote = FALSE,
            row.names = FALSE)
only_entrezID = entrez_id[-1]
write.table(only_entrezID,file = "ENTREZID.txt",sep = '\t',quote = FALSE,
            row.names = FALSE)
#==============================================================================================
# GO analysis ORA(over-representation analysis)
#==============================================================================================
#**************all(Biological Process,Cellular Component,Molecular Function)*******************
ogo_all = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = "org.EcK12.eg.db",        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
head(ogo_all,2)
write.table(ogo_all,"GO_aLL_enrich.txt",row.names =FALSE)
dotplot(ogo_all,title="EnrichmentGO_all_dot")
barplot(ogo_all,showCategory=10,title="EnrichmentGO_all_bar")
#*********************************** MF(Molecular Function)*************************************
ogo_MF = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = "org.EcK12.eg.db",        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
head(ogo_MF)
write.table(ogo_MF,"GO_MF_enrich.txt",row.names =FALSE)
dotplot(ogo_MF,title="EnrichmentGO_MF_dot")
barplot(ogo_MF,showCategory=10,title="EnrichmentGO_MF_bar")
plotGOgraph(ogo_MF)
#.rs.restartR()                    # if occur error 
goplot(ogo_MF)
emapplot(ogo_MF,showCategory = 30)
cnetplot(ogo_MF,showCategory = 5)
#********************************BP(Biological Process)*******************************
ogo_BP = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = "org.EcK12.eg.db",        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
head(ogo_BP)
write.table(ogo_BP,"GO_BP_enrich.txt",row.names =FALSE)
dotplot(ogo_BP,title="EnrichmentGO_BP_dot")
barplot(ogo_BP,showCategory=10,title="EnrichmentGO_BP_bar")
plotGOgraph(ogo_BP)
#.rs.restartR()                    # if occur error 
goplot(ogo_BP)
emapplot(ogo_BP,showCategory = 30)
cnetplot(ogo_BP,showCategory = 5)
#********************************CC(Cellular Component)*******************************
ogo_CC = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = "org.EcK12.eg.db",        
    ont = "CC",                     # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
head(ogo_CC)
write.table(ogo_CC,"GO_CC_enrich.txt",row.names =FALSE)
dotplot(ogo_CC,title="EnrichmentGO_CC_dot")
barplot(ogo_CC,showCategory=10,title="EnrichmentGO_CC_bar")
plotGOgraph(ogo_CC)
#.rs.restartR()                    # if occur error 
goplot(ogo_CC)
emapplot(ogo_CC,showCategory = 30)
cnetplot(ogo_CC,showCategory = 5)
#==============================================================================================
# GO analysis GSEA(gene set enrichment analysis)    supplement later
#==============================================================================================
#ggo_CC = gseGO(
#    geneList = ,           # id & fold change etc.
#    OrgDb = "org.EcK12.eg.db",        
#    ont = "CC",
#    nPerm = 1000,
#    minGSSize = 100,
#    maxGSSize = 500,
#    pvalueCutoff = 0.05,
#    verbose = FALSE
# )
#=====================================================================================
# KEGG analysis
#=====================================================================================
gk = enrichKEGG(
    gene = entrez_id$ENTREZID,
    keyType = "kegg", 
    organism = 'ecj',         #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05 )
bitr_kegg(
  entrez_id$ENTREZID,
  fromType = "ncbi-geneid",
  toType = 'kegg',
  organism='ecj')    



