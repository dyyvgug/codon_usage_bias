#============================================================================================== 
# 2019-7-29.Author:Dong Yingying.Roughly  observe the translation efficiency.
# Translation efficiency is subtracted from the TPM value quantified by a sample of 
# RNAseq transcripts and the corresponding TPM value of the RIBO-seq transcript(different lab).
# And take out genes that express high expression and high translation level.
#==============================================================================================
library(ggplot2)
library(MASS)
species = "C_elegans_Ensl_WBcel235"
KEGG_spe = "Caenorhabditis elegans"
RNAseq_path = "/RNAseq1/experiment2/SRR1056314_abund.out"
RNA = read.table(paste0("/home/hp/Desktop/other_riboseq/",species,RNAseq_path),sep = "\t",header = T,quote = "")
RNA = RNA[,-c(3,4,5,6,7)]
setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment2/aligned"))
dir.create("ribo_num")
ribo_array = list.files(getwd(),pattern = ".out$")
ribo_array
for (i in ribo_array){
  if(FALSE) # examination
    { 
    ribo = read.table("SRR1804340_abund.out",sep = "\t",header = T,quote = "")
    name = "SRR1804340_abund.out"
    name = sub("^([^.]*).*", "\\1",name)
    name = gsub("_abund","",name)
  }
  ribo = read.table(i,sep = "\t",header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i)
  name = sub("_abund","",name)
  ribo = ribo[,-c(2,3,4,5,6,7)]
  names(ribo) = c("Gene.ID","ribo_FPKM","ribo_TPM")
  RNA_ribo = merge(RNA,ribo,by = "Gene.ID",all = T)
  RNA_ribo[RNA_ribo == 0] <-NA
  RNA_ribo = RNA_ribo[complete.cases(RNA_ribo),]
  #RNA_ribo_num = cbind(RNA_ribo,round(RNA_ribo$TPM / RNA_ribo$ribo_TPM,3))
  ribo_num = round(RNA_ribo$ribo_TPM / RNA_ribo$TPM,3)
  RNA_ribo_num = cbind(RNA_ribo,ribo_num)
  RNA_ribo_num = RNA_ribo_num[order(RNA_ribo_num$ribo_num,decreasing = T),]
  RNA_ribo_num$Gene.ID = gsub("_.*", "", RNA_ribo_num[,1])
  write.table(RNA_ribo_num,file = paste0("./ribo_num/",name,"_riboNum.txt"),
              sep = "\t",quote = F,row.names = F)
  #=================================================================================================================
  # Extract genes encoding only proteins 
  #=================================================================================================================
  only_protein = read.table(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref/CBI_CAI.txt"),header = T)
  only_protein = only_protein[,-c(3,4,5)]
  only_protein_num = merge(RNA_ribo_num,only_protein,by.x = "Gene.ID",by.y = "transcription_id",all = T)
  only_protein_num = only_protein_num[complete.cases(only_protein_num),]
  only_protein_num = only_protein_num[order(only_protein_num$ribo_num,decreasing = T),]
  write.table(only_protein_num,file = paste0("./ribo_num/",name,"_ProtRiboNum.txt"),
              sep = '\t',quote = F,row.names = F)
  q = quantile(only_protein_num$ribo_num,probs = seq(0,1,0.01))
  q
  hiTE = only_protein_num[only_protein_num$ribo_num >= q[99],]
  write.table(hiTE,file = paste0("./ribo_num/",name,"_highTE_all.txt"),
              sep = '\t',quote = F,row.names = F)
  write.table(paste(hiTE$Gene.Name,hiTE$ribo_num,sep = "\t"),
              file = paste0("./ribo_num/",name,"_highTE_geneNAME.txt"),
              sep = '\t',quote = F,row.names = F,col.names = F)
  #==================================================================================================================
  # Observe the correlation between RNAseq and RIBOseq
  #==================================================================================================================
  co_RNA_ri = cor(only_protein_num$TPM,only_protein_num$ribo_TPM)
  co_RNA_ri
  svg(paste0(name,"cor_RNA_ri.svg"))
  #plot(only_protein_num$TPM,only_protein_num$ribo_TPM,log = "xy",main = paste0(name,"cor_RNA_ri  ",co_RNA_ri),
  #     xlab="RNA_TPM",ylab="ri_TPM",pch=19,col=rgb(0,0,100,50,maxColorValue=205))
  ggplot(only_protein_num,aes(x = TPM ,y = ribo_TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_RNA_ri    ","r=",co_RNA_ri))+
    #scale_x_continuous(trans='log10')+
    #scale_y_continuous(trans='log10')+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "blue",size = 0.75)+
    theme(
      panel.background = element_rect(fill = "lightblue",
                                      colour = "lightblue",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"), 
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "white")
    )
  dev.off()
  write.table(co_RNA_ri,file = "cor_RNA_ri.txt",sep = '\t',append = T,quote = FALSE,
              row.names = F, col.names = F)
  #==================================================================================================================
  # Take out genes that express high expression and high translation level
  #==================================================================================================================
  df <- data.frame(only_protein_num$Gene.Name,only_protein_num$TPM,only_protein_num$ribo_TPM)
  gene_id = only_protein_num$Gene.Name
  names(df) = c("Gene_name","TPM","ribo_TPM")
  #tmp <- cor(df$TPM,df$ribo_TPM)
  #tmp[upper.tri(tmp)] <- 0
  #data.new <- df[,!apply(tmp,2,function(x) any(x < 0.6))] #something wrong
  ehE_hT <- subset(x = df,subset = TPM>10^3 & ribo_TPM>10^3,select = c(Gene_name,TPM,ribo_TPM))
  hE_hT <- subset(x = df,subset = TPM>10^3 & ribo_TPM>55,select = c(Gene_name,TPM,ribo_TPM))
  lE_lT <- subset(x =df,subset = TPM<.3 & ribo_TPM<.3,select = c(Gene_name,TPM,ribo_TPM))
  write.table(ehE_hT,file = paste0("./ribo_num/",name,"_ehiE_ht_gene.txt"),sep = "\t",quote = FALSE,
              row.names = F)
  write.table(hE_hT,file = paste0("./ribo_num/",name,"_hiE_ht_gene.txt"),sep = "\t",quote = FALSE,
              row.names = F)
  write.table(lE_lT,file = paste0("./ribo_num/",name,"_lE_lT_gene.txt"),sep = "\t",quote = FALSE,
              row.names = F )
  #======================================================================================================
  # GO and KEGG analyze high translation efficiency genes,high RNA level & high translation level genes,
  # low RNA level & low translation level genes.
  #======================================================================================================
  library(topGO)
  library(clusterProfiler)
  library(pathview)
  library(AnnotationHub)
  require(AnnotationHub)
  hub = AnnotationHub()
  unique(hub$species)
  hub$species[which(hub$species== KEGG_spe)]
  query(hub,KEGG_spe)
  hub[hub$species ==  KEGG_spe &hub$rdataclass == 'OrgDb']
  OrgDb = hub[["AH70577"]]
  keytypes(OrgDb)
  #columns(OrgDb)
  #==============================================================================================
  # Import data and conversion id
  #==============================================================================================
  symbol_id = read.table("gene_id.txt",header=FALSE)      
  symbol_id = as.character(symbol_id$V1)                    
  entrez_id = bitr(symbol_id, 'SYMBOL','ENTREZID', Bacillus.OrgDb)
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
    OrgDb = Bacillus.OrgDb,        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
  head(ogo_all,2)
  write.table(ogo_all,"GO_aLL_enrich.txt",row.names =FALSE)
  svg(filename = "EnrichmentGO_ALL_dot.svg")
  dotplot(ogo_all,title = "EnrichmentGO_all_dot")
  dev.off()
  svg(filename = "EnrichmentGO_ALL_bar.svg")
  barplot(ogo_all,showCategory = 10,title = "EnrichmentGO_all_bar")
  dev.off()
  #*********************************** MF(Molecular Function)*************************************
  ogo_MF = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = Bacillus.OrgDb,        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(ogo_MF)
  write.table(ogo_MF,"GO_MF_enrich.txt",row.names =FALSE)
  svg(filename = "EnrichmentGO_MF_dot.svg")
  dotplot(ogo_MF,title = "EnrichmentGO_MF_dot")
  dev.off()
  barplot(ogo_MF,showCategory = 10,title = "EnrichmentGO_MF_bar")
  plotGOgraph(ogo_MF)
  #.rs.restartR()                    # if occur error 
  goplot(ogo_MF)
  emapplot(ogo_MF,showCategory = 30)
  cnetplot(ogo_MF,showCategory = 5)
  #********************************BP(Biological Process)*******************************
  ogo_BP = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = Bacillus.OrgDb,        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(ogo_BP)
  write.table(ogo_BP,"GO_BP_enrich.txt",row.names =FALSE)
  svg(filename = "EnrichmentGO_BP_dot.svg")
  dotplot(ogo_BP,title = "EnrichmentGO_BP_dot")
  dev.off()
  barplot(ogo_BP,showCategory = 10,title = "EnrichmentGO_BP_bar")
  plotGOgraph(ogo_BP)
  #.rs.restartR()                    # if occur error 
  goplot(ogo_BP)
  emapplot(ogo_BP,showCategory = 30)
  cnetplot(ogo_BP,showCategory = 5)
  #********************************CC(Cellular Component)*******************************
  ogo_CC = enrichGO(
    gene = entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = Bacillus.OrgDb,        
    ont = "CC",                     # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           # other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 1,            
    qvalueCutoff = 1,
    readable = TRUE) 
  head(ogo_CC)
  write.table(ogo_CC,"GO_CC_enrich.txt",row.names =FALSE)
  svg(filename = "EnrichmentGO_CC_dot.svg")
  dotplot(ogo_CC,title = "EnrichmentGO_CC_dot")
  dev.off()
  barplot(ogo_CC,showCategory = 10,title = "EnrichmentGO_CC_bar")
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
  #    OrgDb = Bacillus.OrgDb,        
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
  # symbol to entrez,use bioDBnet,https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
  entrez_id = read.table("symbol_entrez_id.txt",sep = "\t",quote = "",header = T)
  KEGG_id = bitr_kegg(
    entrez_id$Gene.ID,
    fromType = "ncbi-geneid",
    toType = 'kegg',
    organism='bsu')    
  head(KEGG_id)
  write.table(KEGG_id,file = "KEGG_id.txt",sep = '\t',quote = FALSE,
              row.names = FALSE)
  ke = enrichKEGG(
    gene = KEGG_id$kegg,
    keyType = "kegg", 
    organism = 'bsu',         #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
    pAdjustMethod = "BH", 
    pvalueCutoff = 1, 
    qvalueCutoff = 1 )
  head(ke)
  write.table(ke,"KEGG_enrich.txt",row.names =FALSE)
  svg(filename = "KEGG_dot.svg")
  dotplot(ke,showCategory = 10,title="KEGG_dot")
  dev.off()
  svg(filename = "KEGG_bar.svg")
  barplot(ke,showCategory = 10,title="KEGG_bar")
  dev.off()
  #.rs.restartR()                    # if occur error 
  emapplot(ke,showCategory = 30)
  cnetplot(ke,showCategory = 5)
  #browseKEGG(ke, "keggid")          # Mark enriched genes on the pathway map
  
}
