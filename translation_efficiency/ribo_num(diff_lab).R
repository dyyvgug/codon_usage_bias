#============================================================================================== 
# 2019-7-29.Modified date:2019-9-5.Author:Dong Yingying.Roughly  observe the translation efficiency.
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
  threshhold <- 1
  df = subset(df, df[,2] > threshhold) 
  df = subset(df, df[,3] > threshhold)
  #tmp <- cor(df$TPM,df$ribo_TPM)
  #tmp[upper.tri(tmp)] <- 0
  #data.new <- df[,!apply(tmp,2,function(x) any(x < 0.6))] #something wrong
  qRNA = quantile(df$TPM,probs = seq(0,1,0.01))
  qRNA
  hiRNA = df[df$TPM > qRNA[97],]    #TOP 4%
  loRNA = df[df$TPM < qRNA[5],]     #BOTTOM 4%
  qRIBO = quantile(df$ribo_TPM,probs = seq(0,1,0.01))
  qRIBO
  hiRIBO = df[df$ribo_TPM > qRIBO[97],]    #TOP 4%
  other_ribo = df[df$ribo_TPM < qRIBO[97],]    #In order to compare the differences between other ribosomal genes and high expression of high translation ribosomal genes. 
  loRIBO = df[df$ribo_TPM < qRIBO[5],]     #BOTTOM 4%
  hE_hT <- merge(hiRNA,hiRIBO,all = F)  #TOP4% mRNA level and top4% RIBOseq level,intersection.
  lE_lT <- merge(loRNA,loRIBO,all = F)
  ehE_hT <- subset(x = df,subset = TPM>10^3 & ribo_TPM>10^3,select = c(Gene_name,TPM,ribo_TPM))
  #hE_hT <- subset(x = df,subset = TPM>10^3 & ribo_TPM>55,select = c(Gene_name,TPM,ribo_TPM))
  #lE_lT <- subset(x =df,subset = TPM<.3 & ribo_TPM<.3,select = c(Gene_name,TPM,ribo_TPM))
  write.table(hE_hT_def,file = paste0("./ribo_num/",name,"_hE_ht_def_gene.txt"),sep = "\t",quote = FALSE,
              row.names = F)
  write.table(lE_lT_def,paste0("./ribo_num/",name,"_lE_lT_def_gene.txt"),sep = "\t",quote = FALSE,
              row.names = F)
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
  ##high translation efficiency genes conversion id 
  hiTE_symbol_id =  hiTE$Gene.Name     
  hiTE_symbol_id = as.character(hiTE_symbol_id)                    
  hiTE_entrez_id = bitr(hiTE_symbol_id, 'SYMBOL','ENTREZID',OrgDb)
  head(hiTE_entrez_id,2)
  write.table(hiTE_entrez_id,file = paste0(name,"_hiTE_sym_entrez.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  TE_only_entrezID = hiTE_entrez_id[-1]
  write.table(TE_only_entrezID,file = paste0(name,"_hiTE_ENTREZID.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  ##extreme high RNA expression level and high translation level genes
  ehE_hT_symID = ehE_hT$Gene_name
  ehE_hT_symID = as.character(ehE_hT_symID)
  ehE_hT_entID = bitr(ehE_hT_symID,'SYMBOL','ENTREZID',OrgDb)
  head(ehE_hT_entID,2)
  write.table(ehE_hT_entID,file = paste0(name,"_ehE_hT_sym_entrez.txt"),sep = '\t',quote = F,row.names = F)
  ehET_only_entID = ehE_hT_entID[-1]
  write.table(ehET_only_entID,file = paste0(name,"_ehET_ENTREZID.txt"),sep = '\t',quote = F,row.names = F)
  ## high RNA expression level and high translation level genes
  hE_hT_symID = hE_hT$Gene_name
  hE_hT_symID = as.character(hE_hT_symID)
  hE_hT_entID = bitr(hE_hT_symID,'SYMBOL','ENTREZID',OrgDb)
  head(hE_hT_entID,2)
  write.table(hE_hT_entID,file = paste0(name,"_hE_hT_sym_entrez.txt"),sep = '\t',quote = F,row.names = F)
  hET_only_entID = hE_hT_entID[-1]
  write.table(hET_only_entID,file = paste0(name,"_hET_ENTREZID.txt"),sep = '\t',quote = F,row.names = F)
  ## low RNA expression level and low translation level genes
  lE_lT_symID = lE_lT$Gene_name
  lE_lT_symID = as.character(lE_lT_symID)
  lE_lT_entID = bitr(lE_lT_symID,'SYMBOL','ENTREZID',OrgDb)
  head(lE_lT_entID,2)
  write.table(lE_lT_entID,file = paste0(name,"_lE_lT_sym_entrez.txt"),sep = '\t',quote = F,row.names = F)
  hET_only_entID = lE_lT_entID[-1]
  write.table(hET_only_entID,file = paste0(name,"_lET_ENTREZID.txt"),sep = '\t',quote = F,row.names = F)
  #==============================================================================================
  # GO analysis ORA(over-representation analysis)
  #==============================================================================================
  #**************all(Biological Process,Cellular Component,Molecular Function)*******************
  if(FALSE) # something wrong,modify it later
  {
  GO_all <- {
    go_all = enrichGO(
      gene = y, 
      keyType = "ENTREZID",
      OrgDb = OrgDb,        
      ont = "ALL",                    # Can also be a kind of CC,BP,MF
      pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
      pvalueCutoff = 0.05,            
      qvalueCutoff = 0.2,
      readable = TRUE)                # ID to Symbol,easy to read
    #head(paste0(go_all,2))
    write.table(go_all,file = paste0(name,"_",x,"_aLL_enrich.txt"),row.names =FALSE)
    svg(filename = paste0(name,"_",x,"_ALLdot.svg"))
    dotplot(go_all,title = "EnrichmentGO_all_dot")
    dev.off()
    svg(filename = paste0(name,"_",x,"_ALLbar.svg"))
    barplot(go_all,showCategory = 10,title = "EnrichmentGO_all_bar")
    dev.off()
  }
  GO_all("hiTE",hiTE_entrez_id$ENTREZID)  ## high TE
  }
  ## high TE
  hiTE_go_all = enrichGO(
    gene = hiTE_entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
  head(hiTE_go_all,2)
  write.table(hiTE_go_all,file = paste0(name,"_hiTE_aLL_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hiTE_ALLdot.svg"))
  dotplot(hiTE_go_all,title = "EnrichmentGO_all_dot")
  dev.off()
  svg(filename = paste0(name,"_hiTE_ALLbar.svg"))
  barplot(hiTE_go_all,showCategory = 10,title = "EnrichmentGO_all_bar")
  dev.off()
  ## extreme high RNA expression level and high translation level genes
  ehE_hT_go_all = enrichGO(
    gene = ehE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
  head(ehE_hT_go_all,2)
  write.table(ehE_hT_go_all,file = paste0(name,"_ehE_hT_aLL_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_ehE_hT_ALLdot.svg"))
  dotplot(ehE_hT_go_all,title = "EnrichmentGO_all_dot")
  dev.off()
  svg(filename = paste0(name,"_ehE_hT_ALLbar.svg"))
  barplot(ehE_hT_go_all,showCategory = 10,title = "EnrichmentGO_all_bar")
  dev.off()
  ## high RNA expression level and high translation level genes
  hE_hT_go_all = enrichGO(
    gene = hE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
  head(hE_hT_go_all,2)
  write.table(hE_hT_go_all,file = paste0(name,"_hE_hT_aLL_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hE_hT_ALLdot.svg"))
  dotplot(hE_hT_go_all,title = "EnrichmentGO_all_dot")
  dev.off()
  svg(filename = paste0(name,"_hE_hT_ALLbar.svg"))
  barplot(hE_hT_go_all,showCategory = 10,title = "EnrichmentGO_all_bar")
  dev.off()
  ## low RNA expression level and low translation level genes
  lE_lT_go_all = enrichGO(
    gene = lE_lT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "ALL",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE)                # ID to Symbol,easy to read
  head(lE_lT_go_all,2)
  write.table(lE_lT_go_all,file = paste0(name,"_lE_lT_aLL_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_lE_lT_ALLdot.svg"))
  dotplot(lE_lT_go_all,title = "EnrichmentGO_all_dot")
  dev.off()
  svg(filename = paste0(name,"_lE_lT_ALLbar.svg"))
  barplot(lE_lT_go_all,showCategory = 10,title = "EnrichmentGO_all_bar")
  dev.off()
  #*********************************** MF(Molecular Function)*************************************
  ## high TE
   hiTE_go_MF = enrichGO(
    gene = hiTE_entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(hiTE_go_MF)
  write.table(hiTE_go_MF,file = paste0(name,"_hiTE_MF_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hiTE_MFdot.svg"))
  dotplot(hiTE_go_MF,title = "EnrichmentGO_MF_dot")
  dev.off()
  svg(filename = paste0(name,"_hiTE_MFbar.svg"))
  barplot(hiTE_go_MF,showCategory = 10,title = "EnrichmentGO_MF_bar")
  dev.off()
  #plotGOgraph(hiTE_go_MF)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hiTE_MFgoplot.svg"))
  goplot(hiTE_go_MF)
  dev.off()
  #emapplot(hiTE_go_MF,showCategory = 30)
  #cnetplot(hiTE_go_MF,showCategory = 5)
  ## extreme high RNA expression level and high translation level genes
  ehE_hT_go_MF = enrichGO(
    gene = ehE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(ehE_hT_go_MF)
  write.table(ehE_hT_go_MF,file = paste0(name,"_ehE_hT_MF_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_ehE_hT_MFdot.svg"))
  dotplot(ehE_hT_go_MF,title = "EnrichmentGO_MF_dot")
  dev.off()
  svg(filename = paste0(name,"_ehE_hT_MFbar.svg"))
  barplot(ehE_hT_go_MF,showCategory = 10,title = "EnrichmentGO_MF_bar")
  dev.off()
  #plotGOgraph(ehE_hT_go_MF)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_ehE_hT_MFgoplot.svg"))
  goplot(ehE_hT_go_MF)
  dev.off()
  ## high RNA expression level and high translation level genes
  hE_hT_go_MF = enrichGO(
    gene = hE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(hE_hT_go_MF)
  write.table(hE_hT_go_MF,file = paste0(name,"_hE_hT_MF_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hE_hT_MFdot.svg"))
  dotplot(hE_hT_go_MF,title = "EnrichmentGO_MF_dot")
  dev.off()
  svg(filename = paste0(name,"_hE_hT_MFbar.svg"))
  barplot(hE_hT_go_MF,showCategory = 10,title = "EnrichmentGO_MF_bar")
  dev.off()
  #plotGOgraph(hE_hT_go_MF)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hE_hT_MFgoplot.svg"))
  goplot(hE_hT_go_MF)
  dev.off()
  ## low RNA expression level and low translation level genes
  lE_lT_go_MF = enrichGO(
    gene = lE_lT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "MF",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 1,            
    qvalueCutoff = 1,
    readable = TRUE) 
  head(lE_lT_go_MF)
  write.table(lE_lT_go_MF,file = paste0(name,"_lE_lT_MF_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_lE_lT_MFdot.svg"))
  dotplot(lE_lT_go_MF,title = "EnrichmentGO_MF_dot")
  dev.off()
  svg(filename = paste0(name,"_lE_lT_MFbar.svg"))
  barplot(lE_lT_go_MF,showCategory = 10,title = "EnrichmentGO_MF_bar")
  dev.off()
  #plotGOgraph(lE_lT_go_MF)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_lE_lT_MFgoplot.svg"))
  goplot(lE_lT_go_MF)
  dev.off()
  #********************************BP(Biological Process)*******************************
  ## high TE
  hiTE_go_BP = enrichGO(
    gene = hiTE_entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(hiTE_go_BP)
  write.table(hiTE_go_BP,file = paste0(name,"_hiTE_BP_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hiTE_BPdot.svg"))
  dotplot(hiTE_go_BP,title = "EnrichGO_BP_dot")
  dev.off()
  svg(filename = paste0(name,"_hiTE_BPbar.svg"))
  barplot(hiTE_go_BP,showCategory = 10,title = "EnrichmentGO_BP_bar")
  dev.off()
  #plotGOgraph(hiTE_go_BP)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hiTE_BPgoplot.svg"))
  goplot(hiTE_go_BP)
  dev.off()
  #emapplot(hiTE_go_BP,showCategory = 30)
  #cnetplot(hiTE_go_BP,showCategory = 5)
  ## extreme high
  ehE_hT_go_BP = enrichGO(
    gene = ehE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(ehE_hT_go_BP)
  write.table(ehE_hT_go_BP,file = paste0(name,"_ehE_hT_BP_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_ehE_hT_BPdot.svg"))
  dotplot(ehE_hT_go_BP,title = "EnrichGO_BP_dot")
  dev.off()
  svg(filename = paste0(name,"_ehE_hT_BPbar.svg"))
  barplot(ehE_hT_go_BP,showCategory = 10,title = "EnrichmentGO_BP_bar")
  dev.off()
  #plotGOgraph(ehE_hT_go_BP)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_ehE_hT_BPgoplot.svg"))
  goplot(ehE_hT_go_BP)
  dev.off()
  ## high expression and high translation
  hE_hT_go_BP = enrichGO(
    gene = hE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(hE_hT_go_BP)
  write.table(hE_hT_go_BP,file = paste0(name,"_hE_hT_BP_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hE_hT_BPdot.svg"))
  dotplot(hE_hT_go_BP,title = "EnrichGO_BP_dot")
  dev.off()
  svg(filename = paste0(name,"_hE_hT_BPbar.svg"))
  barplot(hE_hT_go_BP,showCategory = 10,title = "EnrichmentGO_BP_bar")
  dev.off()
  #plotGOgraph(hE_hT_go_BP)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hE_hT_BPgoplot.svg"))
  goplot(hE_hT_go_BP)
  dev.off()
  ## low expression and low translation
  lE_lT_go_BP = enrichGO(
    gene = lE_lT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "BP",                    # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(lE_lT_go_BP)
  write.table(lE_lT_go_BP,file = paste0(name,"_lE_lT_BP_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_lE_lT_BPdot.svg"))
  dotplot(lE_lT_go_BP,title = "EnrichGO_BP_dot")
  dev.off()
  svg(filename = paste0(name,"_lE_lT_BPbar.svg"))
  barplot(lE_lT_go_BP,showCategory = 10,title = "EnrichmentGO_BP_bar")
  dev.off()
  #plotGOgraph(lE_lT_go_BP)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_lE_lT_BPgoplot.svg"))
  goplot(lE_lT_go_BP)
  dev.off()
  #********************************CC(Cellular Component)*******************************
  ## high TE
  hiTE_go_CC = enrichGO(
    gene = hiTE_entrez_id$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "CC",                     # Can also be a kind of CC,BP,MF
    pAdjustMethod = "BH",           # other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 1,            
    qvalueCutoff = 1,
    readable = TRUE) 
  head(hiTE_go_CC)
  write.table(hiTE_go_CC,paste0(name,"_hiTE_CC_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hiTE_CCdot.svg"))
  dotplot(hiTE_go_CC,title = "EnrichmentGO_CC_dot")
  dev.off()
  svg(filename = paste0(name,"_hiTE_CCbar.svg"))
  barplot(hiTE_go_CC,showCategory = 10,title = "EnrichmentGO_CC_bar")
  dev.off()
  #plotGOgraph(hiTE_go_CC)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hiTE_CCgoplot.svg"))
  goplot(hiTE_go_CC)
  dev.off()
  #emapplot(hiTE_go_CC,showCategory = 30)
  #cnetplot(hiTE_go_CC,showCategory = 5)
  ## extreme high
  ehE_hT_go_CC = enrichGO(
    gene = ehE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "CC",                    # Can also be a kind of CC,CC,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(ehE_hT_go_CC)
  write.table(ehE_hT_go_CC,file = paste0(name,"_ehE_hT_CC_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_ehE_hT_CCdot.svg"))
  dotplot(ehE_hT_go_CC,title = "EnrichGO_CC_dot")
  dev.off()
  svg(filename = paste0(name,"_ehE_hT_CCbar.svg"))
  barplot(ehE_hT_go_CC,showCategory = 10,title = "EnrichmentGO_CC_bar")
  dev.off()
  #plotGOgraph(ehE_hT_go_CC)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_ehE_hT_CCgoplot.svg"))
  goplot(ehE_hT_go_CC)
  dev.off()
  ## high expression and high translation
  hE_hT_go_CC = enrichGO(
    gene = hE_hT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "CC",                    # Can also be a kind of CC,CC,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(hE_hT_go_CC)
  write.table(hE_hT_go_CC,file = paste0(name,"_hE_hT_CC_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hE_hT_CCdot.svg"))
  dotplot(hE_hT_go_CC,title = "EnrichGO_CC_dot")
  dev.off()
  svg(filename = paste0(name,"_hE_hT_CCbar.svg"))
  barplot(hE_hT_go_CC,showCategory = 10,title = "EnrichmentGO_CC_bar")
  dev.off()
  #plotGOgraph(hE_hT_go_CC)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_hE_hT_CCgoplot.svg"))
  goplot(hE_hT_go_CC)
  dev.off()
  ## low expression and low translation
  lE_lT_go_CC = enrichGO(
    gene = lE_lT_entID$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = OrgDb,        
    ont = "CC",                    # Can also be a kind of CC,CC,MF
    pAdjustMethod = "BH",           #other correction methods: holm,hochberg,hommel,bonferroni,BH,BY,fdr,none
    pvalueCutoff = 0.05,            
    qvalueCutoff = 0.2,
    readable = TRUE) 
  head(lE_lT_go_CC)
  write.table(lE_lT_go_CC,file = paste0(name,"_lE_lT_CC_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_lE_lT_CCdot.svg"))
  dotplot(lE_lT_go_CC,title = "EnrichGO_CC_dot")
  dev.off()
  svg(filename = paste0(name,"_lE_lT_CCbar.svg"))
  barplot(lE_lT_go_CC,showCategory = 10,title = "EnrichmentGO_CC_bar")
  dev.off()
  #plotGOgraph(lE_lT_go_CC)
  #.rs.restartR()                    # if occur error 
  svg(filename = paste0(name,"_lE_lT_CCgoplot.svg"))
  goplot(lE_lT_go_CC)
  dev.off()
  #=====================================================================================
  # KEGG analysis
  #=====================================================================================
  ## high TE 
  hiTE_KEGG_id = bitr_kegg(
    hiTE_entrez_id$ENTREZID,
    fromType = "ncbi-geneid",
    toType = 'kegg',
    organism='cel')           #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
  head(hiTE_KEGG_id)
  write.table(hiTE_KEGG_id,file = paste0(name,"_hiTE_KEGGid.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  hiTE_ke = enrichKEGG(
    gene = hiTE_KEGG_id$kegg,
    keyType = "kegg", 
    organism = 'cel',         
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.2 )
  head(hiTE_ke)
  write.table(hiTE_ke,paste0(name,"_hiTE_KEGG_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hiTE_KEGGdot.svg"))
  dotplot(hiTE_ke,showCategory = 10,title="KEGG_dot")
  dev.off()
  svg(filename = paste0(name,"_hiTE_KEGGbar.svg"))
  barplot(hiTE_ke,showCategory = 10,title="KEGG_bar")
  dev.off()
  #.rs.restartR()                    # if occur error 
  emapplot(hiTE_ke,showCategory = 30)
  cnetplot(hiTE_ke,showCategory = 5)
  #browseKEGG(ke, "keggid")          # Mark enriched genes on the pathway map
  ## extreme high RNA expression level and high translation level genes
  ehE_hT_KEGG_id = bitr_kegg(
    ehE_hT_entID$ENTREZID,
    fromType = "ncbi-geneid",
    toType = 'kegg',
    organism='cel')           #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
  head(ehE_hT_KEGG_id)
  write.table(ehE_hT_KEGG_id,file = paste0(name,"_ehE_hT_KEGGid.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  ehE_hT_ke = enrichKEGG(
    gene = ehE_hT_KEGG_id$kegg,
    keyType = "kegg", 
    organism = 'cel',         
    pAdjustMethod = "BH", 
    pvalueCutoff = 1, 
    qvalueCutoff = 1 )
  head(ehE_hT_ke)
  write.table(ehE_hT_ke,paste0(name,"_ehE_hT_KEGG_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_ehE_hT_KEGGdot.svg"))
  dotplot(ehE_hT_ke,showCategory = 10,title="KEGG_dot")
  dev.off()
  svg(filename = paste0(name,"_ehE_hT_KEGGbar.svg"))
  barplot(ehE_hT_ke,showCategory = 10,title="KEGG_bar")
  dev.off()
  #.rs.restartR()                    # if occur error 
  emapplot(ehE_hT_ke,showCategory = 30)
  cnetplot(ehE_hT_ke,showCategory = 5)
  ## high RNA expression level and high translation level genes
  hE_hT_KEGG_id = bitr_kegg(
    hE_hT_entID$ENTREZID,
    fromType = "ncbi-geneid",
    toType = 'kegg',
    organism='cel')           #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
  head(hE_hT_KEGG_id)
  write.table(hE_hT_KEGG_id,file = paste0(name,"_hE_hT_KEGGid.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  hE_hT_ke = enrichKEGG(
    gene = hE_hT_KEGG_id$kegg,
    keyType = "kegg", 
    organism = 'cel',         
    pAdjustMethod = "BH", 
    pvalueCutoff = 1, 
    qvalueCutoff = 1 )
  head(hE_hT_ke)
  write.table(hE_hT_ke,paste0(name,"_hE_hT_KEGG_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_hE_hT_KEGGdot.svg"))
  dotplot(hE_hT_ke,showCategory = 10,title="KEGG_dot")
  dev.off()
  svg(filename = paste0(name,"_hE_hT_KEGGbar.svg"))
  barplot(hE_hT_ke,showCategory = 10,title="KEGG_bar")
  dev.off()
  #.rs.restartR()                    # if occur error 
  emapplot(hE_hT_ke,showCategory = 30)
  cnetplot(hE_hT_ke,showCategory = 5)
  ## low RNA expression level and low translation level genes
  lE_lT_KEGG_id = bitr_kegg(
    lE_lT_entID$ENTREZID,
    fromType = "ncbi-geneid",
    toType = 'kegg',
    organism='cel')           #abbreviation https://www.genome.jp/kegg/catalog/org_list.html
  head(lE_lT_KEGG_id)
  write.table(lE_lT_KEGG_id,file = paste0(name,"_lE_lT_KEGGid.txt"),sep = '\t',quote = FALSE,
              row.names = FALSE)
  lE_lT_ke = enrichKEGG(
    gene = lE_lT_KEGG_id$kegg,
    keyType = "kegg", 
    organism = 'cel',         
    pAdjustMethod = "BH", 
    pvalueCutoff = 1, 
    qvalueCutoff = 1 )
  head(lE_lT_ke)
  write.table(lE_lT_ke,paste0(name,"_lE_lT_KEGG_enrich.txt"),row.names =FALSE)
  svg(filename = paste0(name,"_lE_lT_KEGGdot.svg"))
  dotplot(lE_lT_ke,showCategory = 10,title="KEGG_dot")
  dev.off()
  svg(filename = paste0(name,"_lE_lT_KEGGbar.svg"))
  barplot(lE_lT_ke,showCategory = 10,title="KEGG_bar")
  dev.off()
  #.rs.restartR()                    # if occur error 
  emapplot(lE_lT_ke,showCategory = 30)
  cnetplot(lE_lT_ke,showCategory = 5)
}
