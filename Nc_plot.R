library('ggplot2')
species = "hg38"
species_abb = "Hg"
TPM_file = "/experiment2/SRR8215976.gtf.T"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref"))
codonW = read.table(file = 'CBI_CAI_bycodonW.txt',sep = '\t',header = T,quote = "")
codonW = codonW[-16]
codonW[codonW == '*****'] <-NA
codonW = codonW[complete.cases(codonW),]
codonW$Nc = as.numeric(as.character(codonW$Nc))
#=================================================================================
# Nc-plot
#=================================================================================
# Custom standard curve function
st.fun <- function(gc3){
  2 + gc3 + 29/(gc3**2 + (1-gc3)**2)
}
# Scatter plot + standard curve
svg(file = "Nc_plot.svg")
#jpeg(file = "Nc_plot.jpg",width=1000,height=600)
ggplot(codonW, aes(x = GC3s, y = Nc)) + 
  geom_point(shape = 16, size = 0.1) + 
  stat_function(fun = st.fun, size = 1,color = 'blue') +
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65)) +
  theme(axis.line = element_line(colour = 'black')) +
  labs(title = paste0(species_abb,"_Nc_GC3s"))+
  xlab('GC3s')
dev.off()
#===================================================================================
# Nc and mRNA levels
#===================================================================================
mRNA_level = read.table(paste0("/media/hp/Katniss/DYY/aligned/",species,TPM_file),sep = '\t',header = F,quote = "")
names(mRNA_level) = c("rna_id","gene_name","FPKM","TPM")
Nc_TPM = merge(codonW,mRNA_level,by.x = "transcription_id",by.y = "rna_id",all = T)
Nc_TPM = Nc_TPM[Nc_TPM$TPM > 1,]
Nc_TPM$TPM = log2(Nc_TPM$TPM)
Nc_TPM = Nc_TPM[complete.cases(Nc_TPM),] # Or df = na.omit(df)
if(FALSE) # Multi-line comment
{
  mid<-mean(Nc_TPM$TPM)
  Nc_TPM <- ggplot(Nc_TPM,aes(x = GC3s,y = Nc,color = TPM))+
    geom_point(shape = 16,size = 0.1)+
    stat_function(fun = st.fun,size = 1,color = "black")+
    scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
    scale_color_gradient(low="green", high="red")+
    theme(axis.line = element_line(colour = 'black')) +
    labs(title = paste0(species_abb,"_Nc_GC3s"))+
    xlab('GC3s')
  Nc_TPM
}
sum(Nc_TPM$TPM == 0)
Nc_TPM = Nc_TPM[Nc_TPM$TPM > 0,] 
#Nc_TPM[] <- lapply(Nc_TPM, function(x) ifelse(x > 500, 500, x))
svg(file = "Nc_TPM_ex2.svg")
Nc_TPM_plot <- ggplot(Nc_TPM,aes(x = GC3s,y = Nc,color = TPM))+
  geom_point(shape = 16,size = 1)+
  stat_function(fun = st.fun,size = 1,color = "black")+
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  scale_color_gradientn(colours = rainbow(5))+
  theme(axis.line = element_line(colour = 'black')) +
  labs(title = paste0(species_abb,"_Nc_GC3s"))+
  xlab('GC3s')+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"))
Nc_TPM_plot
dev.off()


