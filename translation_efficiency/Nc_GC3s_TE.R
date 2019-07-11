#====================================================================================== 
# 2019-7-11.Author:Dong Yingying.Observe the relationship between translation efficiency 
# and codon usage bais index --Nc value & GC3s value.
#======================================================================================
library('ggplot2')
species = "C_elegans_Ensl_WBcel235"
species_abb = "Ce"
TEfile = "SRR1804340_ProtRiboNum.txt"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/ref"))
codonW = read.table(file = 'CBI_CAI_bycodonW.txt',sep = '\t',header = T,quote = "")
codonW = codonW[-16]
codonW[codonW == '*****'] <-NA
codonW = codonW[complete.cases(codonW),]
codonW$Nc = as.numeric(as.character(codonW$Nc))
#=================================================================================
# Drawing Nc-plot
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
# Nc and TE(translation efficiency) levels
#===================================================================================
TE_level = read.table(TEfile,sep = '\t',header = T,quote = "")
Nc_TPM_TE = merge(codonW,TE_level,by.x = "transcription_id",by.y = "Gene.ID",all = T)
Nc_TPM_TE = Nc_TPM_TE[Nc_TPM_TE$TPM > 1,]
Nc_TPM_TE = Nc_TPM_TE[complete.cases(Nc_TPM_TE),] # Or df = na.omit(df)
Nc_TPM_TE$TPM = log2(Nc_TPM_TE$TPM)
Nc_TPM_TE$ribo_TPM = log2(Nc_TPM_TE$ribo_TPM)
Nc_TPM_TE <- Nc_TPM_TE[!duplicated(Nc_TPM_TE[, c("TPM", "ribo_TPM")]), ] # Remove duplicate lines
#=========================================================================================
# Transcript Quantitative TPM Value from Mapping of SRR1804340 Samples RNAseq 
#  vs. Nc, GC3s of Species Genome
#=========================================================================================
svg(file = "Nc_RNA_TPM.svg")
Nc_RNA_TPM_plot <- ggplot(Nc_TPM_TE,aes(x = GC3s,y = Nc,color = TPM))+
  geom_point(shape = 16,size = 0.5)+
  stat_function(fun = st.fun,size = 1,color = "black")+
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  scale_color_gradientn(colours = rainbow(7))+
  labs(title = paste0(species_abb,"_Nc_RNA_TPM"))+
  xlab('GC3s')+
  theme_bw() +       
  theme(axis.line = element_line(colour = "black"))
Nc_RNA_TPM_plot
dev.off()
#=========================================================================================
# Protected Transcript Quantitative TPM Value from Mapping of SRR1804340 Samples Ribo-seq 
#  vs. Nc, GC3s of Species Genome
#=========================================================================================
svg(file = "Nc_RNA_riTPM.svg")
Nc_RNA_riTPM_plot <- ggplot(Nc_TPM_TE,aes(x = GC3s,y = Nc,color = ribo_TPM))+
  geom_point(shape = 16,size = 0.5)+
  stat_function(fun = st.fun,size = 1,color = "black")+
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  scale_color_gradientn(colours = rainbow(7))+
  labs(title = paste0(species_abb,"_Nc_RNA_riTPM"))+
  xlab('GC3s')+
  theme_bw() +       
  theme(axis.line = element_line(colour = "black"))
Nc_RNA_riTPM_plot
dev.off()
#=========================================================================================
#  Roughly calculated translation efficiency vs. Nc, GC3s of Species Genome
#=========================================================================================
df <- Nc_TPM_TE[c(1,9,10,21)]
df <- cbind(df[1], apply(df[2:4],2, function(x) ifelse(x > 200, 200, x)))
#df = df[df$ribo_num_log2 > 1,]
#df$ribo_num_log2 = log2(df$ribo_num+1)
hist(df$ribo_num)
df = df[df$ribo_num > 1,]
svg(file = "Nc_RNA_TE.svg")
Nc_TE_plot <- ggplot(df,aes(x = GC3s,y = Nc,color = ribo_num))+
  geom_point(shape = 16,size = 1)+
  stat_function(fun = st.fun,size = 1,color = "black")+
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(title = paste0(species_abb,"_Nc_TE"))+
  xlab('GC3s')+
  theme_bw() +       # remove background
  theme(axis.line = element_line(colour = "black"))
Nc_TE_plot
dev.off()

if(FALSE) # Multi-line comment
{
  mid<-mean(Nc_TPM_TE$TPM)
  sum(Nc_TPM$TPM == 0)
  Nc_TPM = Nc_TPM[Nc_TPM$TPM > 0,] 
  Nc_TPM[] <- lapply(Nc_TPM, function(x) ifelse(x > 500, 500, x))
  svg(file = "Nc_TPM_ex1.svg")
  Nc_TPM_plot <- ggplot(Nc_TPM,aes(x = GC3s,y = Nc,color = TPM))+
  geom_point(shape = 16,size = 1)+
  stat_function(fun = st.fun,size = 1,color = "black")+
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+       # Three colours gradientn
  theme(axis.line = element_line(colour = 'black')) +
  labs(title = paste0(species_abb,"_Nc_GC3s"))+
  xlab('GC3s')+
  theme_bw() +       # remove background
  theme(axis.line = element_line(colour = "black"))
  Nc_TPM_plot
  dev.off()
  Nc_TPM_plot + theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  ) # lightblue background
}

