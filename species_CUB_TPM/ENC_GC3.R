library('ggplot2')
library('reshape2')
library('corrplot')
library('readr')
library('ggpubr')
setwd("/media/hp/disk2/linchenghao/hg38/ref")
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
  labs(title = "Hg_Nc_GC3s")+
  xlab('GC3s')
dev.off()
