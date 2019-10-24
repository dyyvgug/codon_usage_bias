spec = "hg19"
ex = 2
setwd(paste0("/media/hp/Katniss/DYY/aligned/", spec, "/experiment", ex))

cbi_cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                        spec, "/ref/CBI_CAI.txt", sep = ""), sep = "\t", header = T)
#setwd(paste("/media/hp/Katniss/DYY/aligned/",species,sep = ""))
#gtf_array = list.files(getwd(),pattern = "\\d.+[txt]$") 
#gtf_array

dir.create(paste0("/media/hp/disk1/DYY/reference/annotation/", spec,"/picture"))

draw.dotplot = function(species, exp, sample)    {

  gtf = read.table(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment",exp,"/", 
                          sample), sep = "\t", header=F)
  name = sub(".gtf.T","",gtf)
  names(gtf) = c("transcription_id", "gene_id", "FPKM", "TPM")
  df = merge(cbi_cai, gtf, by.x = "transcription_id", by.y = "transcription_id", all = T)
  #df[is.na(df)] = 0
  df = df[complete.cases(df),] # Or df = na.omit(df)
  df = df[df$length > 0,]
  
  GC_cor = cor.test(df[,"GC"],log2(df$TPM+1))

  
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", spec,"/picture/", sample, "_GC_TPM.jpg"))
  plot(df$GC,log2(df$TPM+1),ylab="TPM",
       xlab="GC",pch=19,col=rgb(0,0,100,50,maxColorValue=255),
       main = paste(species," cor =",round(GC_cor$estimate,3)) )
  
  dev.off()
  
  #CAI_cor = cor(df[,"CAI"],log2(df$TPM+1))
  cor = cor.test(df$CAI, log2(df$FPKM+1))
  
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", spec,"/picture/", sample, "_CAI_TPM.jeg"))
  plot(df[,"CAI"],log2(df$TPM+1),ylab="TPM", 
       xlab="CAI",pch=19,col=rgb(0,0,100,50,maxColorValue=255),
        main = paste(species," cor =",round(cor$estimate,3)))
  dev.off()
  
  
  CBI_cor = cor(df[,"CBI"],log2(df$TPM+1))
  CBI_cor
  
  jpeg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", spec,"/picture/", sample, "_CBI_TPM.jpg"))
  
  plot(df[,"CBI"],log2(df$TPM+1),ylab="TPM", 
       xlab="CBI",pch=19,col=rgb(0,0,100,50,maxColorValue=255),
       main = paste(species," cor =",round(CBI_cor,3)))
  #write.table(paste(GC_cor,CAI_cor,CBI_cor,sep = '\t'),file = "97correlation",
              #quote = FALSE,row.names = "correlation", col.names = "\tGC_TPM\tCAI_TPM\tCBI_TPM")
  dev.off()
}

#setwd(paste0("/media/hp/Katniss/DYY/aligned/", spec, "/experiment", ex))
list.files(getwd())
#t = list.files(getwd(),pattern = "SRR\\d.+[gtf].T$")
t = list.files(getwd(),pattern = "SRR\\d.+.T.txt$")
t

species = spec
exp = ex
sample = s


for (s in t)  {
  draw.dotplot(spec, ex, s)
  
}

