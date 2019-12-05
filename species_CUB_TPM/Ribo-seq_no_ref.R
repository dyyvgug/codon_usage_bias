#========================================================================================================
# 2019-12-4.Author:Yingying Dong.Correlation analysis of modified CAI and various 
#  gene expression levels in Ribo-seq all samples.The weight from ribosomal protein genes.
#========================================================================================================
library(getopt)
library(ggplot2)
library(MASS)
library(scales)
command=matrix(c("species","s",1,"character",
                 "experiment","e",1,"numeric",
                 "species_abb","A",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$species) || is.null(args$experiment) || is.null(args$species_abb)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

species = args$species
exp = args$experiment
speA_exp = args$species_abb

if(FALSE){
  species = "C_elegans_Ensl_WBcel235"
  exp = "3"
  speA_exp = "Ce"
}
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                        species, "/ref/CBI_CAI.txt", sep = ""), sep = "\t", header = T,quote = '')
cai = data.frame(cai$transcription_id,cai$CAI)
dir.create(paste0("picture_riboseq_byGLOBAL",exp))
dir.create(paste0("correlation_riboseq_byGLOBAL",exp))

setwd(paste0("~/Desktop/other_riboseq/",species,"/experiment",exp,"/aligned_ri/"))
ribo_array = list.files(getwd(),pattern = ".out$")
ribo_array
for (i in ribo_array) {
  if (FALSE){
    gtf = read.table("SRR1804340_abund.out", sep = "\t", header=T,quote = "")
    name = "SRR1804340"
  }
  gtf = read.table(i, sep = "\t", header = T,quote = "",fill = T)
  name = sub("^([^.]*).*", "\\1",i) 
  name = sub("_abund","",name)
  
  gtf = gtf[,-c(2,3,4,5,6,7)]
  names(cai) = c("gene_id","CAI")
  cai$gene_id = sub(" ","",cai$gene_id)
  df = merge(cai, gtf, by.x = "gene_id",by.y = "Gene.ID", all = T)
  df[df == 0] <- NA
  df = df[complete.cases(df),] # Or df = na.omit(df)

  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_riboseq_byGLOBAL",exp,"/",name,"_CAI_TPM.svg"))
  CAI_cor = cor(df$CAI,df$TPM)
  CAI_cor
  CAI_cor_p = cor.test(df$CAI,df$TPM)
  
  p <- ggplot(df,aes(x = df$CAI ,y = df$TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_mCAI_TPM    ","r=",round(CAI_cor,3),"  p=",round(CAI_cor_p$p.value,5)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("CAI(no reference set) value")+
    ylab("Ribo-seq (TPM)")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p)
  dev.off()
  write.table(CAI_cor,file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_riboseq_byGLOBAL",exp,"/","CAI_cor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
}
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_riboseq_byGLOBAL",exp,"/" ))
a = read.table("CAI_cor.txt",sep = '\t',header = F)

aveCAI = mean(a$V1)
aveCAI

write.table(aveCAI,file = "mean_CAI.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM")

