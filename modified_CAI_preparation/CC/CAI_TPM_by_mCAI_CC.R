#========================================================================================================
# 2019-10-28.Modified date:2019-10-29.Author:Yingying Dong.Correlation analysis of modified CAI and 
#  global gene expression in all samples.
#========================================================================================================
library(getopt)
library(ggplot2)
library(MASS)
library(scales)
command=matrix(c("species","s",1,"character",
                 "experiment","e",1,"numeric",
                 "species_abb_num","A",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$species) || is.null(args$experiment) || is.null(args$species_abb_num)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}
species = args$species
exp = args$experiment
speA_exp = args$species_abb_num

#species = "C_elegans_Ensl_WBcel235"
#exp = "1"
#speA_exp = "Ce"

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/", species))
cai  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                            species, "/ref/",speA_exp,"_mCAI_CC.txt", sep = ""), sep = "\t", header = T)
cai$gene_id = gsub(">gene-", "",cai$gene_id)
dir.create(paste0("picture_bymCAI_CC",exp))
dir.create(paste0("correlation_bymCAI_CC",exp))
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species,"/experiment",exp,"/"))
gtf_array = list.files(getwd(),pattern = "[SE]RR\\d.+out$") 
gtf_array
#graphics.off()
for (i in gtf_array) {
 if (FALSE){
   gtf = read.table("SRR7160566_abund.out", sep = "\t", header=T,quote = "")
   name = "SRR7160566"
 }
  
  gtf = read.table(i, sep = "\t", header = T,quote = "")
  name = sub("^([^.]*).*", "\\1",i) 
  name = sub("_abund","",name)
  
  gtf = gtf[,-c(2,3,4,5,6,7)]
  cai_tpm = merge(cai, gtf, by.x = "gene_id", by.y = "Gene.ID", all = T)
  cai_tpm[cai_tpm == 0] <- NA
  cai_tpm = cai_tpm[complete.cases(cai_tpm),] # Or cai_tpm = na.omit(cai_tpm)
#===========================================================================
# Correlation_CAI_TPM 
#===========================================================================
  svg(file=paste0("/media/hp/disk1/DYY/reference/annotation/", species,"/picture_bymCAI_CC",exp,"/",name,"_CAI_TPM.svg"))
  CAI_cor = cor(cai_tpm$mCAI_value,cai_tpm$TPM)
  CAI_cor
  p <- ggplot(cai_tpm,aes(x = cai_tpm$mCAI_value ,y = cai_tpm$TPM))+
    geom_point(shape = 16,size = 0.75)+
    labs(title = paste0(name,"cor_mCAI_TPM    ","r=",CAI_cor))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl")+
    stat_smooth(method="lm", se=FALSE,linetype="dashed", color = "red",size = 0.75)+
    theme_bw()+
    xlab("mCAI value")+
    ylab("RNAseq TPM")
  print(p)
  dev.off()
#===============================================================================
# Write out data,convenient to calculate the average
#===============================================================================   
  write.table(CAI_cor,file =  paste0
              ("/media/hp/disk1/DYY/reference/annotation/", species,
                "/correlation_bymCAI_CC",exp,"/","CAI_cor.txt"),append = T,quote = FALSE,
              row.names = F, col.names = F)
  }

setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_bymCAI_CC",exp,"/" ))
a = read.table("CAI_cor.txt",sep = '\t',header = F)

aveCAI = mean(a$V1)
aveCAI

write.table(aveCAI,file = "mean_CAI.txt",
            quote = FALSE,row.names = "mean_correlation", 
            col.names = "\tCAI_TPM")

