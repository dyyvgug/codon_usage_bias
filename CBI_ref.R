#!/usr/bin/env Rscript
setwd("/media/hp/disk1/DYY/reference/annotation/")
## take the arguments from environment
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

species = args[1]
species = "Arabidopsis_thaliana"
df = read.table(paste("/media/hp/disk1/DYY/reference/annotation/",species,
                      "/ref/codon_frequency.txt",sep = ""), sep = "\t", header = F )

a = aggregate(V3~V2, df, max)

df = merge(df, a , by.x = "V2", by.y = "V2")

df$CAI = round(df$V3.x/df$V3.y, 3)

category = function(x)   {
  if (x > 0.8)  {
    y = 3
  }
  else if (x > 0.3 & x <= 0.8)  {
    y = 2
  }
  else {
    y = 1
  }
  return(y)
}

df$cat = sapply(1:64, function(x) category(df[x,"CAI"]))
df = df[,c(-3, -4)]
names(df) = c("AA", "codon", "CAI", "rand")

write.table(df, paste("/media/hp/disk1/DYY/reference/annotation/",species,
                      "/ref/CBI_ref.txt",sep =''), sep = "\t", row.names = F, quote = F)