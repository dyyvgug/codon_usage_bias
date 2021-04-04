#!/usr/bin/env Rscript
setwd("/media/hp/disk1/DYY/reference/annotation/")
## take the arguments from environment
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

species = args[1]

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



##### run perl script ############
system("perl xxxx")



## --------------------------------
## merge CAI and FPKM
## --------------------------------
setwd("/media/hp/disk1/DYY/reference/annotation/")
outcomes = function (species, expression = "V4")  {
  cbi  = read.table(paste("/media/hp/disk1/DYY/reference/annotation/", 
                          species, "/ref/CBI_CAI.txt", sep = ""), sep = "\t", header = T)
  
  setwd(paste("/media/hp/Katniss/DYY/aligned/",species,sep = ""))
  gtf_array = list.files(getwd(),pattern = "\\d.+[txt]$") 
  
  df = cbi
  
  for (i in gtf_array) {
    gtf = read.table(i, sep = "\t", header=F)
    name = sub(".txt", "", i)
    names(gtf) = c("transcription_id", "gene_id", paste(name, c("FPKM", "TPM"), sep = "_"))
                   
    gtf = gtf[,c(-2, -3)]
    df = merge(df, gtf, by.x = "transcription_id", 
             by.y = "transcription_id", all = T)
   # df=cbi,than df can cumulative merge result.
     }
    
    #gtf = gtf [ , apply(gtf, 2, function(x) !any(is.na(x)))]  # remove columns have NA
    
  #change all NA to 0
  df[is.na(df)] = 0
  
  
  # find max and min for each row
  df$max = apply(df[,5:ncol(df)], 1, max)
  df$min = apply(df[,5:ncol(df)], 1, min)
  df$change = log2(df$max/df$min)
  df = df[df$max>1,] 
    #df = df[df$V3 != 0,]
  
  q = quantile(df$change[df$change> 0], prob = seq(0,1,0.2))
  
  q
   
  for ( j in 1:(length(q)-2) )  {
    
  
  
  df1 = df[df$change >= q[j] & df$change < q[j+1],]
  write.table(df1,paste0("/media/hp/disk1/DYY/reference/TPM_CAI/",species,
              species, " range ", round(q[j],2), "-", round(q[j+1],2),".txt"),
              sep = "\t")
  
  for ( i in 5:(ncol(df)-3) )  {
    
      cor = cor.test(df1$avg_CAI, log(df1[[get('i')]]))
   
    jpeg(paste0("/media/hp/disk1/DYY/reference/CAI_CBI_FT/",species,"/FPKM/", 
                colnames(df1)[i], " range ", round(q[j],2), "-", round(q[j+1],2),".jpg"), 
                width = 480, height = 480, quality = 100)
    #par(new = T)
    plot(df1$avg_CAI[df1[[get('i')]]> 0], df1[[get('i')]][df1[[get('i')]] > 0], log = "y" , 
         pch = 16, col = "cyan",cex = 0.5, main = round(cor$estimate, 2),
         xlab = "CAI",ylab = "TPM",xlim = c(0.6,1))
    
    dev.off()
  }
  }
}



