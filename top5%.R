### Take the top 5 percent TPM value
species = "Apis_mellifera"
setwd(paste0("/media/hp/Katniss/DYY/aligned/",species))
list.files(getwd())
tpm_array = list.files(getwd(),pattern = "\\d.+[txt]$")
tpm_array
tpm1 = read.table("SRR6833955.T.txt",sep = "\t",header =FALSE)

name = sub(".T.txt", "", "SRR6833955.T.txt")

names(tpm1) = c("transcription_id","gene_id", paste(name, c("FPKM", "TPM"), sep = "_"))
tpm1 = tpm1[,c(-2, -3)]
tpm1

q = quantile(tpm1$SRR6833955_TPM, prob = seq(0.05,1,0.05))
q
df = tpm1[tpm1$SRR6833955_TPM >= q[19],]
df
df2 = tpm1[tpm1$SRR6833955_TPM >= q[18] & tpm1$SRR6833955_TPM < q[19] ,]
df2
df_fi = tpm1[tpm1$SRR6833955_TPM < q[2] ,]
df_fi
write.table(df,file = "top.txt",sep = "\t",quote = FALSE,row.names = FALSE)

###  Correspondence to the top 5 percent TPM value and CAI
tpm1 = read.table("SRR6833955.T.txt",sep = "\t",header =FALSE)

name = sub(".T.txt", "", "SRR6833955.T.txt")

names(tpm1) = c("transcription_id","gene_id", paste(name, c("FPKM", "TPM"), sep = "_"))

tpm1 = tpm1[,c(-2, -3)]
tpm1
cbi  = read.table(paste("CBI_CAI.txt", sep = ""), sep = "\t", header = T)
cor = merge(tpm1, cbi, by.x = "transcription_id", 
             by.y = "transcription_id", all = T)
cor = na.omit(cor)
plot(cor$CAI, cor$SRR6833955_FPKM,log = "y" , 
     pch = 16, col = "cyan",cex = 0.5, main = "CAI_FPKM",
      xlab = "CAI",ylab = "FPKM",xlim = c(0,1))

