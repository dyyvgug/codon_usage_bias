setwd("~/Documents/IUPred2")
list.files(getwd())
t_array = list.files(getwd(),pattern = ".+[csv]$")
t_array

df = read.table("HnRNPA1.fa.CAI.csv", sep = "\t", header = F )
