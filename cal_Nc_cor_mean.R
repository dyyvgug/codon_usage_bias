species = "Apis_mellifera"
setwd(paste0("/media/hp/disk1/DYY/reference/annotation/",species,"/correlation_Nc1/" ))
a = read.table("allcor.txt",sep = '\t',header = F)
names(a) = c("all_Nc_cor","top_Nc_cor","high_Nc_cor","low_Nc_cor","def_Nc_cor")
write.table(paste(mean(a$all_Nc_cor),mean(a$top_Nc_cor),mean(a$high_Nc_cor),
                  mean(a$low_Nc_cor),mean(a$def_Nc_cor),sep = '\t'),file = "mean_cor.txt",
            quote = FALSE,row.names = "mean_Nc_correlation", 
            col.names = "all_Nc_cor\ttop_Nc_cor\thigh_Nc_cor\tlow_Nc_cor\tdef_Nc_cor")
