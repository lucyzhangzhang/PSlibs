# El-Soda et al. 2019
# QTLs mapped to Arabidopsis that respond to Pi supply

wd <- "~/R/Eutrema/PS/QTL"
setwd(wd)

library(xlsx)
library(openxlsx)

AraPHO <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "PHO")

AraPUE <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "PUE")

PHO <- AraPHO$Gene[!is.na(AraPHO$Gene)]
PUE <- AraPUE$Gene[!is.na(AraPUE$Gene)]

intersect(PHO, PUE)

homologs <- function(listgenes, name = "ListGenes") {
    write(unique(listgenes), file = paste0(name, ".names", sep = "\n"))
    return(blastres)
}
write(unique(PHO), file = "PHO.ara.names", sep = "\n")

EuAnnot <- read.csv("../FPKMS_LengthAdjusted.csv", header = T, row.names = 1)
colnames(EuAnnot)[6:17] <- c(paste0("ps", c(1, 2, 3)),
                    paste0("Ps", c(1, 2, 3)),
                    paste0("pS", c(1, 2, 3)),
                    paste0("PS", c(1, 2, 3)))

