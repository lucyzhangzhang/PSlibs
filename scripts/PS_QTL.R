# El-Soda et al. 2019
# QTLs mapped to Arabidopsis that respond to Pi supply


wd <- "~/R/Eutrema/PS/QTL"
setwd(wd)

library(xlsx)
library(openxlsx)

# leaf phosphate concentration
AraPHO <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "PHO")

# phosphate use efficiency
AraPUE <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "PUE")

# leaf sulfate concentration
AraSUL <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "SUL")

# leaf sulfate concentration
AraIP6 <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "IP6")

intersect(PHO, PUE)

EuAnnot <- read.csv("../FPKMS_LengthAdjusted.csv", header = T, row.names = 1)
colnames(EuAnnot)[6:17] <- c(paste0("ps", c(1, 2, 3)),
                    paste0("Ps", c(1, 2, 3)),
                    paste0("pS", c(1, 2, 3)),
                    paste0("PS", c(1, 2, 3)))

# run python script
PhytoAPI <- function(name) {
    genes <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = "PHO")
    names <- genes$Gene[!is.na(genes$Gene)]
    names
    write(unique(names), file = paste0(name,".names", sep = "\n"))
    scriptLocale <- "../scripts/phytozomeQTL.py"
    command <- paste("python3", scriptLocale, "-l", paste0(name, ".names"), "-n", "'A. thaliana'")
    system(command)
    results <- list()
    results$matches <- read.table("matches.tab")
    blastCommand <- paste("blastn", "-query", "unmatched.fa", 
                          "-db ~/scratch/misc/2transcripts.fa", "-out blast.out",
                          "-max_target_seqs 1 -outfmt 6")
    system(blastCommand)
    results$blast <- read.table("blast.out")
    results$unmatched <- setdiff(names, union(results$matches[,1], results$blast[,1]))
    return(results)
}
    
PHO.res <- PhytoAPI("PHO")
PHO.res
PUE.res <- PhytoAPI("PUE")
PUE.res
SUL.res <- PhytoAPI("SUL")
SUL.res
IP6.res <- PhytoAPI("IP6")
IP6.res
