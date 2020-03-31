# El-Soda et al. 2019
# QTLs mapped to Arabidopsis that respond to Pi supply


wd <- "~/R/Eutrema/PS/QTL"
setwd(wd)

library(xlsx)
library(openxlsx)

EuAnnot <- read.csv("../FPKMS_LengthAdjusted.csv", header = T, row.names = 1)
colnames(EuAnnot)[6:17] <- c(paste0("ps", c(1, 2, 3)),
                    paste0("Ps", c(1, 2, 3)),
                    paste0("pS", c(1, 2, 3)),
                    paste0("PS", c(1, 2, 3)))

# run python script
PhytoAPI <- function(name) {
    # get table of genes
    genes <- read.xlsx("ListOfGenesMapped.xlsx",
                    colNames = T, rowNames = F, sheet = name)

    # get names from table
    names <- genes$Gene[!is.na(genes$Gene)]
    write(unique(names), file = paste0(name,".names"), sep = "\n")

    # run API
    scriptLocale <- "../scripts/phytozomeQTL.py"
    command <- paste("python3", scriptLocale, "-l", paste0(name, ".names"), "-n", "'A. thaliana'")
    system(command)

    # store results, R doesn't support multiple return values, have to
    # put everything in a list in order to call it
    results <- list()
    results$matches <- read.table("matches.tab", sep = "\t")
    colnames(results$matches) <- c("AT", "Gene", "Organism", "relationship") 

    # blast the unmatched entries
    blastCommand <- paste("blastn", "-query", "unmatched.fa", 
                          "-db ~/scratch/misc/2transcripts.fa", "-out blast.out",
                          "-max_target_seqs 1 -outfmt 6")
    system(blastCommand)
    results$blast <- read.table("blast.out")

    # label the columns of outfmt 6
    colnames(results$blast) <- c("AT", "Gene", "pident", "length",
                                 "mismatch", "gapopen", "qstart",
                                 "qend", "sstart", "send",
                                 "evalue", "bitscore")

    # determine the leftovers, the order of x and y in setdiff(x,y) matters
    # i.e. the values of x that are not in y
    results$unmatched <- setdiff(names, union(results$matches[,1], results$blast[,1]))

    return(results)
}
    
# leaf phosphate content
PHO.res <- PhytoAPI("PHO")
PHO.res

# phosphate use efficiency
PUE.res <- PhytoAPI("PUE")
PUE.res

# leaf sulfate content
SUL.res <- PhytoAPI("SUL")
SUL.res

# leaf phytate content
IP6.res <- PhytoAPI("IP6")
IP6.res

library(reshape2)
library(ggplot2)

# get FPKMs, generated from ~/scratch/FPKM/script
fpkm <- read.table("PSPiCombinedFPKM.tab", header = T)
colnames(fpkm)[12:23] <- paste0(rep(c("ps", "Ps", "pS", "PS"), each = 3), c(1, 2, 3))


PHO.dat <- fpkm[which(fpkm$Gene %in% union(PHO.res$matches[,2], PHO.res$blast[,2])),]
PHO.melt <- melt(PHO.dat, id.vars = "Gene")


{PHOplot <- ggplot(PHO.melt, aes( x = variable, y = Gene, fill = log2(value + 1))) +
    geom_tile() +
    labs(x = "Condition", title = "PHO QTL group") + 
    scale_fill_gradient(low = "yellow", high = "red3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
PHOplot}

heatplot <- function(name, res) {
    dat <- fpkm[which(fpkm$Gene %in% union(res$matches[,2], res$blast[,2])),]
    dat.melt <- melt(dat, id.vars = "Gene")

    ggplot(dat.melt, aes( x = variable, y = Gene, fill = log2(value + 1))) +
        geom_tile() +
        labs(x = "Condition", title = paste(name, "QTL group")) + 
        scale_fill_gradient(low = "yellow", high = "red3") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
}

QTLgroups <- c("PHO", "PUE", "SUL", "IP6")
resList <- list(PHO=PHO.res, PUE=PUE.res, SUL=SUL.res, IP6=IP6.res)

plots <- mapply(heatplot, name = QTLgroups, res = resList,
                SIMPLIFY = F, USE.NAMES = F)

library(grid)
library(cowplot)
library(gridExtra)

# adding a legend grob by itself
# https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots

# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}

# mylegend <- g_legend(PHOplot)

plotted <- plot_grid(plotlist = plots, ncol = 2, align = "v")


pdf("QTLheatmaps.pdf", width = 14, height = 10)
plotted <- grid.arrange(arrangeGrob(plotted, left = textGrob("Gene", rot = 90),
                                    bottom = textGrob("Condition")))
plotted
dev.off()

binder <- function(name, res, Ara, App = F, excelName = "output.xslx") {
    links <- rbind(res$matches[,1:2], res$blast[,1:2])
    resSubset <- Ara[match(links[,1], Ara$Gene),]
    resSubset <- resSubset[,8]
    bound <- cbind(links, resSubset)
    colnames(bound) <- c("Arabidopsis", "Eutrema", "Description")
    write.xlsx2(bound, excelName, row.names = F, sheetName = name, append = App)
}

AraList <- list(AraPHO, AraPUE, AraSUL, AraIP6)

for (i in 1:length(QTLgroups)) {
    if (i == 1) {
        App <- F
    } else {
        App <- T
    }
    binder(QTLgroups[i], resList[[i]], AraList[[i]], App = App, excelName = "QTLs.xlsx")
}

# unmatched Arabidopsis entries

for (i in 1:length(QTLgroups)) {
    if (i == 1) {
        App <- F
    } else {
        App <- T
    }
    unmatched <- resList[[i]]$unmatched
    Ara <- AraList[[i]]
    Out <- Ara[match(unmatched, Ara$Gene),]
    write.xlsx2(Out, file = "UnmatchedQTL.xlsx", 
                row.names = F, sheetName = QTLgroups[i], append = App)
}
