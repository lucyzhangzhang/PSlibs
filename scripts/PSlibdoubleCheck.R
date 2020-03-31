library(DESeq2)                        #DEG analysis
library(gplots)                        #graphing and visualization
library(RColorBrewer)                  #graphing and visualization
library(genefilter)                    #GFF maker
library(dplyr)                         #graphing
library(geneplotter)                   #graphing
library(ggplot2)                       #graphing
library(GenomicFeatures)               #GFF functions makeTxdbFromGff
library(reshape2)                      #graphing
library(apeglm)                        #better than ashr
# library(venn)                          #used for venn
# library(topGO)                         #used for GO enrichment (results not informative)
# library(Rgraphviz)                     #used for GO enrichment
# library(pheatmap)                      #used for visualization of DESeq outputs

setwd("~/R/Eutrema/PS") #stuff in here

#####################SETUP###################
# condition <- c("F", "F", "F", rep(c("ps", "Ps", "pS", "PS"), each = 3))
condition <- c(rep(c("ps", "Ps", "pS", "PS"), each = 3))
# condition <- c("F", "F", "F", "ps", "ps", rep(c("Ps", "pS", "PS"), each = 3))
# sampleNO <- c("Y-F2003-1",
#               "Y-F2003-2",
#               "cracker",
sampleNO <- c("ps12_S3",
              #               "ps12_S3",
              "ps48_S10",
              "ps49_S11",
              "ps11_S2",
              "ps40_S6",
              "ps46_S9",
              "ps10_S1",
              "ps38_S5",
              "ps44_S8",
              "ps37_S4",
              "ps41_S7",
              "ps50_S12")

# condition <- sampleN
metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 4)))

metadata <- mutate(metadata, countFile = paste0(metadata$libraryName, "/QC.geneCounts.formatted.for.DESeq.txt"))
sampleTable <- data.frame(sampleName = metadata$libraryName,
                          fileName = metadata$countFile,
                          condition = metadata$condition,
                          sampleNO = metadata$sampleNO,
                          reps = metadata$reps)


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = "QoRTs",
                                          design = ~ condition)

#specifying the reference level
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]
#correct for batch effect
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2, 3))
# DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "PS")

dds <- DESeq(DESeq2Table)
vsd <- vst(dds, blind = F)
plotPCA(vsd, intgroup = c("condition"))
rld <- rlogTransformation(dds, blind = F)
{plotPCA(rld, intgroup=c("condition")) +
    geom_label(aes(label = sampleNO))}


GTF2 <- makeTxDbFromGFF("~/R/Eutrema/PS/drought.gtf", format = "gtf")

GTF2 <- exonsBy(GTF2, by = "gene")
head(GTF2)
lengths <- sum(width(reduce(GTF2)))
lengths <- as.data.frame(lengths)
inter2 <- intersect(rownames(lengths), rownames(dds))
length(inter2)
lens <- lengths[row.names(lengths) %in% inter2,]
mcols(dds)$basepairs <- lens
#big list of fpkms
fpkms <- fpkm(dds, robust = T)
sampNames <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
colnames(fpkms) <- sampNames

# test <- prcomp(fpkms, center = T, scale = T)
test <- prcomp(log2(fpkms + 1), center = T, scale = T)
test$rotation
plot(test$rotation[,1:2])
PCs <- data.frame(test$rotation)
PCs$ID <- rownames(PCs)
PCs$Cond <- condition
PCs$sampleNO <- sampleNO
PCs %>% ggplot(aes(x = PC2, y = PC3, color = condition)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1))
##########################################################################
library(rsem)
# Caitlin's test data (w/o the DROUGHT transcripts)
Csamples <- paste0(rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3), c(1, 2, 3))
Cfpkm <- data.frame()

for (i in Csamples) {
    ifpkm <- read.table(paste0("~/R/Eutrema/PS/caitlinRSEM/", i, ".genes.results"), header = T)
    ifpkm <- ifpkm[,c(1,7)]
    colnames(ifpkm) <- c("gene_id", i)
    if (length(Cfpkm) == 0) {
        Cfpkm <- ifpkm
    } else {
        Cfpkm <- merge(Cfpkm, ifpkm, by = "gene_id")
    }
}

rownames(Cfpkm) <- Cfpkm[,1]
Cfpkm <- Cfpkm[,2:ncol(Cfpkm)]
colnames(Cfpkm) <- sampNames
# Ctest <- prcomp(Cfpkm, center = T, scale = T)
Ctest <- prcomp(log2(Cfpkm + 1), center = T, scale = T)
Ctest$rotation
plot(Ctest$rotation[,1:2])
CPCs <- data.frame(Ctest$rotation)
CPCs$ID <- rownames(CPCs)
CPCs$Cond <- condition

library(cowplot)
Caitlin <- CPCs %>% ggplot(aes(x = PC2, y = PC4, color = condition)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1))  + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + geom_label(aes(label = sampNames))
Me <- PCs %>% ggplot(aes(x = PC2, y = PC4, color = condition)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed")+ geom_label(aes(label = sampleNO))

 ggsave("pics/CaitlinCompare.pdf", plot_grid(Caitlin, Me, nrow = 1), dpi = 250)
