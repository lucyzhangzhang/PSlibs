library(DESeq2) #DEG analysis
library(dplyr) #graphing
library(GenomicFeatures) #GFF functions makeTxdbFromGff
library(reshape2) #graphing
library(VennDiagram) #used for venn

setwd("~/R/Eutrema/PS") #stuff in here

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

# don't filter for minimum number of reads
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = "QoRTs",
                                          design = ~ reps + condition)

dds <- DESeq(DESeq2Table)
