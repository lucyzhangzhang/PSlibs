library(DESeq2) #DEG analysis
library(dplyr) #graphing
library(GenomicFeatures) #GFF functions makeTxdbFromGff
library(reshape2) #graphing
# library(VennDiagram) #used for venn
library(venn)

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

#####################FPKM#################
GTF2 <- makeTxDbFromGFF("~/R/Eutrema/PS/drought.gtf", format = "gtf")

GTF2 <- exonsBy(GTF2, by = "gene")
head(GTF2)
lengths <- sum(width(reduce(GTF2)))
lengths <- as.data.frame(lengths)
lengths
inter2 <- intersect(rownames(lengths), rownames(dds))
length(inter2)
lens <- lengths[row.names(lengths) %in% inter2,]
mcols(dds)$basepairs <- lens
#big list of fpkms
FPKM <- fpkm(dds, robust = T)
sampNames <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
colnames(FPKM) <- sampNames
#writing FPKMs to file of all fpkms
write.table(data.frame("Gene" = rownames(FPKM), FPKM), file = "FPKMS.tab", quote = F, row.names = F, sep = "\t")

######################FPKM median length#################
lenMedian <- median(lens)
ddsLenAdj <- dds
mcols(ddsLenAdj)$basepairs <- lenMedian
fpkmLenAdj <- fpkm(ddsLenAdj, robust = T)
head(fpkmLenAdj)
fpkmLenAdj[row.names(fpkmLenAdj) %in% "nXLOC_008023", ]
colnames(fpkmLenAdj) <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
write.table(data.frame("Gene" = rownames(fpkmLenAdj), fpkmLenAdj), file = "FPKMS_LengthAdjusted2.tab", quote = F, row.names = F, sep = "\t")

# data breakdown
ps <- FPKM[,1:3]
Ps <- FPKM[,4:6]
pS <- FPKM[,7:9]
PS <- FPKM[,10:12]

# get the averages of samples treatment-wise
pS <- (pS=rowMeans(pS))
Ps <- (Ps=rowMeans(Ps))
ps <- (ps=rowMeans(ps))
PS <- (PS=rowMeans(PS))


ps.mean <- data.frame(ps, pS, Ps, PS)
ps.filt <- ps.mean
ps.filt <- ps.filt[(ps.mean$ps >= 1) | (ps.mean$Ps >= 1) | (ps.mean$pS >= 1) | (ps.mean$PS >= 1),]
colnames(ps.filt) <- colnames(ps.mean)
for (i in 1:nrow(ps.mean)) {
    row <- ps.mean[i,]
    rownames(row) <- rownames(ps.mean[i,]) 
    if (any(row >= 1)) {
        ps.filt <- rbind(ps.filt, row)
    }        
}

sets <- list(ps=rownames(ps.filt[which(ps.filt$ps != 0),]),
             Ps=rownames(ps.filt[which(ps.filt$Ps != 0),]),
             pS=rownames(ps.filt[which(ps.filt$pS != 0),]),
             PS=rownames(ps.filt[which(ps.filt$PS != 0),]))
pdf("PSIntersect.pdf", width = 5, height = 5)
venn(sets, zcolor = "style")
dev.off()

ps <- rownames(ps.filt[which(ps.filt$ps != 0),])
Ps <- rownames(ps.filt[which(ps.filt$Ps != 0),])
pS <- rownames(ps.filt[which(ps.filt$pS != 0),])
PS <- rownames(ps.filt[which(ps.filt$PS != 0),])
