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

# removing rows with under 1 FPKM expression in any treatment 
ps.filt <- ps.mean

ps.filt <- ps.filt[(ps.mean$ps >= 1) | (ps.mean$Ps >= 1) | (ps.mean$pS >= 1) | (ps.mean$PS >= 1),]
colnames(ps.filt) <- colnames(ps.mean)


# for (i in 1:nrow(ps.mean)) {
#     row <- ps.mean[i,]
#     rownames(row) <- rownames(ps.mean[i,]) 
#     if (any(row >= 1)) {
#         ps.filt <- rbind(ps.filt, row)
#     }        
# }

sets <- list(ps=rownames(ps.filt[which(ps.filt$ps != 0),]), # 
             Ps=rownames(ps.filt[which(ps.filt$Ps != 0),]),
             pS=rownames(ps.filt[which(ps.filt$pS != 0),]),
             PS=rownames(ps.filt[which(ps.filt$PS != 0),]))
pdf("pics/PSIntersect.pdf", width = 5, height = 5) # 
venn(sets, zcolor="style", cexil = 1, cexsn = 1.2)
zeroset <- matrix(1000*c(0,1,1,0,0,0,0,1,1,0), ncol = 2)
lines(zeroset, col='white', lwd=5)
dev.off()

psn <- rownames(ps.filt[which(ps.filt$ps != 0),])
Psn <- rownames(ps.filt[which(ps.filt$Ps != 0),])
pSn <- rownames(ps.filt[which(ps.filt$pS != 0),])
PSn <- rownames(ps.filt[which(ps.filt$PS != 0),])
# refer to venn diagram sets
n1100 <- setdiff(intersect(Psn, psn), union(pSn, PSn))
n1110 <- setdiff(intersect(intersect(psn, Psn), pSn), PSn)
n0110 <- setdiff(intersect(Psn, pSn), union(psn, PSn))
n0011 <- setdiff(intersect(pSn, PSn), union(Psn, psn))
n0111 <- setdiff(intersect(intersect(PSn, Psn), pSn), psn)
n1010 <- setdiff(intersect(psn, pSn), union(Psn, PSn))
n1011 <- setdiff(intersect(intersect(psn, PSn), pSn), Psn)
n1001 <- setdiff(intersect(psn, PSn), union(Psn, pSn))
n1101 <- setdiff(intersect(intersect(psn, Psn), PSn), pSn)
n0101 <- setdiff(intersect(Psn, PSn), union(pSn, psn))
n0001 <- setdiff(PSn, union(union(psn, Psn), pSn))
n1100 <- data.frame(Genes=n1100, set=rep("n1100", length(n1100)))
n1110 <- data.frame(Genes=n1110, set=rep("n1110", length(n1110)))
n0110 <- data.frame(Genes=n0110, set=rep("n0110", length(n0110)))
n0011 <- data.frame(Genes=n0011, set=rep("n0011", length(n0011)))
n0111 <- data.frame(Genes=n0111, set=rep("n0111", length(n0111)))
n1010 <- data.frame(Genes=n1010, set=rep("n1010", length(n1010)))
n1011 <- data.frame(Genes=n1011, set=rep("n1011", length(n1011)))
n1001 <- data.frame(Genes=n1001, set=rep("n1001", length(n1001)))
n1101 <- data.frame(Genes=n1101, set=rep("n1101", length(n1101)))
n0101 <- data.frame(Genes=n0101, set=rep("n0101", length(n0101)))
n0001 <- data.frame(Genes=n0001, set=rep("n0001", length(n0001)))

setTest <- function(x) {
    setID <- rep(0, 4)
    if (x %in% psn) {
        setID[1] <- 1
    }
    if (x %in% Psn) {
        setID[2] <- 1
    }
    if (x %in% pSn) {
        setID[3] <- 1
    }
    if (x %in% PSn) {
        setID[4] <- 1
    }
   return(paste(setID, collapse = ""))
}

library(ggplot2)

allGenes <- data.frame(Genes = union(union(psn, PSn), union(pSn, Psn)))
allGenes <- mutate(allGenes, setID = sapply(allGenes$Genes, setTest))
someGenes <- allGenes %>% filter(setID != "1111") %>% arrange(., setID)

FPKMsomeGenes <- ps.mean[rownames(ps.mean) %in% someGenes$Genes,]
FPKMsomeGenes <- data.frame(Names=rownames(FPKMsomeGenes), FPKMsomeGenes, stringsAsFactors = F)

# reorder data based on column of another data set
FPKMsomeGenes <- FPKMsomeGenes[order(match(FPKMsomeGenes$Names, someGenes$Genes)), ]
FPKMsomeGenes <- data.frame(setID=someGenes$setID, FPKMsomeGenes, stringsAsFactors = F)

write(unique(FPKMsomeGenes$Names), "someGenes", sep = "\n")
#############DO CREMA SHENANIGANS HERE####################
#predict some genes
setwd("/home/lucy/R/Eutrema/PS/crema")

some.predictions <- read.csv("final_ensemble_predictions.csv", header = T, row.names = 1)
some.predictions <- some.predictions[order(match(toupper(rownames(some.predictions)), toupper(FPKMsomeGenes$Names))),] 

FPKMsomeGenes <- cbind(FPKMsomeGenes, prediction = factor(some.predictions$prediction, levels = c(0, 1)), stringsAsFactors = F)

FPsomeG.melt <- melt(FPKMsomeGenes, id = c("Names", "setID", "prediction"))
FPsomeG.melt$variable <- substr(FPsomeG.melt$variable,1, 2) 
FPsomeG.melt <- FPsomeG.melt[FPsomeG.melt$value != 0,]

some.plot <- {ggplot(FPsomeG.melt, aes(x = variable, y = value, color = setID, shape = prediction, size = prediction)) + 
    geom_point(position = "jitter")
}
    some.plot
ggsave("someFPKMs.pdf", some.plot, dpi = 300, height = 8, width = 6)

PCA <- prcomp(FPKMsomeGenes[,3:6], center = T, scale = T)
plot(PCA)
biplot(PCA)
summary(PCA)
exprVals<-data.frame(PCA$x)
sampleVals<-data.frame(PCA$rotation)

ggsave("PCAintersect.pdf", 
{ggplot() +
    geom_point(data=exprVals, aes(x = PC1, y = PC2, color = FPKMsomeGenes$setID)) +
    geom_segment(data = sampleVals, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(), show.legend = T) + 
    scale_color_discrete(name = "Grouping")
}, dpi = 300, height = 8, width = 9)

######ADDING THE FIELD LIBRARIES##########
setwd("~/R/Eutrema/PS")

#specifying the reference level
 condition <- c("F", "F", "F", "ps", "ps", rep(c("Ps", "pS", "PS"), each = 3))
 sampleF <- c("Y-F2003-1",
               "Y-F2003-2",
               "cracker",
               # sampleNO <- c("ps12_S3",
               #                "ps12_S3",
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

metadata <- data.frame(sample = sampleF,
                       condition = condition,
                       libraryName = sampleF)

metadata <- mutate(metadata, countFile = paste0(metadata$libraryName, "/QC.geneCounts.formatted.for.DESeq.txt"))
sampleTable <- data.frame(sampleName = metadata$libraryName,
                          fileName = metadata$countFile,
                          condition = metadata$condition,
                          sample = metadata$sample)


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = "QoRTs",
                                          design = ~condition)

dds.f <- DESeq(DESeq2Table)
lens <- lengths[row.names(lengths) %in% inter2,]

######################FPKM median length#################
lenMedian <- median(lens)
mcols(dds.f)$basepairs <- lenMedian
fpkm.f <- fpkm(dds.f, robust = T)
fpkm.f[row.names(fpkm.f) %in% "nXLOC_008023", ]
colnames(fpkm.f) <- c("2003-1", "2003-2", "CC", "ps1", "ps2",  "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")

F2003.1 <- rownames(fpkm.f[which(fpkm.f[,1] != 0),])
F2003.2 <- rownames(fpkm.f[which(fpkm.f[,2] != 0),])
ccn <- rownames(fpkm.f[which(fpkm.f[,3] != 0),])
allGeneN <- union(union(psn, PSn), union(pSn, Psn))

fSets <- list(cab=allGeneN,
              F2003.1=F2003.1,
              F2003.2=F2003.2,
              CC=ccn)
pdf("fieldVenn.pdf", width = 5, height = 5)
venn(fSets)
dev.off()

field.union <- union(union(F2003.1, F2003.2), ccn)
write(field.union, "fieldunion.names", sep = "\n")

field.exp.raw <- read.csv("crema/fieldG/final_ensemble_predictions.csv", header = T)
field.specific <- field.exp.raw[!toupper(field.exp.raw[,1]) %in% toupper(allGeneN), ]
field.spec.pred <- field.specific %>% filter(prediction == 1)
# write.table(field.spec.pred, "test", sep = "\t", quote = F)

# expr_annot doesn't have all the field expressed genes because it's generated from lib
full.annot <- read.csv("~/Eutrema/FPKM/2018-10-15-drought_geneFPKM.csv", header = T)
# remove duplicated rows
full.annot <- full.annot[!duplicated(full.annot[,1]),]
rownames(full.annot) <- full.annot[,1]
fsp.annot <- expr_annot[toupper(rownames(full.annot)) %in% toupper(field.spec.pred[,1]),]
write.table(fsp.annot, "fsp.annot.tab", sep = "\t", quote = F)

####################STRING DB#####################

library(STRINGdb)
species <- get_STRING_species(version="11.0", species_name=NULL)
species[grepl(".*Eutrema.*", species$compact_name),]
stringdb <- STRINGdb$new(version = "10", species=3702)

setnames <- rbind(n1100,n1110,n0110,n0011,n0111,
                  n1010,n1011,n1001,n1101,n0101,n0001)

setnames$At <- c(NA,"AT3G26110","AT1G16950",NA, NA,
                 "AT3G41950","AT4G17030", "AT3G55515",
                 "AT5G50720", "AT3G58270", NA, "AT5G45860",
                 "AT5G38197", NA)

means <- data.frame(Genes=rownames(ps.mean), ps.mean)
meansfilt <- means[which(rownames(means) %in% setnames$Genes),]

setmeans <- merge(meansfilt, setnames, by = "Genes")
setmeans <- setmeans[order(setmeans$set),]

nodes <- stringdb$map(setmeans, "At", removeUnmappedRows = F)
stringdb$plot_network(nodes$STRING_id)

AT3G26110 <- stringdb$mp("AT3G26110")
stringdb$get_neighbors(AT3G26110)

# what about the differentially expressed genes?
###################HTTP REQ####################

####################HEATMAPS###################

library(gplots)
# library(dendextend)

distTest <- t(ps.mean)
distTest <- as.matrix(dist(distTest))
heatmap(distTest)

log2fpkm <- log2(FPKM + 1)
distTest <- t(log2fpkm)
distTest <- as.matrix(dist(distTest))
rowCols <- rep(c("darkred", "forestgreen", "orange", "blue"), each = 3)
pdf("pics/fpkmHeatmap.pdf")
heatmap.2(distTest, tracecol = NULL, colRow = rowCols, colCol = rowCols)
dev.off()

# un-reordered
heatmap.2(distTest, dendrogram = "none", Rowv = NULL, Colv = NULL,  
          tracecol = NULL, colRow = rowCols, colCol = rowCols)
