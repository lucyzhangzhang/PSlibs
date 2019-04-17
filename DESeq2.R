library(DESeq2)
library(plyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
#library(biomaRt)
library(dplyr)
#library(EDASeq) #doesn't work
library(fdrtool)
library(geneplotter)
library(ggplot2)
library(GenomicFeatures)
library(reshape2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(Rgraphviz)

setwd("~/R/Eutrema/PS")

condition <- rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3)
sampleNO <- c("ps10_S1",
             "ps11_S2",
             "ps12_S3",
             "ps37_S4",
             "ps38_S5",
             "ps40_S6",
             "ps41_S7",
             "ps44_S8",
             "ps46_S9",
             "ps48_S10",
             "ps49_S11",
             "ps50_S12")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO)

metadata <- mutate(metadata, countFile = paste0(metadata$libraryName, "/QC.geneCounts.formatted.for.DESeq.txt"))
sampleTable <- data.frame(sampleName = metadata$libraryName,
                fileName = metadata$countFile,
                condition = metadata$condition,
                sampleNO = metadata$sampleNO)


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = "QoRTs",
                                      design = ~ condition)

as.data.frame( colData(DESeq2Table) )
# dds <- DESeq(DESeq2Table)
# res <- results(dds)
# res
# sum( res$pvalue < 0.05, na.rm=TRUE )
# sum( res$padj < 0.05, na.rm=TRUE )


#see how many genes have non-zero expression
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]

GeneCounts <- counts(DESeq2Table)

idx.nz <- apply(GeneCounts, 1, function(x)
{
  all(x > 0)
})
sum(idx.nz)

DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

#Empirical Cumulative Distribution Function
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1500))

#densities
multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))


rld <- rlogTransformation(DESeq2Table, blind=TRUE)

#make heatmap plots
png("HeatmapPlots.png")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

png("PCA_samplesno0.png", height = 4, width = 6, units = "in", res = 300)
DESeq2::plotPCA(rld, intgroup=c("condition"))
dev.off()

#remove the low read count sample...but it wasn't that one causing the problems
#DESeq2::plotPCA(rlogTransformation(DESeq2Table[,-3], blind = T))

#estimate dispersions
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <- nbinomWaldTest(DESeq2Table)

#eeeeh why are none of the ajusted p-values are less than 0.05 whaat
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")


table(DESeq2Res$padj < 0.05)
png("plot_MA.png", height = 4, width = 6, res = 300, units = "in")
plotMA(DESeq2Res, ylim = c(-4, 4))
dev.off()
#attr(DESeq2Res,"filterThreshold")
metadata(DESeq2Res)$filterThreshold



plot(metadata(DESeq2Res)$filterNumRej,type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

hist(DESeq2Res$pvalue, col = "lavender",
     main = "pvals", xlab = "p-values")

# #remove NA pval
# DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
# #remove NA padj
# DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
# #remove padj, add fdr results later on
# DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
#
# #FDR estimates
# FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
#
# #null model variance
# FDR.DESeq2Res$param[1, "sd"]
#
# #add new padj
# DESeq2Res[,"padj"] <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
#
#
# hist(FDR.DESeq2Res$pval, col = "royalblue4",
#      main = "Padj", xlab = "CORRECTED p-values")
#Try to get FPKMs

diff <- which(DESeq2Res[,"padj"] < 0.05)
diffName <- DESeq2Res[diff,]
diffName <- rownames(diffName)
diffName
myGTF <- read.table("~/scratch/v1.0/drought.gtf")
myGTF <- myGTF[,-c(9, 11, 12, 14)]
colnames(myGTF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "geneID", "transcriptID")
myGTF <- mutate(myGTF, length = abs(myGTF$start - myGTF$end))
myGTF <- myGTF %>% filter( grepl("transcript", feature)) %>% arrange(geneID)
head(myGTF)

GTF2 <- makeTxDbFromGFF("~/scratch/v1.0/drought.gtf", format = "gtf")
GTF2 <- exonsBy(GTF2, by = "gene")
head(GTF2)
lengths <- sum(width(reduce(GTF2)))
lengths <- as.data.frame(lengths)
lengths

#GTF
# inter <- intersect(myGTF$geneID, rownames(DESeq2Table))
# head(inter)
# smolGTF <- myGTF %>% filter( geneID %in% inter )

#GFF
inter2 <- intersect(rownames(lengths), rownames(DESeq2Table))
length(inter2)
lens <- lengths[row.names(lengths) %in% inter2,]
mcols(DESeq2Table)$basepairs <- lens

#mcols(DESeq2Table)$basepairs <- myGTF[,ncol(myGTF)]
fpkms <- fpkm(DESeq2Table, robust = T)
diffDP <- fpkms[row.names(fpkms) %in% diffName, ]
diffDP <- as.data.frame(diffDP)
write.table(data.frame("Gene" = rownames(diffDP), diffDP), file = "diffDP.tab", quote = F, row.names = F)

graphs <- list()
for (i in 1:nrow(diffDP)) {
  local({
data <- cbind(melt(diffDP[i,]), condition)
graph <- ggplot(data, aes(x = variable, y = value, fill = factor(condition), color = factor(condition))) +
  geom_bar(stat = "identity", position = "dodge") +
 # ylim(0, 520) +
  labs(x = paste(rownames(diffDP)[i], "expression", sep = " ")) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank())
  graphs[[i]] <<- graph
  })
}


oneg <- ggplot(cbind(melt(diffDP[1,]), condition), aes(x = variable, y = value, fill = factor(condition), color = factor(condition))) +
  geom_bar(stat = "identity", position = "dodge") +
  # ylim(0, 520) +
  labs(x = paste(rownames(diffDP)[i], "expression", sep = " ")) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank())
oneg
# graphs

grid <- plot_grid(plotlist = graphs, ncol = 5)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


legend <- g_legend(oneg)
legend

pluselg <- plot_grid(grid, legend, ncol = 2, rel_widths = c(1, 0.05))
finalPlot <- grid.arrange(arrangeGrob(pluselg, left = textGrob("FPKM", rot = 90, vjust = 1)))
# ggdraw(add_sub(pluselg, "FPKM", angle = 90, fontfamily = "serif", size = 10, hjust = -15, vjust = -85))

ggsave("Unmerged-Analysis.png", finalPlot, height = 8, width = 12, dpi = 250)




#########################################################################################################3
DESeq2Table2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = "QoRTs/smolMerged2",
                                          design = ~ condition)

as.data.frame( colData(DESeq2Table2) )

keep <- rowSums(counts(DESeq2Table2)) >= 10
DESeq2Table2 <- DESeq2Table2[keep,]
GeneCounts <- counts(DESeq2Table2)

DESeq2Table2 <- estimateSizeFactors(DESeq2Table2)
sizeFactors(DESeq2Table2)

#Empirical Cumulative Distribution Function
multiecdf( counts(DESeq2Table2, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1500))

#densities
multidensity( counts(DESeq2Table2, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))


rld2 <- rlogTransformation(DESeq2Table2, blind=TRUE)

#make heatmap plots
png("HeatmapPlots.png")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

png("PCA_samplessmol.png", height = 4, width = 6, units = "in", res = 300)
DESeq2::plotPCA(rld2, intgroup=c("condition"))
dev.off()

#remove the low read count sample...but it wasn't that one causing the problems
#DESeq2::plotPCA(rlogTransformation(DESeq2Table[,-3], blind = T))

#estimate dispersions
DESeq2Table2 <- estimateDispersions(DESeq2Table2)
plotDispEsts(DESeq2Table2)

DESeq2Table2 <- nbinomWaldTest(DESeq2Table2)

#eeeeh why are none of the ajusted p-values are less than 0.05 whaat
DESeq2Res2 <- results(DESeq2Table2, pAdjustMethod = "BH")


table(DESeq2Res$padj < 0.05)
png("plot_MA.png", height = 4, width = 6, res = 300, units = "in")
plotMA(DESeq2Res2, ylim = c(-4, 4))
dev.off()
#attr(DESeq2Res,"filterThreshold")
metadata(DESeq2Res2)$filterThreshold



plot(metadata(DESeq2Res2)$filterNumRej,type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

hist(DESeq2Res2$pvalue, col = "lavender",
     main = "pvals", xlab = "p-values")

diff2 <- which(DESeq2Res2[,"padj"] < 0.05)
diffName2 <- DESeq2Res2[diff2,]
diffName2 <- rownames(diffName2)
diffName2
myGTF <- read.table("~/scratch/v1.0/drought.gtf")
myGTF <- myGTF[,-c(9, 11, 12, 14)]
colnames(myGTF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "geneID", "transcriptID")
myGTF <- mutate(myGTF, length = abs(myGTF$start - myGTF$end))
myGTF <- myGTF %>% filter( grepl("transcript", feature)) %>% arrange(geneID)
head(myGTF)

GTF2 <- makeTxDbFromGFF("~/scratch/v1.0/drought.gtf", format = "gtf")
GTF2 <- exonsBy(GTF2, by = "gene")
head(GTF2)
lengths <- sum(width(reduce(GTF2)))
lengths <- as.data.frame(lengths)
lengths

#GTF
# inter <- intersect(myGTF$geneID, rownames(DESeq2Table))
# head(inter)
# smolGTF <- myGTF %>% filter( geneID %in% inter )

#GFF
inter2 <- intersect(rownames(lengths), rownames(DESeq2Table2))
length(inter2)
lens <- lengths[row.names(lengths) %in% inter2,]
mcols(DESeq2Table2)$basepairs <- lens

#mcols(DESeq2Table)$basepairs <- myGTF[,ncol(myGTF)]
fpkms2 <- fpkm(DESeq2Table2, robust = T)
diffDP2 <- fpkms2[row.names(fpkms2) %in% diffName2, ]
diffDP2 <- as.data.frame(diffDP2)
write.table(data.frame("Gene" = rownames(diffDP2), diffDP2), file = "diffDP2.tab", quote = F, row.names = F)

graphs <- list()
for (i in 1:nrow(diffDP2)) {
  local({
data <- cbind(melt(diffDP2[i,]), condition)
graph <- ggplot(data, aes(x = variable, y = value, fill = factor(condition), color = factor(condition))) +
  geom_bar(stat = "identity", position = "dodge") +
 # ylim(0, 520) +
  labs(x = paste(rownames(diffDP2)[i], "expression", sep = " ")) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank())
  graphs[[i]] <<- graph
  })
}
oneg2 <- ggplot(cbind(melt(diffDP2[1,]), condition), aes(x = variable, y = value, fill = factor(condition), color = factor(condition))) +
  geom_bar(stat = "identity", position = "dodge") +
  # ylim(0, 520) +
  labs(x = paste(rownames(diffDP2)[1], "expression", sep = " ")) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank())
oneg2
# graphs

grid <- plot_grid(plotlist = graphs, ncol = 3)


legend <- g_legend(oneg2)
legend

pluselg <- plot_grid(grid, legend, ncol = 2, rel_widths = c(1, 0.05))
finalPlot <- grid.arrange(arrangeGrob(pluselg, left = textGrob("FPKM", rot = 90, vjust = 1)))
# ggdraw(add_sub(pluselg, "FPKM", angle = 90, fontfamily = "serif", size = 10, hjust = -15, vjust = -85))
ggsave("Merged-Analysis.png", finalPlot, height = 8, width = 12, dpi = 250)

IPS <- fpkms2[row.names(fpkms2) %in% c("XLOC_008023"), ]
IPS_melt <- melt(IPS)

IPS_plot <- ggplot(IPS_melt,
                        aes(x = row.names(IPS_melt), y = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Treatment", y = "FPKM") +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_fill_discrete(breaks = c("Thhalv10015137m.g"),
  #                     labels = c("IPS2", "IPS1")) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

print(IPS_plot)

genes <- read.table("Genes")

graph <- function(gene) {
  dat <- fpkms2[row.names(fpkms2) %in% gene,]
  melt <- melt(dat)
  plot <- ggplot(melt, aes( x = row.names(melt), y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = gene) +
    theme(text = element_text(size=10,
                              family="serif"),
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "top",
          legend.title = element_blank())
}

thing <- sapply(genes[,2], graph, USE.NAMES = F, simplify = F)
head(thing)

plot_grid(plotlist = thing, ncol = 4)


