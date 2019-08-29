library(DESeq2) #DEG analysis
library(plyr) #data manipulation
library(gplots) #graphing and visualization
# library(RColorBrewer) #graphing and visualization
library(genefilter) #GFF maker
library(dplyr) #graphing
library(geneplotter) #graphing
# library(ggplot2) #graphing
library(GenomicFeatures) #GFF functions makeTxdbFromGff
library(reshape2) #graphing
# library(apeglm) #better than ashr
library(venn) #used for venn
library(VennDiagram) #used for venn
# library(topGO)  #used for GO enrichment (results not informative)
# library(Rgraphviz)  #used for GO enrichment
# library(pheatmap) #used for visualization of DESeq outputs

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
                                          design = ~ reps + condition)

#specifying the reference level
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]
#correct for batch effect
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2, 3))
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "PS")

dds <- DESeq(DESeq2Table)
# res <- results(dds, pAdjustMethod = "BH")
# res
# rld <- rlogTransformation(dds, blind = F)
# png("PCA_samples.png", height = 600, width = 1000, res = 125)
# plotPCA.san(rld, intgroup=c("condition"))
# dev.off()
# vsd <- vst(dds, blind = F)
# plotPCA.san(vsd, intgroup = c("condition"))
# DESeq2::plotPCA(vsd, intgroup = c("condition"))
#https://www.biostars.org/p/243695/
library(ggrepel)
# plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
# {
#   rv <- rowVars(assay(object))
#   select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
#                                                      length(rv)))]
#   pca <- prcomp(t(assay(object)[select, ]))
#   percentVar <- pca$sdev^2/sum(pca$sdev^2)
#   if (!all(intgroup %in% names(colData(object)))) {
#     stop("the argument 'intgroup' should specify columns of colData(dds)")
#   }
#   intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
#   group <- if (length(intgroup) > 1) {
#     factor(apply(intgroup.df, 1, paste, collapse = " : "))
#   }
#   else {
#     colData(object)[[intgroup]]
#   }
#   d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group,
#                   intgroup.df, name = colData(rld)[,1])
#   if (returnData) {
#     attr(d, "percentVar") <- percentVar[2:3]
#     return(d)
#   }
#   ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3)
# 
# }


# rldS <- rld[, rld$condition %in% c("ps", "Ps", "pS", "PS")]

#how many genes are significantly expressed
png("PCA_samplesno.png", height = 700, width = 1000, res = 125)
DESeq2::plotPCA(rldS, intgroup=c("condition"))
dev.off()
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_PS"  "condition_Ps_vs_PS" "condition_pS_vs_PS" "condition_ps_vs_PS"

res1_2 <- lfcShrink(dds, coef = "condition_Ps_vs_PS", type = "apeglm")
res1_3 <- lfcShrink(dds, coef = "condition_pS_vs_PS", type = "apeglm")
res1_4 <- lfcShrink(dds, coef = "condition_ps_vs_PS", type = "apeglm")
# res1_I <- lfcShrink(dds, coef = "Intercept", type = "apeglm")
# res1_F <- lfcShrink(dds, coef = "condition_F_vs_PS", type = "apeglm")
res.Ps <- results(dds, contrast = c("condition", "Ps", "PS"))
res.pS <- results(dds, contrast = c("condition", "pS", "PS"))
res.ps <- results(dds, contrast = c("condition", "ps", "PS"))
#########################PLOTTING###############################################
# r1_I <- cbind.data.frame(row.names(res1_I), rep(0, nrow(res1_I)), res1_I$lfcSE, rep("PS", nrow(res1_I)))
# colnames(r1_I) <- c("gene", "lfc", "se", "sample")
# r1_2 <- cbind.data.frame(row.names(res1_2), res1_2$log2FoldChange, res1_2$lfcSE, rep("Ps", nrow(res1_2)))
# colnames(r1_2) <- c("gene", "lfc", "se", "sample")
# r1_3 <- cbind.data.frame(row.names(res1_3), res1_3$log2FoldChange, res1_3$lfcSE, rep("pS", nrow(res1_3)))
# colnames(r1_3) <- c("gene", "lfc", "se", "sample")
# r1_4 <- cbind.data.frame(row.names(res1_4), res1_4$log2FoldChange, res1_4$lfcSE, rep("ps", nrow(res1_4)))
# colnames(r1_4) <- c("gene", "lfc", "se", "sample")
r1_2 <- cbind.data.frame(row.names(res.Ps), res.Ps$log2FoldChange, res.Ps$lfcSE, rep("Ps", nrow(res.Ps)))
colnames(r1_2) <- c("gene", "lfc", "se", "sample")
r1_3 <- cbind.data.frame(row.names(res.pS), res.pS$log2FoldChange, res.pS$lfcSE, rep("pS", nrow(res.pS)))
colnames(r1_3) <- c("gene", "lfc", "se", "sample")
r1_4 <- cbind.data.frame(row.names(res.ps), res.ps$log2FoldChange, res.ps$lfcSE, rep("ps", nrow(res.ps)))
colnames(r1_4) <- c("gene", "lfc", "se", "sample")
# r1_F <- cbind.data.frame(row.names(res1_F), res1_F$log2FoldChange, res1_F$lfcSE, rep("F", nrow(res1_F)))
# colnames(r1_F) <- c("gene", "lfc", "se", "sample")
lists <- rbind(r1_2, r1_3, r1_4)
IPS <- lists[lists$gene %in% "Thhalv10015137m.g", ]

ggplot(IPS, aes(x = sample, y = lfc)) +
  geom_point() +
  xlab(paste0(as.character(IPS$gene[1]), " expression")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(x = sample, ymax = lfc + se, ymin = lfc - se), width = 0.1) +
  ylim(-11,11) +
    theme(text = element_text(size=10,
                              family="serif"),
          axis.text.x = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank())

genes <- read.table("Genes")
genes2 <- genes[genes[,2] %in% lists$gene, ]
oldNames <- genes2[,2]

library(grid) #used for venn and graphs
library(gridExtra) #used for Venn and graphs
library(cowplot)  #used for log2 expression graphs
#very important function which plots everything in a list of gene names
graph <- function(gene) {
  #gets the data from big "melted" list
  dat <- lists[lists$gene %in% gene,]
  #graphs the data
  ggplot(dat, aes(x = sample, y = lfc)) +
    geom_point() +
    labs(x = genes[genes$V2 %in% gene,]$V1) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(x = sample, ymax = lfc + se, ymin = lfc - se), width = 0.1) +
    ylim(-11,11) +
    # ylim(min(dat$lfc-dat$se-1),max(dat$lfc+dat$se+1)) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank())
}

thing <- sapply(oldNames, graph, USE.NAMES = F, simplify = F)
plotted <- plot_grid(plotlist = thing,  ncol = 4, align =  "v")
# plotted
plotted <- grid.arrange(arrangeGrob(plotted, left = textGrob("log2-fold Change", rot = 90, vjust = 1)))
# ggsave("BeforeDotComp.png", plotted, height = 12, width = 11, dpi = 300)
ggsave("DEGnoShrink.png", plotted, height = 12, width = 11, dpi = 300)

################################TREES###################################

png("FourVennvs2.png", width = 1000, height = 1000, res = 200)
venn(list(Ps=rownames(res.Ps[which(res.Ps[,"padj"] < 0.05),]),
pS=rownames(res1_3[which(res1_3[,"padj"] < 0.05),]),
ps=rownames(res1_4[which(res1_4[,"padj"] < 0.05),]),
F=rownames(res1_F[which(res1_F[,"padj"] < 0.05),])), zcolor = "style")
dev.off()


png("smolvenn.png", width = 800, height = 800, res = 250)
venn(list(Ps=rownames(res.Ps[which(res.Ps[,"padj"] < 0.05),]),
                       pS=rownames(res.pS[which(res.pS[,"padj"] < 0.05),]),
                       ps=rownames(res.ps[which(res.ps[,"padj"] < 0.05),])), zcolor = "style")
dev.off()


upPS <- venn.diagram(list(Ps=rownames(res.Ps[which(res.Ps[,"padj"] < 0.05 & res.Ps[,"log2FoldChange"] > 0),]),
                  pS=rownames(res.pS[which(res.pS[,"padj"] < 0.05 & res.pS[,"log2FoldChange"] > 0),]),
                  ps=rownames(res.ps[which(res.ps[,"padj"] < 0.05 & res.ps[,"log2FoldChange"] > 0),])),
                  filename = NULL, col = "red")
downPS <-venn.diagram(list(Ps=rownames(res.Ps[which(res.Ps[,"padj"] < 0.05 & res.Ps[,"log2FoldChange"] < 0),]),
                   pS=rownames(res.pS[which(res.pS[,"padj"] < 0.05 & res.pS[,"log2FoldChange"] < 0),]),
                   ps=rownames(res.ps[which(res.ps[,"padj"] < 0.05 & res.ps[,"log2FoldChange"] < 0),])),
#                   F=rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] < 0),])),
                   filename = NULL, col = "navyblue")

png("vennDiagramJuly.png", height = 800, width = 1200, res = 125)
grid.arrange(gTree(children = gList(textGrob("Upregulated", gp = gpar(fontfamily = "HersheySerif")))),
  gTree(children = gList(textGrob("Downregulated", gp = gpar(fontfamily = "HersheySerif")))),
gTree(children = upPS),
gTree(children = downPS),
ncol = 2, widths = c(1, 1), heights = c(0.04, 1))
dev.off()



#############################GO Enrichment##################################

#read in GO relations
GOs <- readMappings("Drought.go", sep = "\t", IDsep = ";")

#Get results from DESeq

GOres1_2 <- as.data.frame(res1_2@listData)
rownames(GOres1_2) <- res1_2@rownames
GOres1_3 <- as.data.frame(res1_3@listData)
rownames(GOres1_3) <- res1_3@rownames
GOres1_4 <- as.data.frame(res1_4@listData)
rownames(GOres1_4) <- res1_4@rownames
GOres1_F <- as.data.frame(res1_F@listData)
rownames(GOres1_F) <- res1_F@rownames

dat <- list(Ps=GOres1_2, pS=GOres1_3, ps=GOres1_4, F=GOres1_F)
names <- c("Ps", "pS", "ps", "F")

i=1

topGoPipe <- for (data in dat) {
#gene filter to get significant padj genes
  gene_filter <- function(allScore) {
    return(allScore < 0.05)
  }

  name = names[i]

  data <- as.data.frame(dat[[i]], stringAsFactors=F)

  #Up regulated genes
  Fup <- data[data$log2FoldChange > 0,]$padj
  names(Fup) <- rownames(data[data$log2FoldChange > 0, ])
  Fup <- Fup[complete.cases(Fup)]

  pos_genes <- new("topGOdata",
                   description = paste0("PS vs. ", name),
                   ontology = "BP",
                   allGenes = Fup,
                   geneSel = gene_filter,
                   nodeSize = 10,
                   annotationFun = annFUN.gene2GO,
                   gene2GO = GOs)

  pos <- topGO::genes(pos_genes)

  sig <- sigGenes(pos_genes)

  testFisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher")

  resFtest <- getSigGroups(pos_genes, testFisher)

  posTable <- GenTable(pos_genes, Fisher = resFtest, topNodes = 200)
  write.table(posTable[,c(1,ncol(posTable))], file = paste0("smol", name, "_POS.tab"), sep = "\t", quote = F, row.names = F)
  png(paste0(name, "GOPos.png"), height = 1500, width = 1500, res = 200)
  showSigOfNodes(pos_genes, score(resFtest), firstSigNodes = 10, useInfo = "all")
  dev.off()

  #Down regulated genes
  Fdown <- data[data$log2FoldChange < 0,]$padj
  names(Fdown) <- rownames(data[data$log2FoldChange < 0, ])
  Fdown <- Fdown[complete.cases(Fdown)]

  neg_genes <- new("topGOdata",
                   description = paste0("PS vs. ", name),
                   ontology = "BP",
                   allGenes = Fdown,
                   geneSel = gene_filter,
                   nodeSize = 10,
                   annotationFun = annFUN.gene2GO,
                   gene2GO = GOs)

  neg <- topGO::genes(neg_genes)

  sigNeg <- sigGenes(neg_genes)

  testFisherN <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher")

  resFtestN <- getSigGroups(neg_genes, testFisherN)

  negTable <- GenTable(neg_genes, Fisher = resFtestN, topNodes = 200)
  write.table(negTable[,c(1,ncol(negTable))], file = paste0("smol", name, "_NEG.tab"), sep = "\t", quote = F, row.names = F)

  png(paste0(name, "GONeg.png"), height = 1500, width = 1500, res = 200)
  showSigOfNodes(neg_genes, score(resFtestN), firstSigNodes = 10, useInfo = "all")
  dev.off()
  i <- i + 1
}

name <- "ps"
data <- GOres1_4
  Fup <- data[data$log2FoldChange > 0,]$padj
  names(Fup) <- rownames(data[data$log2FoldChange > 0, ])
  Fup <- Fup[complete.cases(Fup)]

  pos_genes <- new("topGOdata",
                   description = paste0("PS vs. ", name),
                   ontology = "BP",
                   allGenes = Fup,
                   geneSel = gene_filter,
                   nodeSize = 10,
                   annotationFun = annFUN.gene2GO,
                   gene2GO = GOs)

  pos <- topGO::genes(pos_genes)

  sig <- sigGenes(pos_genes)

  testFisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher")

  resFtest <- getSigGroups(pos_genes, testFisher)

  posTable <- GenTable(pos_genes, Fisher = resFtest, topNodes = 200)
  write.table(posTable, file = paste0(name, "_POSlong.tab"), sep = "\t", quote = F)
  png(paste0(name, "GOPos.png"), height = 1500, width = 1500, res = 200)
  showSigOfNodes(pos_genes, score(resFtest), firstSigNodes = 10, useInfo = "all")
  dev.off()

###########################GET LIST OF GENE NAMES #################


rownames(res.Ps[which(res.Ps[,"padj"] < 0.05 & res.Ps[,"log2FoldChange"] > 0),])
rownames(res.Ps[which(res.Ps[,"padj"] < 0.05 & res.Ps[,"log2FoldChange"] < 0),])
rownames(res.pS[which(res.pS[,"padj"] < 0.05 & res.pS[,"log2FoldChange"] > 0),])
rownames(res.pS[which(res.pS[,"padj"] < 0.05 & res.pS[,"log2FoldChange"] < 0),])
rownames(res.ps[which(res.ps[,"padj"] < 0.05 & res.ps[,"log2FoldChange"] > 0),])
rownames(res.ps[which(res.ps[,"padj"] < 0.05 & res.ps[,"log2FoldChange"] < 0),])
rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] > 0),])
rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] < 0),])


# rownames(res1_2[which(res1_2[,"padj"] < 0.2 & res1_2[,"log2FoldChange"] > 0),])
# rownames(res1_2[which(res1_2[,"padj"] < 0.2 & res1_2[,"log2FoldChange"] < 0),])
# rownames(res1_3[which(res1_3[,"padj"] < 0.2 & res1_3[,"log2FoldChange"] > 0),])
# rownames(res1_3[which(res1_3[,"padj"] < 0.2 & res1_3[,"log2FoldChange"] < 0),])
# rownames(res1_4[which(res1_4[,"padj"] < 0.2 & res1_4[,"log2FoldChange"] > 0),])
# rownames(res1_4[which(res1_4[,"padj"] < 0.2 & res1_4[,"log2FoldChange"] < 0),])
# rownames(res1_F[which(res1_F[,"padj"] < 0.2 & res1_F[,"log2FoldChange"] > 0),])
# rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] < 0),])

liPs <- (res1_2[which(res1_2[,"padj"] < 0.05),])
lipS <- (res1_3[which(res1_3[,"padj"] < 0.05),])
lips <- (res1_4[which(res1_4[,"padj"] < 0.05),])
liF <- (res1_F[which(res1_F[,"padj"] < 0.05),])
liPs
lipS
lips
liF

#############################FPKM#######################
myGTF <- read.table("~/R/Eutrema/PS/drought.gtf")
myGTF <- myGTF[,-c(9, 11, 12, 14)]
colnames(myGTF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "geneID", "transcriptID")
myGTF <- mutate(myGTF, length = abs(myGTF$start - myGTF$end))
myGTF <- myGTF %>% filter( grepl("transcript", feature)) %>% arrange(geneID)
head(myGTF)

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
fpkms2 <- fpkm(dds, robust = T)
sampNames <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
colnames(fpkms2) <- sampNames
#writing FPKMs to file of all fpkms
write.table(data.frame("Gene" = rownames(fpkms2), fpkms2), file = "FPKMS.tab", quote = F, row.names = F, sep = "\t")

# lists <- list(Ps=res1_2, pS=res1_3, ps=res1_4, Field=res1_F)
lists <- list(Ps=res1_2, pS=res1_3, ps=res1_4)

getFPKMs <- function(i) {
  diff2 <- i[which(i[,"padj"] < 0.05),]
  diffName2 <- rownames(diff2)

  diffDP2 <- fpkms2[row.names(fpkms2) %in% diffName2, ]
  diffDP2 <- as.data.frame(diffDP2)
  return(diffDP2)
}
for (i in 1:length(lists)) {
fpkmdf <-  getFPKMs(lists[[i]])
  write.table(data.frame("Gene" = rownames(fpkmdf), fpkmdf), file = paste0("diff", names(lists)[i], ".tab"), quote = F, row.names = F, sep = "\t")
}

######################FPKM median length#################
lenMedian <- median(lens)
ddsLenAdj <- dds
mcols(ddsLenAdj)$basepairs <- lenMedian
fpkmLenAdj <- fpkm(ddsLenAdj, robust = T)
head(fpkmLenAdj)
fpkmLenAdj[row.names(fpkmLenAdj) %in% "nXLOC_008023", ]
colnames(fpkmLenAdj) <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
write.table(data.frame("Gene" = rownames(fpkmLenAdj), fpkmLenAdj), file = "FPKMS_LengthAdjusted2.tab", quote = F, row.names = F, sep = "\t")
###############PHEATMAP#############
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition")])
#must have this or else have error msg: Error in check.length("fill") :
#                                         'gpar' element 'fill' must not be length 0
  rownames(df) <- colnames(assay(ntd)[select,])
png("heatmap2.png", height = 1200, width = 1600, res = 200)
pheatmap::pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()


distsRL <- dist(t(assay(rld))) # Calculate distances using transformed (and normalized) counts
mat <- as.matrix(distsRL) # convert to matrix
colnames(mat) <- paste0(condition, ".", colnames(mat))
rownames(mat) <- colnames(mat) # set rownames in the matrix
colnames(mat) = NULL # remove column names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
png("heatmap.png", height = 1200, width = 1600, res = 200)
pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors, cluster_rows = F, cluster_cols = F, legend_breaks = c(0, 50, 100, 150, 200, max(mat)),
         legend_labels = c("0", "50", "100", "150", "200", "Distance (count)"),
         labels_col = rownames(mat))
dev.off()

###################PHYLOSEQ####################
library(phyloseq)
library(ape)
unloadNamespace("Rgraphviz")
unloadNamespace("VennDiagram")
unloadNamespace("gridExtra")
unloadNamespace("grid")
unloadNamespace("cowplot")
vst <- varianceStabilizingTransformation(dds)
vst <- assay(vst)

bigOtu <- otu_table(vst, taxa_are_rows = T)

sampledat <- sample_data(data.frame(row.names=sample_names(bigOtu), condition = condition, stringsAsFactors = F,
                                    sampleID = sampleNO))
physeq <- phyloseq(bigOtu, sampledat)
# plot_bar(physeq)
physeq
ord <- ordinate(physeq, method = "PCoA")
pcoa <- plot_ordination(physeq, ord, type = "samples", color = "condition")
pcoa
ggsave("vst-pcoa.png", pcoa, height = 4, width = 6, dpi = 125)

###Subsetteds dds of cabinet plants only

vst <- varianceStabilizingTransformation(dds)
subsetvst <- vst[, vst$condition %in% c("Ps", "pS", "ps", "PS")]
subsetvst <- assay(subsetvst)
smallOtu <- otu_table(subsetvst, taxa_are_rows = T)
sampleSmall <- c("ps12_S3",
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

sampleDatS <- sample_data(data.frame(row.names=sample_names(smallOtu), condition = rep(c("ps", "Ps", "pS", "PS"), each = 3),
                                     sampleID =sampleSmall, stringsAsFactors = F))
phySmall <- phyloseq(smallOtu, sampleDatS)
ordS <- ordinate(phySmall, method = "PCoA")
pocaS <- plot_ordination(phySmall, ordS, type = "samples", color = "condition")
pocaS
ggsave("vst-pcoaS.png", pocaS, height = 4, width = 6, dpi = 125)


##################PCA#################
library(FactoMineR)
library(factoextra)


# fpkmPCA <- prcomp(fpkms2, center = T,  scale = T)
# nonAdj <- fviz_pca_var(fpkmPCA,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# fviz_pca_var(fpkmPCA,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(1, 2))
# ggsave("fpkm23", nonAdj, dpi = 125, height = 7, width = 7)
# fpkmPCA <- prcomp(fpkms2, center = T,  scale = T)
# fviz_pca_var(fpkmPCA,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# ggsave("fpkm23.png", nonAdj, dpi = 125, height = 7, width = 7)

# fpkmPCA2 <- princomp(fpkms2, cor = F,  scores = T)

#median lengths adjusted FPKMs
fpkmlog <- log2(fpkmLenAdj+1)
fpkmLnAdj <- prcomp(fpkmlog, center = T, scale = T)
# fviz_pca_var(fpkmLnAjPCA,geom = c("point", "text") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(1, 2))f
# fpkmLnAjPCA2 <- princomp(fpkmLenAdj, cor = T, scores = T)
# fviz_pca_var(fpkmLnAjPCA2,geom = c("point", "text") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# PCAdecomp23 <- fviz_pca_var(fpkmLnAjPCA2,geom = c("point", "text") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# ggsave("PSPCA23Decomp.png", PCAdecomp23, dpi = 125, width = 7, height = 7)
# PCAdecomp23nt <- fviz_pca_var(fpkmLnAjPCA2,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# PCAdecomp23nt
# ggsave("PSPCA23Decompnt.png", PCAdecomp23nt, dpi = 125, width = 7, height = 7)
# PCAdecomp12nt <- fviz_pca_var(fpkmLnAjPCA2,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(1, 2))
# PCAdecomp12nt
# ggsave("PSPCA12Decompnt.png", PCAdecomp23nt, dpi = 125, width = 7, height = 7)


# fpkmLnAjPCA2unscaled <- princomp(fpkmLenAdj, cor = F, scores = T)
# PCAdecomp23unscaled <- fviz_pca_var(fpkmLnAjPCA2unscaled,geom = c("point", "text") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# PCAdecomp23unscaled
# ggsave("PSPCA23Decompunscaled.png", PCAdecomp23unsclaed, dpi = 125, width = 7, height = 7)

# fpkmLnAjPCA <- prcomp(fpkmLenAdj, center = F, scale = F)
# fviz_pca_var(fpkmLnAjPCA,geom = c("point", "text") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))

# fpkmLnAjPCAscale <- prcomp(fpkmLenAdj, center = T, scale = T)
# scaled23 <- fviz_pca_var(fpkmLnAjPCAscale,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(2, 3))
# ggsave("fpkmLnAdj23.png", scaled23, dpi = 125, width = 7, height = 7)
# scaled24 <-
fviz_pca_biplot(fpkmLnAdj, geom = c("point"), geom.var = c("point", "text"), axes = c(2, 4), alpha.ind = 0.01, repel = T)
fviz_pca_var(fpkmLnAdj, geom = c("point"), habillage = factor(condition), var.size = 1,
             axes = c(2,4), rug = TRUE, mean.point = T)
# fviz_pca_var(fpkmLnAdj, geom = c("point"), habillage = factor(condition), var.size = 0.1, mean.point.size = 3,
#              axes = c(2,3), mean.point = T, rug = TRUE, xlim = c(-0.5,0.5))
screeplot(fpkmLnAdj)
# ggsave("fpkmLnAdj24.png", scaled23, dpi = 125, width = 7, height = 7)
# scaled12 <- fviz_pca_var(fpkmLnAjPCAscale,geom = c("point") ,habillage = factor(condition),repel = T, col.var = "contrib", axes = c(1, 2))
# ggsave("fpkmLnAdj12.png", scaled12, dpi = 125, width = 7, height = 7)
# plotPCA(fpkmLnAjPCAscale)
# ggbiplot(fpkmPCA,choices = c(2,3), alpha = 0.1)
# ggbiplot(fpkmPCA2,choices = c(2,3), alpha = 0.1)
# plot.PCA(FM_PCA)
#
# plot(fpkmPCA)
# dat <- fpkmLenAdj %>% data.frame
# pca1 <- PCA(dat)
# dat$pc1 <- pca1$ind$coord[, 1]
# dat$pc2 <- pca1$ind$coord[, 2]
# dat$pc3 <- pca1$ind$coord[, 3]
# dat$pc4 <- pca1$ind$coord[, 4]
#
# pca.vars <- pca1$var$coord %>% data.frame
#
# pca.vars$vars <- rownames(pca.vars)
#
# pca.vars.m <- melt(pca.vars, id.vars = "vars")
# p <- ggplot(dat, aes(x = pc2, y = pc4)) +
#   geom_point(alpha=0.2)
# p
# pca <- plotPCA(vsd, intgroup=c("condition"), returnData = T)
# percentVar <- round(100 * attr(pca, "percentVar"))
# ggplot(pca, aes(PC1, PC2, color = condition)) +
#   geom_point(size = 3) +
#   labs(x = paste0("PC1: ", percentVar[1], "% variance"),
#        y = paste0("PC2: ", percentVar[2], "% variance")) +
#   coord_fixed()
#
# fviz_pca_var(fpkmPCA,geom = c("point", "text") ,habillage = factor(condition),repel = F, col.var = "contrib", axes = c(2, 3))
#
# fviz_pca_var(fpkmPCA2,geom = c("point", "text") ,habillage = factor(condition),repel = F, col.var = "contrib", axes = c(2, 3))
##################CAITLIN'S METHOD################

head(fpkmLenAdj)
dim(fpkmLenAdj)
fpkmLA_log <- log2(fpkmLenAdj+1)

pcalog <- prcomp(t(fpkmLA_log), center = T, scale = T)

plot(pcalog, type = "l")
x.var <- pcalog$sdev^2
x.pvar <- x.var/sum(x.var)

plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')

summary(pcalog)
screeplot(pcalog)
fviz_screeplot(pcalog)
fviz_pca_ind(pcalog, habillage = factor(condition), addEllipses = TRUE, ellipse.type = "confidence", axes = c(4, 3))
# var_habillage <- sapply(r1_4$lfc, function(x){ifelse(abs(x) >= 1, "Differential", "Normal")})
# fviz_pca_var(pcalog, geom = "point", alpha = 0.75) #uhhhh...
confplot <- fviz_pca_biplot(pcalog, geom.var = "pooint", axes = c(2, 4), habillage = factor(condition), addEllipses = T,ellipse.type = "confidence")
confplot
ggsave("PS_fpkmMeanPlot.png", confplot, dpi = 125, width = 10, height = 7)


pcalog2 <- prcomp(fpkmLA_log, center = T, scale = T)

plot(pcalog2, type = "l")
x.var <- pcalog2$sdev^2
x.pvar <- x.var/sum(x.var)

plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')

summary(pcalog2)
screeplot(pcalog2)
fviz_screeplot(pcalog2)
fviz_pca_ind(pcalog2, addEllipses = TRUE, ellipse.type = "confidence")
var_habillage <- sapply(r1_4$lfc, function(x){ifelse(abs(x) >= 1, "Differential", "Normal")})
fviz_pca_var(pcalog2, geom = "point", alpha = 0.75)
confplot <- fviz_pca_biplot(pcalog2, alpha.ind = 0.1, geom.ind = "pooint", geom.var = c("point", "text"),
                            axes = c(2, 4), addEllipses = T,ellipse.type = "confidence",
                            repel = T, col.var = "#043d09", col.ind = rep(c("#d4d8c5", "#d4d8c5", "#d4d8c5", "#d4d8c5"), each = 3),
                            addEllipses = T, ellipse.type = "confidence") +
  xlim(-1.1,1.1) + ylim(-1.1,1.1)
confplot

var_plot <- fviz_pca_var(pcalog2, geom = "point", mean.point = T, alpha = 0.75, habillage = factor(condition), axes = c(2, 4),
                         addEllipses = T, ellipse.type = "confidence", title = "Variables PCA, 95% confidence") +
  xlim(-0.5,0.5) + ylim(-0.5,0.5)
var_plot

#############with normal non-adjusted fpkms#################3
head(fpkms2)
dim(fpkms2)
fpkmLog <- log2(fpkms2+1)

pcaFP <- prcomp(fpkmLog, center = T, scale = T)
fviz_screeplot(pcaFP)
var_plot2 <- fviz_pca_var(pcaFP, geom = "point", mean.point = T, alpha = 0.75, habillage = factor(condition), axes = c(2, 3),
                         addEllipses = T, ellipse.type = "confidence", title = "Variables PCA, 95% confidence") #+
#  xlim(-0.5,0.5) + ylim(-0.5,0.5)
var_plot2

##################Another way of analizing the workflow########
library(xlsx)
### PCA gene averages
library(openxlsx)
# setwd("~/scratch/PS/") #working directory goes here
setwd("~/R/Eutrema/PS")
#setwd("~/Documents/mcmaster/phd/rscripts/pca/pca_gene_avg/") # hard drive failure...getting scripts from Dropbox currently...please save this i guess??!?!?!
fpkm <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", row.names=1, header = T) #scp from McMaster cluster
head(fpkm)

#drought <- fpkm[,c(13:28,44:45)]
fpkm.val <- fpkm[,c(1:12)]

ps <- fpkm.val[,1:3]
Ps <- fpkm.val[,4:6]
pS <- fpkm.val[,7:9]
PS <- fpkm.val[,10:12]

pS <- apply(pS, 1, median)
Ps <- apply(Ps, 1, median)
ps <- apply(ps, 1, median)
PS <- apply(PS, 1, median)

ps_med <- cbind(pS, Ps, ps, PS)

ps.log <- log2(ps_med+1)
## all genes
pca<- prcomp(ps.log, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)
#log.expr <- log2(expr+1)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)
#rv <- rowVars(as.matrix(log.expr))
#names(rv) <- rownames(expr)
#select <- order(rv, decreasing = TRUE)[1:2000]
#pca <- prcomp(log.expr[select,], .scale=TRUE, center=TRUE)

dim(exprVals)
dim(sampleVals)

#using covariance
ps.cov <- cov(ps.log)
pca.cov <- prcomp(ps.cov, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca.cov$rotation  # The loadings are here

sdev <- pca.cov$sdev

screeplot(pca.cov)
#log.expr <- log2(expr+1)

plot(pca.cov, type = "l")
summary(pca.cov)
exprVals<-data.frame(pca.cov$x)
sampleVals<-data.frame(pca.cov$rotation)
#rv <- rowVars(as.matrix(log.expr))
#names(rv) <- rownames(expr)
#select <- order(rv, decreasing = TRUE)[1:2000]
#pca <- prcomp(log.expr[select,], .scale=TRUE, center=TRUE)

dim(exprVals)
dim(sampleVals)
expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1, header = T)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)
colnames(FPKM)  <- sampNames
##################
# Write to excel #
##################

pca_list <- list(neg2neg4_genes,  neg2pos4_genes, pos2neg4_genes, pos2pos4_genes)
names(pca_list) = c("neg2neg4",  "neg2pos4", "pos2neg4", "pos2pos4")
#write.xlsx(pca_list, file = "pca_genelists_updated2.xlsx", rowNames=TRUE)
#####
## Visualizing PCAs...

#samples <-colnames(ps_med)
samples <- c("pS", "Ps", "ps", "PS")
coords <- data.frame(X=rep(0, 4), Y=rep(0, 4),sampleVals, Samples = samples)
coords$Treatment <- factor(coords$Samples, c("pS", "Ps", "ps", "PS"))

library(ggplot2)
library(ggrepel)
theme_set(theme_bw())
exprVals
labels_pos2neg4 <- ifelse(rownames(exprVals) %in% rownames(pos2neg4_58), "Yes", "No")
labels_pos2pos4 <- ifelse(rownames(exprVals) %in% rownames(pos2pos4_58), "Yes", "No")
labels_neg2neg4 <- ifelse(rownames(exprVals) %in% rownames(neg2neg4_58), "Yes", "No")
exprVals_lab <- cbind(exprVals, labels_pos2neg4, labels_pos2pos4, labels_neg2neg4)

### Labelling plot!
#covariance matrix...is stupid don't do it
#PC2, PC4 nothing good
{ggplot(exprVals, aes(x = PC2, y = PC4))  + 
    geom_segment(data = coords, aes(x = X, y = Y, xend = PC2, yend = PC4, color = Treatment), 
                 arrow = arrow(), size = 1) +
#    geom_point(data = pca$x, aes(x = PC2, y = PC4), shape=19, alpha=0.3) +
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))}

#PC1, PC2 ... nothing good
    {ggplot(exprVals, aes(x = PC1, y = PC2))  + 
    geom_segment(data = coords, aes(x = X, y = Y, xend = PC1, yend = PC2, color = Treatment), 
                 arrow = arrow(), size = 1) +
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))}

#PC3, PC2 not responsive to low sulfur
    {ggplot(exprVals, aes(x = PC3, y = PC2))  + 
    geom_segment(data = coords, aes(x = X, y = Y, xend = PC3, yend = PC2, color = Treatment), 
                 arrow = arrow(), size = 1) +
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))}

#no covariance
median <- (pc24plot_labelled <- ggplot(exprVals, aes(x = PC2, y = PC4)) +
    geom_point(shape=19, alpha=0.3) +
    geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
    geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
    # geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Treatment), arrow=arrow(), size=1) +
    geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=1.5, colour="black") +
    geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), hjust=1.5, colour="white") +
    # geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
    xlab("PC2 (0.74%)") + ylab("PC4 (0.18%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)))
median
ggsave("PSMed.png", median, dpi = 100, height = 4, width = 6)
#png("figs/pc24_labelled_recoloured.png", bg="transparent")
#print(pc24plot_labelled)
#dev.off()

(pc24plot_pos2neg4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
    # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
    # geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
    geom_point(data = exprVals_lab[exprVals_lab$labels_pos2neg4 == "Yes",], colour = "mediumpurple1", alpha = 0.5) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size =1) +
    # geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
    # geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
    # geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
    xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
    scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    theme_set(theme_bw(base_size=16)))

(pc24plot_pos2pos4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
    #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
    #  geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
    geom_point(data = exprVals_lab[exprVals_lab$labels_pos2pos4 == "Yes",], colour = "mediumpurple1", alpha = 0.7) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
    #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
    #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
    #  geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
    xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
    scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    theme_set(theme_bw(base_size=16)))

(pc24plot_neg2neg4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
    #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
    #  geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
    geom_point(data = exprVals_lab[exprVals_lab$labels_neg2neg4 == "Yes",], colour = "mediumpurple1", alpha = 0.4) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
    #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
    #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
    #  geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
    xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
    scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    theme_set(theme_bw(base_size=16)))

#ggsave("pos2pos4_labelled.pdf", pc24plot_pos2pos4lab)
#ggsave("pos2neg4_labelled.pdf", pc24plot_pos2neg4lab)
#ggsave("neg2neg4_labelled.pdf", pc24plot_neg2neg4lab)



(pc24plot_recol <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="#3366FF", size=4) +
    # geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=2) +
    geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#3366FF", size=4) +
    geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
    geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), vjust=-1, colour="white") +
    xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
    #scale_colour_manual(values=c("#fd8a19", "#fcfc35", "#f91301", "#febfcb"))
    #scale_colour_manual(values=c( "#febfcb", "#fcfc35", "#f91301", "#fd8a19" ))
    #scale_colour_manual(values=c( "coral2", "#fcfc35", "#f91301", "#fd8a19" ))
    scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    theme_set(theme_bw(base_size=16)))
#ggsave(pc24plot_recol, file="pc24_recoloured.pdf")

pc24plot_nolab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
  # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="#3366FF", size=4) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=2) +
  #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#3366FF", size=4) +
  #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
  #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), vjust=-1, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
  #scale_colour_manual(values=c("#fd8a19", "#fcfc35", "#f91301", "#febfcb"))
  #scale_colour_manual(values=c( "#febfcb", "#fcfc35", "#f91301", "#fd8a19" ))
  #scale_colour_manual(values=c( "coral2", "#fcfc35", "#f91301", "#fd8a19" ))
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16))
pc24plot_nolab
ggsave(pc24plot_nolab, file="pc24_recoloured_noLabels.pdf")





## loop for all PC1 combinations
#for (tryThis in 2:4){
# toUse<-noquote(paste("PC", tryThis, sep=""))
#fileName<-paste("figs/PC1",toUse, ".pdf",sep="")
#pdf(fileName)
# print(pc12plot <- ggplot(coords, aes_string("PC1", toUse)) +
# geom_point(aes(color = Samples), size=5) +
    #         geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=toUse, colour=Samples), arrow=arrow(), size=1))# +
# coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)))
#dev.off()
#}

pc12plot <- ggplot(coords, aes_string("PC1", "PC2")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC2, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))
pc12plot

pc13plot <- ggplot(coords, aes_string("PC1", "PC3")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC3, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))
pc13plot

pc14plot <- ggplot(coords, aes_string("PC1", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))
pc14plot

pc23plot <- ggplot(coords, aes_string("PC2", "PC3")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC3, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))
pc23plot

pc24plot <- ggplot(coords, aes_string("PC2", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))
pc24plot

pc34plot <- ggplot(coords, aes_string("PC3", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC3, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))

#pdf("figs/PC1PC2.pdf")
#print(pc12plot)
#dev.off()
#
#pdf("figs/PC1PC3.pdf")
#print(pc13plot)
#dev.off()
#
#pdf("figs/PC1PC4.pdf")
#print(pc14plot)
#dev.off()
#
#pdf("figs/PC2PC3.pdf")
#print(pc23plot)
#dev.off()
#
#pdf("figs/PC2PC4.pdf")
#print(pc24plot)
#dev.off()
#
#pdf("figs/PC3PC4.pdf")
#print(pc34plot)
#dev.off()

write.xlsx(p2n4, file = "mediantopPCA_loadings.xlsx", sheetName = "Pos2Neg4", append = F)
write.xlsx(p2p4, file = "mediantopPCA_loadings.xlsx", sheetName = "Pos2Pos4", append = T)
write.xlsx(n2p4, file = "mediantopPCA_loadings.xlsx", sheetName = "Neg2Pos4", append = T)
write.xlsx(n2n4, file = "mediantopPCA_loadings.xlsx", sheetName = "Neg2Neg4", append = T)

#########################GET MEANS PCA############################
fpkm.val <- fpkm[,c(1:12)]

ps <- fpkm.val[,1:3]
Ps <- fpkm.val[,4:6]
pS <- fpkm.val[,7:9]
PS <- fpkm.val[,10:12]

pS <- (pS=rowMeans(pS))
Ps <- (Ps=rowMeans(Ps))
ps <- (ps=rowMeans(ps))
PS <- (PS=rowMeans(PS))

ps_mean <- cbind(ps, pS, Ps, PS)

ps.log <- log2(ps_mean+1)

## all genes
pca<- prcomp(ps.log, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)
#log.expr <- log2(expr+1)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)
#rv <- rowVars(as.matrix(log.expr))
#names(rv) <- rownames(expr)
#select <- order(rv, decreasing = TRUE)[1:2000]
#pca <- prcomp(log.expr[select,], .scale=TRUE, center=TRUE)

dim(exprVals)
dim(sampleVals)

expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1)
#

#Top 50 loading genes for each principle component instead of intersections (which doesn't make much sense becuase
#they would be some kind of middle ground for each PC instead of the most contribution to each PC)
pos2 <- exprVals[which(exprVals[,2] > 0), c(2,4)]
pos2_50 <- pos2[order(pos2$PC2, decreasing = T)[1:200],]
pos2_fpkm <- FPKM[rownames(FPKM) %in% rownames(pos2_50),]  
pos2_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2_50),]
pos2_genes <- merge(pos2_50, pos2_genes[,18:ncol(pos2_genes)], by = 0)
rownames(pos2_genes) <- pos2_genes$Row.names
pos2_genes <- subset(pos2_genes, select = -c(Row.names))
p2 <- cbind(pos2_fpkm, pos2_genes)
p2 <- p2[order(p2$PC2, decreasing = T), ]


neg2 <- exprVals[which(exprVals[,2] < 0), c(2,4)]
neg2_50 <- neg2[order(neg2$PC2, decreasing = T)[1:200],]
neg2_fpkm <- FPKM[rownames(FPKM) %in% rownames(neg2_50),]  
neg2_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2_50),]
neg2_genes <- merge(neg2_50, neg2_genes[,18:ncol(neg2_genes)], by = 0)
rownames(neg2_genes) <- neg2_genes$Row.names
neg2_genes <- subset(neg2_genes, select = -c(Row.names))
n2 <- cbind(neg2_fpkm, neg2_genes)
n2 <- n2[order(n2$PC2, decreasing = T),]

pos4 <- exprVals[which(exprVals[,4] > 0), c(2,4)]
pos4_50 <- pos4[order(pos4$PC4, decreasing = T)[1:200],]
pos4_fpkm <- FPKM[rownames(FPKM) %in% rownames(pos4_50),]  
pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos4_50),]
pos4_genes <- merge(pos4_50, pos4_genes[,18:ncol(pos4_genes)], by = 0)
rownames(pos4_genes) <- pos4_genes$Row.names
pos4_genes <- subset(pos4_genes, select = -c(Row.names))
p4 <- cbind(pos4_fpkm, pos4_genes)
p4 <- p4[order(p4$PC4, decreasing = T),]

neg4 <- exprVals[which(exprVals[,4] < 0), c(2,4)]
neg4_50 <- neg4[order(neg4$PC4, decreasing = T)[1:200],]
neg4_fpkm <- FPKM[rownames(FPKM) %in% rownames(neg4_50),]  
neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg4_50),]
neg4_genes <- merge(neg4_50, neg4_genes[,18:ncol(neg4_genes)], by = 0)
rownames(neg4_genes) <- neg4_genes$Row.names
neg4_genes <- subset(neg4_genes, select = -c(Row.names))
n4 <- cbind(neg4_fpkm, neg4_genes)
n4 <- n4[order(n4$PC4, decreasing = T),]

write.xlsx2(p2, file = "meanPCA_loadings.xlsx", sheetName = "Pos2", append = F, row.names = T)
write.xlsx2(n2, file = "meanPCA_loadings.xlsx", sheetName = "Neg2", append = T, row.names = T)
write.xlsx2(p4, file = "meanPCA_loadings.xlsx", sheetName = "Pos4", append = T, row.names = T)
write.xlsx2(n4, file = "meanPCA_loadings.xlsx", sheetName = "Neg4", append = T, row.names = T)
####################
# Pos PC2, Neg PC4 #
####################
# top 50?
pos2neg4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] < 0), c(2,4)]

#top loading in each
pos2_200 <-pos2neg4[order(pos2neg4$PC2, decreasing=TRUE)[1:200],]
neg4_200 <- pos2neg4[order(pos2neg4$PC4, decreasing = FALSE)[1:200],]

#loading value
pos2neg4_58 <- pos2_200[rownames(pos2_200) %in% rownames(neg4_200),]

#fpkm values
pos2neg4_fpkm <- FPKM[rownames(FPKM) %in% rownames(pos2neg4_58),]
#annotation fields
pos2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2neg4_58),]
#pos2neg4_genes <- cbind(pos2neg4_58, pos2neg4_genes[,7:ncol(pos2neg4_genes)])
pos2neg4_genes <- merge(pos2neg4_58, pos2neg4_genes[,18:ncol(pos2neg4_genes)], by=0)
rownames(pos2neg4_genes) <- pos2neg4_genes$Row.names
pos2neg4_genes <- subset(pos2neg4_genes, select = -c(Row.names))
p2n4 <- cbind(pos2neg4_fpkm, pos2neg4_genes)
#write.csv(pos2neg4_genes, file="pos2neg4_genes.csv")

####################
# Pos PC2, Pos PC4 #
####################
# top 200?
# 
pos2pos4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] > 0), c(2,4)]

#top loading in each
pos2_200 <-pos2pos4[order(pos2pos4$PC2, decreasing=TRUE)[1:200],]
pos4_200 <- pos2pos4[order(pos2pos4$PC4, decreasing = FALSE)[1:200],]

#loading value
pos2pos4_58 <- pos2_200[rownames(pos2_200) %in% rownames(pos4_200),]

#fpkm values
pos2pos4_fpkm <- FPKM[rownames(FPKM) %in% rownames(pos2pos4_58),]
#annotation fields
pos2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2pos4_58),]
#pos2pos4_genes <- cbind(pos2pos4_58, pos2pos4_genes[,7:ncol(pos2pos4_genes)])
pos2pos4_genes <- merge(pos2pos4_58, pos2pos4_genes[,18:ncol(pos2pos4_genes)], by=0)
rownames(pos2pos4_genes) <- pos2pos4_genes$Row.names
pos2pos4_genes <- subset(pos2pos4_genes, select = -c(Row.names))
p2p4 <- cbind(pos2pos4_fpkm, pos2pos4_genes)
# 
# write.csv(pos2pos4_genes, file="pos2pos4_genes.csv")
# ####################
# # Neg PC2, Pos PC4 #
# ####################
# # top 200?
# 
neg2pos4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] > 0), c(2,4)]

#top loading in each
neg2_200 <-neg2pos4[order(neg2pos4$PC2, decreasing=TRUE)[1:200],]
pos4_200 <- neg2pos4[order(neg2pos4$PC4, decreasing = FALSE)[1:200],]

#loading value
neg2pos4_58 <- neg2_200[rownames(neg2_200) %in% rownames(pos4_200),]

#fpkm values
neg2pos4_fpkm <- FPKM[rownames(FPKM) %in% rownames(neg2pos4_58),]
#annotation fields
neg2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2pos4_58),]
#neg2pos4_genes <- cbind(neg2pos4_58, neg2pos4_genes[,7:ncol(neg2pos4_genes)])
neg2pos4_genes <- merge(neg2pos4_58, neg2pos4_genes[,18:ncol(neg2pos4_genes)], by=0)
rownames(neg2pos4_genes) <- neg2pos4_genes$Row.names
neg2pos4_genes <- subset(neg2pos4_genes, select = -c(Row.names))
n2p4 <- cbind(neg2pos4_fpkm, neg2pos4_genes)
# 
# #write.csv(neg2pos4_genes, file="neg2pos4_genes.csv")
# ####################
# # Neg PC2, Neg PC4 #
# ####################
# # top 200?
# 
neg2neg4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] < 0), c(2,4)]

#top loading in each
neg2_200 <-neg2neg4[order(neg2neg4$PC2, decreasing=TRUE)[1:200],]
neg4_200 <- neg2neg4[order(neg2neg4$PC4, decreasing = FALSE)[1:200],]

#loading value
neg2neg4_58 <- neg2_200[rownames(neg2_200) %in% rownames(neg4_200),]

#fpkm values
neg2neg4_fpkm <- FPKM[rownames(FPKM) %in% rownames(neg2neg4_58),]
#annotation fields
neg2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2neg4_58),]
#neg2neg4_genes <- cbind(neg2neg4_58, neg2neg4_genes[,7:ncol(neg2neg4_genes)])
neg2neg4_genes <- merge(neg2neg4_58, neg2neg4_genes[,18:ncol(neg2neg4_genes)], by=0)
rownames(neg2neg4_genes) <- neg2neg4_genes$Row.names
neg2neg4_genes <- subset(neg2neg4_genes, select = -c(Row.names))
n2n4 <- cbind(neg2neg4_fpkm, neg2neg4_genes)
# 

pos2neg4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] < 0), c(2,4)]

pos2_200 <-pos2neg4[order(pos2neg4$PC2, decreasing=TRUE)[1:200],]
neg4_200 <- pos2neg4[order(pos2neg4$PC4, decreasing = FALSE)[1:200],]

pos2neg4_58 <- pos2_200[rownames(pos2_200) %in% rownames(neg4_200),]
#pos2neg4_58 <- neg4_200[na.omit(match(rownames(pos2_200), rownames(neg4_200))),] # locations in pos2_200
#na.omit(match(rownames(pos2neg4_58), rownames(expr_annot)))
#pos2neg4_genes <- expr_annot[na.omit(match(rownames(pos2neg4_58), rownames(expr_annot))),]
pos2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2neg4_58),]
#pos2neg4_genes <- cbind(pos2neg4_58, pos2neg4_genes[,7:ncol(pos2neg4_genes)])
pos2neg4_genes <- merge(pos2neg4_58, pos2neg4_genes[,7:ncol(pos2neg4_genes)], by=0)
rownames(pos2neg4_genes) <- pos2neg4_genes$Row.names
pos2neg4_genes <- subset(pos2neg4_genes, select = -c(Row.names))
#write.csv(pos2neg4_genes, file="pos2neg4_genes.csv")

####################
# Pos PC2, Pos PC4 #
####################
# top 50?

pos2pos4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] > 0), c(2,4)]

pos2_200 <-pos2pos4[order(pos2pos4$PC2, decreasing=TRUE)[1:300],]
pos4_200 <- pos2pos4[order(pos2pos4$PC4, decreasing = TRUE)[1:300],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#pos2pos4_top <- pos2_200[na.omit(match(rownames(pos4_200), rownames(pos2_200))),] # 33
pos2pos4_top <- pos2_200[rownames(pos2_200) %in% rownames(pos4_200),]
pos2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2pos4_top),]
#pos2pos4_genes <- cbind(pos2pos4_top, pos2pos4_genes[,7:ncol(pos2pos4_genes)])
pos2pos4_genes <- merge(pos2pos4_top, pos2pos4_genes[,7:ncol(pos2pos4_genes)], by=0)
rownames(pos2pos4_genes) <- pos2pos4_genes$Row.names
pos2pos4_genes <- subset(pos2pos4_genes, select = -c(Row.names))

write.csv(pos2pos4_genes, file="pos2pos4_genes.csv")
####################
# Neg PC2, Pos PC4 #
####################
# top 50?

neg2pos4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] > 0), c(2,4)]

neg2_200 <-neg2pos4[order(neg2pos4$PC2, decreasing=FALSE)[1:200],]
pos4_200 <- neg2pos4[order(neg2pos4$PC4, decreasing = TRUE)[1:200],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#neg2pos4_top <- neg2_200[na.omit(match(rownames(pos4_200), rownames(neg2_200))),] # 78
neg2pos4_top <- neg2_200[rownames(neg2_200) %in% rownames(pos4_200),]

neg2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2pos4_top),]
#neg2pos4_genes <- cbind(neg2pos4_top, neg2pos4_genes[,7:ncol(neg2pos4_genes)])
neg2pos4_genes <- merge(neg2pos4_top, neg2pos4_genes[,7:ncol(neg2pos4_genes)], by=0)
row.names(neg2pos4_genes) <- neg2pos4_genes$Row.names
neg2pos4_genes <- subset(neg2pos4_genes, select = -c(Row.names))

#write.csv(neg2pos4_genes, file="neg2pos4_genes.csv")
####################
# Neg PC2, Neg PC4 #
####################
# top 50?

neg2neg4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] < 0), c(2,4)]

neg2_300 <-neg2neg4[order(neg2neg4$PC2, decreasing=FALSE)[1:300],]
neg4_300 <- neg2neg4[order(neg2neg4$PC4, decreasing = FALSE)[1:300],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#neg2neg4_top <- neg2_200[na.omit(match(rownames(neg4_200), rownames(neg2_200))),] # 78
neg2neg4_top <- neg2_300[rownames(neg2_300) %in% rownames(neg4_300),]

# neg2_200[rownames(neg2_200) %in% rownames(pos4_200),]
neg2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2neg4_top),]
#neg2neg4_genes <- cbind(neg2neg4_top, neg2neg4_genes[,7:ncol(neg2neg4_genes)])
neg2neg4_genes <- merge(neg2neg4_top, neg2neg4_genes[,7:ncol(neg2neg4_genes)], by=0)
row.names(neg2neg4_genes) <- neg2neg4_genes$Row.names
neg2neg4_genes <- subset(neg2neg4_genes, select = -c(Row.names))


#write.csv(neg2neg4_genes, file="neg2neg4_genes.csv")

##################
# Write to excel #
##################

pca_list <- list(neg2neg4_genes,  neg2pos4_genes, pos2neg4_genes, pos2pos4_genes)
names(pca_list) = c("neg2neg4",  "neg2pos4", "pos2neg4", "pos2pos4")
#write.xlsx(pca_list, file = "pca_genelists_updated2.xlsx", rowNames=TRUE)
#####
## Visualizing PCAs...

#samples <-colnames(ps_med)
samples <- c("ps", "pS", "Ps", "PS")
coords <- data.frame(X=rep(0, 4), Y=rep(0, 4),sampleVals, Samples = samples)
coords$Treatment <- factor(coords$Samples, c("ps", "pS", "Ps", "PS"))

theme_set(theme_bw())

labels_pos2neg4 <- ifelse(rownames(exprVals) %in% rownames(pos2neg4_58), "Yes", "No")
labels_pos2pos4 <- ifelse(rownames(exprVals) %in% rownames(pos2pos4_top), "Yes", "No")
labels_neg2neg4 <- ifelse(rownames(exprVals) %in% rownames(neg2neg4_top), "Yes", "No")
exprVals_lab <- cbind(exprVals, labels_pos2neg4, labels_pos2pos4, labels_neg2neg4)

### Labelling plot!
library(ggrepel)
library(shadowtext)

set.seed(51)
(mean <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Treatment), arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
    plot.points(points, exprVals) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), nudge_y = 0.07) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), nudge_y = 0.07) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10018244m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10018244m.g" ,], aes(PC2,PC4, label = "SULTR1;2"), nudge_y = 0.07) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10005068m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10005068m.g" ,], aes(PC2,PC4, label = "Rhodanese"), nudge_y = 0.07) +
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    xlab("PC2 (0.78%)") + ylab("PC4 (0.22%)") +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)))
mean
median
# ggsave("PSMean.pdf", mean, dpi = 250, height = 6, width = 8)
ggsave("PSMeannoLabel.pdf", mean, dpi = 250, height = 6, width = 8)
write.xlsx(p2n4, file = "meantopPCA_loadings.xlsx", sheetName = "Pos2Neg4", append = F)
write.xlsx(p2p4, file = "meantopPCA_loadings.xlsx", sheetName = "Pos2Pos4", append = T)
write.xlsx(n2p4, file = "meantopPCA_loadings.xlsx", sheetName = "Neg2Pos4", append = T)
write.xlsx(n2n4, file = "meantopPCA_loadings.xlsx", sheetName = "Neg2Neg4", append = T)
####################FIND Top 50 loading genes###################

medPCA <- PCA(log2(ps_med + 1), graph = F)
meanpCA <- PCA(log2(ps_mean + 1), graph = F)

fviz_pca_ind(medPCA, col.ind = "cos2", 
alpha.ind = 0.5, label = "none", axes = c(2,4)) +
scale_color_gradient2(low = "yellow", mid = "blue", high = "red", midpoint = 0.5)

head(medPCA$ind$cos2)
head(medPCA$ind$contrib)

fviz_contrib(medPCA, choice = "ind", axes = 2:4, top = 50)

