library(DESeq2)                        #DEG analysis
                                       # library(plyr) #data manipulation
                                       # library(gplots) #graphing and visualization
                                       # library(RColorBrewer) #graphing and visualization
library(genefilter)                    #GFF maker
library(dplyr)                         #graphing
library(geneplotter)                   #graphing
library(ggplot2)                       #graphing
library(GenomicFeatures)               #GFF functions makeTxdbFromGff
library(reshape2)                      #graphing
library(apeglm)                        #better than ashr
library(venn)                          #used for venn
                                       # library(VennDiagram) #used for venn
library(topGO)                         #used for GO enrichment (results not informative)
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
# rld <- rlogTransformation(dds, blindmeanComparison.pdf", plotted, dpi = 250, height = 14, width = 22) png("PCA_samples.png", height = 600, width = 1000, res = 125)
# plotPCA.san(rld, intgroup=c("condition"))
# dev.off()
# vsd <- vst(dds, blind = F)
# plotPCA.san(vsd, intgroup = c("condition"))
# DESeq2::plotPCA(vsd, intgroup = c("condition"))
#https://www.biostars.org/p/243695/
library(ggrepel)

#how many genes are significantly expressed
# png("PCA_samplesno.png", height = 700, width = 1000, res = 125)
# DESeq2::plotPCA(rldS, intgroup=c("condition"))
# dev.off()
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_PS"  "condition_Ps_vs_PS" "condition_pS_vs_PS" "condition_ps_vs_PS"

res1_2 <- lfcShrink(dds, coef = "condition_Ps_vs_PS", type = "apeglm")
res1_3 <- lfcShrink(dds, coef = "condition_pS_vs_PS", type = "apeglm")
res1_4 <- lfcShrink(dds, coef = "condition_ps_vs_PS", type = "apeglm")
res.Ps <- results(dds, contrast = c("condition", "Ps", "PS"))
res.pS <- results(dds, contrast = c("condition", "pS", "PS"))
res.ps <- results(dds, contrast = c("condition", "ps", "PS"))
#########################PLOTTING###############################################
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

{ggplot(IPS, aes(x = sample, y = lfc)) +
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
          axis.title.y = element_blank())}

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
plotted
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

dat <- list(Ps=GOres1_2, pS=GOres1_3, ps=GOres1_4)
names <- c("Ps", "pS", "ps")

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


liPs <- (res1_2[which(res1_2[,"padj"] < 0.05),])
lipS <- (res1_3[which(res1_3[,"padj"] < 0.05),])
lips <- (res1_4[which(res1_4[,"padj"] < 0.05),])
liF <- (res1_F[which(res1_F[,"padj"] < 0.05),])
liPs
lipS
lips
liF

#############################FPKM#######################
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
# for (i in 1:length(lists)) {
# fpkmdf <-  getFPKMs(lists[[i]])
#   write.table(data.frame("Gene" = rownames(fpkmdf), fpkmdf), file = paste0("diff", names(lists)[i], ".tab"), quote = F, row.names = F, sep = "\t")
# }

######################FPKM median length#################
lenMedian <- median(lens)
ddsLenAdj <- dds
mcols(ddsLenAdj)$basepairs <- lenMedian
fpkmLenAdj <- fpkm(ddsLenAdj, robust = T)
head(fpkmLenAdj)
fpkmLenAdj[row.names(fpkmLenAdj) %in% "nXLOC_008023", ]
colnames(fpkmLenAdj) <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
# write.table(data.frame("Gene" = rownames(fpkmLenAdj), fpkmLenAdj), file = "FPKMS_LengthAdjusted2.tab", quote = F, row.names = F, sep = "\t")
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
###   CAITLIN'S METHOD
library(xlsx)
### PCA gene averages
library(openxlsx)
# setwd("~/scratch/PS/") #working directory goes here
setwd("~/R/Eutrema/PS")
#setwd("~/Documents/mcmaster/phd/rscripts/pca/pca_gene_avg/") # hard drive failure...getting scripts from Dropbox currently...please save this i guess??!?!?!
# fpkm <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", row.names=1, header = T) #scp from McMaster cluster
fpkm <- fpkms2
head(fpkm)

#########################GET MEANS PCA############################
fpkm.val <- fpkm[,c(1:12)]             # 

fpkm.val <- log2(fpkm.val + 1)
ps <- fpkm.val[,1:3]
Ps <- fpkm.val[,4:6]
pS <- fpkm.val[,7:9]
PS <- fpkm.val[,10:12]

pS <- (pS=rowMeans(pS))
Ps <- (Ps=rowMeans(Ps))
ps <- (ps=rowMeans(ps))
PS <- (PS=rowMeans(PS))

ps_mean <- cbind(ps, pS, Ps, PS)


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
neg2_50 <- neg2[order(neg2$PC2)[1:200],]
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
neg4_50 <- neg4[order(neg4$PC4)[1:200],]
neg4_fpkm <- FPKM[rownames(FPKM) %in% rownames(neg4_50),]  
neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg4_50),]
neg4_genes <- merge(neg4_50, neg4_genes[,18:ncol(neg4_genes)], by = 0)
rownames(neg4_genes) <- neg4_genes$Row.names
neg4_genes <- subset(neg4_genes, select = -c(Row.names))
n4 <- cbind(neg4_fpkm, neg4_genes)
n4 <- n4[order(n4$PC4, decreasing = T),]

# write.xlsx2(p2, file = "meanPCA_loadings.xlsx", sheetName = "Pos2", append = F, row.names = T)
# write.xlsx2(n2, file = "meanPCA_loadings.xlsx", sheetName = "Neg2", append = T, row.names = T)
# write.xlsx2(p4, file = "meanPCA_loadings.xlsx", sheetName = "Pos4", append = T, row.names = T)
# write.xlsx2(n4, file = "meanPCA_loadings.xlsx", sheetName = "Neg4", append = T, row.names = T)
# 
# write(rownames(p2), "p2.names", sep = "\n")
# write(rownames(n2), "n2.names", sep = "\n")
# write(rownames(p4), "p4.names", sep = "\n")
# write(rownames(n4), "n4.names", sep = "\n")

### Labelling plot!
library(ggrepel)
library(shadowtext)
library(stringr)

# preidction lists
pred.names <- scan("~/R/Eutrema/PS/crema/pcaLoading/pred.names", what = "character")
filtName <- scan("~/R/Eutrema/PS/filter.name", what = "character")
allGpredEq1 <- read.csv("~/R/Eutrema/PS/crema/allG/final_ensemble_predictions.csv", header = T, row.name = 1)

coords <- data.frame(X=rep(0,4), Y=rep(0,4), sampleVals, Treatment=c("ps", "pS", "Ps", "PS"))
araConvert <- read.table("~/Eutrema/FPKM/eutremaToArabidopsis.names")
names(araConvert) <- c("Eutr", "AraT")
araTr <- araConvert
araTr <- araConvert[complete.cases(araConvert),]
#araTr$Ara <- gsub("[:punct:]+\\d", "", araTr$AraT)
araTr$Ara <- gsub("\\.[0-9]*", "", araTr$AraT)


Sassim <- read.xlsx("~/R/Eutrema/Arabi2016/GSLandSassim.xlsx", sheet = 1)
Sassim$Ara <- toupper(Sassim$Ara)
Sassim <- data.frame(sapply(Sassim, str_trim))
Sassim <- merge(Sassim, araTr, by = "Ara")
# write.table(Sassim, file = "/home/lucy/R/Eutrema/MapMan/data/Sassim.txt", quote = T, sep = "\t")
GSL <- read.xlsx("~/R/Eutrema/Arabi2016/GSLandSassim.xlsx", sheet = 2)
GSL$Ara <- toupper(GSL$Ara)
GSL <- data.frame(sapply(GSL, str_trim))
GSL <- merge(GSL, araTr, by = "Ara")
# write.table(GSL, "/home/lucy/R/Eutrema/MapMan/data/GSL.txt", quote = T, sep = "\t")

exprVals$Type <- ifelse(rownames(exprVals) %in% GSL$Eutr, "GSL", ifelse(rownames(exprVals) %in% Sassim$Eutr, "S assimilation", "other"))
exprVals$Type <- factor(exprVals$Type, levels = c("other", "GSL", "S assimilation"))
rbindGS <- rbind(GSL, Sassim)
exprVals <- merge(exprVals, rbindGS, all.x = T, by.x = 0, by.y = "Eutr")
# rownames(exprVals) <- exprVals$Row.names

exprVals <- exprVals[order(exprVals$Type),]

# lncRNA stuff
allGlnc <- rownames(allGpredEq1[which(allGpredEq1$prediction == 1),])
pcaSum <- as.data.frame(summary(pca)$importance)
set.seed(51)
{meanplot <- ggplot(exprVals, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = Type), shape=21, color = "None") +
    scale_fill_manual(values = c( "gray", "coral3", "chartreuse4")) +
    #     geom_point(data = exprVals[which(toupper(rownames(exprVals)) %in% toupper(allGlnc)) ,], 
    #                colour="#9fcc2e", size=2, alpha = 0.7) +    
    #     geom_point(data = exprValsGS, aes(x = PC2, y = PC4, group = class), color = c("red", "blue"), size = 2, alpha = 0.7) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    geom_label_repel(data = exprVals[which(abs(exprVals$PC2) >= 0.1 | abs(exprVals$PC4) >= 0.1),], aes(label = gene_name, fill = Type), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
    xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",pcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,2), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}
# ggsave("PSMean.pdf", mean, dpi = 250, height = 6, width = 8)
# ggsave("~/R/Eutrema/PS/pics/PSMeanPred.pdf", mean, dpi = 250, height = 6, width = 8)
ggsave("~/R/Eutrema/PS/pics/PSMeanPC24GSLSassim.pdf", meanplot, dpi = 150, height = 10, width = 12)  # 
# anti aliasing of Cairo type image raster
# ggsave("~/R/Eutrema/PS/pics/PSMeanPred.png", meanplot, dpi = 450, height = 4, width = 5, type="cairo") 
# write.xlsx(p2n4, file = "meantopPCA_loadings.xlsx", sheetName = "Pos2Neg4", append = F)
# write.xlsx(p2p4, file = "meantopPCA_loadings.xlsx", sheetName = "Pos2Pos4", append = T)
# write.xlsx(n2p4, file = "meantopPCA_loadings.xlsx", sheetName = "Neg2Pos4", append = T)
# write.xlsx(n2n4, file = "meantopPCA_loadings.xlsx", sheetName = "Neg2Neg4", append = T)

# looking at the lipid metabolism genes in the PS libraries?
# grep for MGD DGD SQD and acyl-transferase
lipidGenes <- read.csv("~/R/Eutrema/PS/lipids.csv", header = F)
colnames(lipidGenes) <- c("Eutr", "transcript_id", "scaffold", "begin", "end", "strand", 
                          "AraT", "name","class","family")

lipidGenes$Ara <- gsub("\\.[0-9]*", "", lipidGenes$AraT)
exprVals <- merge(exprVals, lipidGenes[,c(1,8:10)], by.x = "Row.names", by.y = "Eutr", all.x = T)

exprVals <- exprVals[rev(order(exprVals$family)),]
{meanplot <- ggplot(exprVals, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = family), shape=21, color = "None") +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    geom_label_repel(data = exprVals[which(abs(exprVals$PC2) >= 0.1 | abs(exprVals$PC4) >= 0.1),], aes(label =name, fill = family), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
    xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",pcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,2), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}
("~/R/Eutrema/PS/pics/PSMeanPC24lipids.pdf", meanplot, dpi = 150, height = 12, width = 14)  # 

theGenes <- read.table("~/R/Eutrema/PS/Genes", sep = "\t")
colnames(theGenes) <- c("family", "name", "gene_id", "Function")

exprPlusGenes <- merge(exprVals, theGenes, by.x = "Row.names", by.y = "gene_id", all.x = T)
exprPlusGenes <- exprPlusGenes[rev(order(exprPlusGenes$family.y)),]
{
    meanplot <- ggplot(exprPlusGenes, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = Function.y), shape=21, color = "None") +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  geom_label_repel(data = exprPlusGenes[which(abs(exprPlusGenes$PC2) >= 0.05 | abs(exprPlusGenes$PC4) >= 0.05),], aes(label = name.y, fill = Function.y), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
  xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",pcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,2), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}
ggsave("~/R/Eutrema/PS/pics/PSMeanPC24phossulf.pdf", meanplot, dpi = 150, height = 12, width = 14)
###################################################
lncRNA tree

psMean <- as.data.frame(ps_mean, StringsAsFactors = F)
psPred <- psMean[psMean$ps != 0,]
psPred <- rownames(psPred[which(toupper(rownames(psPred)) %in% toupper(allGlnc)),])
PsPred <- psMean[psMean$ps != 0,]
PsPred <- rownames(PsPred[which(toupper(rownames(PsPred)) %in% toupper(allGlnc)),])
pSPred <- psMean[psMean$ps != 0,]
pSPred <- rownames(pSPred[which(toupper(rownames(pSPred)) %in% toupper(allGlnc)),])
PSPred <- psMean[psMean$ps != 0,]
PSPred <- rownames(PSPred[which(toupper(rownames(PSPred)) %in% toupper(allGlnc)),])

predSet <- list(ps=psPred, pS=pSPred, Ps=PsPred, PS=PSPred)

pdf("pics/lncVenn.pdf", width = 5, height = 5)
venn(predSet, zcolor = "style")
dev.off()


# DEG lncRNA
# lists = r-binded

lncPs <- res.Ps[!is.na(res.Ps$padj),]
lncPs <- lncPs[lncPs$padj <=0.05,]
lncPs <- lncPs[which(toupper(rownames(lncPs)) %in% toupper(allGlnc)),]
lncpS <- res.pS[!is.na(res.Ps$padj),]
lncpS <- lncpS[lncpS$padj <=0.05,]
lncPs <- lncPs[which(toupper(rownames(lncPs)) %in% toupper(allGlnc)),]
lncps <- res.ps[!is.na(res.Ps$padj),]
lncps <- lncps[lncps$padj <=0.05,]
lncPs <- lncPs[which(toupper(rownames(lncPs)) %in% toupper(allGlnc)),]

lncDEG <- list(lncPs=lncPs, lncpS=lncpS, lncps=lncps)

#################################################################################
# arcsinh transformation

# results <- results(dds)
countsNorm <- counts(dds, normalized = T)
colnames(countsNorm)  <- sampNames
countasinh <- asinh(countsNorm)
pcaall <- prcomp(countasinh, center  = T, scale = T)

aexprVals <- data.frame(pcaall$x)
asampleVals  <-  data.frame(pcaall$rotation)

singleTreatment <- paste0(rep(c("ps", "pS", "Ps", "PS"), each = 3), 1:3)
singleTreatment <- factor(singleTreatment, levels = singleTreatment)
apcaSum <- as.data.frame(summary(pcaall)$importance)
acoords <- data.frame(asampleVals, Treatment = singleTreatment, X=rep(0,12), Y=rep(0,12))


# combining the biological replicates
aps <- countasinh[,1:3]
aPs <- countasinh[,4:6]
apS <- countasinh[,7:9]
aPS <- countasinh[,10:12]

apS <- (pS=rowMeans(apS))
aPs <- (Ps=rowMeans(aPs))
aps <- (ps=rowMeans(aps))
aPS <- (PS=rowMeans(aPS))

ps_mean <- cbind(aps, apS, aPs, aPS)

library(grid)
library(gridExtra)
library(cowplot)
pcaCount <- prcomp(ps_mean, center=TRUE, scale=TRUE)  # PCA with centering and scaling
countRot <- data.frame(pcaCount$rotation)  # The loadings are here
plotlists <- {
PCplots <- list()
count <- 1
for (a in 1:4) {
    for (b in 1:4) {
        if (a < b) {
            plots <- ggplot(aexprVals, aes_string(paste0("PC",a), paste0("PC", b))) +
    	    geom_point(shape=19, color = "black", alpha = 0.2, fill = "None") +
            geom_segment(data=acoords, aes_string(x="X", y="Y", xend=paste0("PC",a), yend=paste0("PC", b), color="Treatment"), 
                         arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
            scale_color_manual(values=rep(c( "#febfcb", "gold", "#f91301", "#fd8a19" ), each=3), name = "") +
  #   geom_label_repel(data = aexprVals[which(abs(aexprVals$PC3) >= 0.1 | abs(aexprVals$PC4) >= 0.1),], aes(label =name, fill = family), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
            labs(title = paste0("PC", a, " vs. PC", b)) +
            xlab(paste0("PC", a, " ", apcaSum[2,a]*100,"%")) + 
            ylab(paste0("PC", b, " ", apcaSum[2,b]*100,"%")) +
            coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1,1)) +
            theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
                  axis.text=element_text(size=12)) +
	    theme_bw()
            PCplots[[count]] <- plots
            count <- count + 1

        }
    }
}
}

plotted <- plot_grid(plotlist = PCplots,  ncol = 3, align =  "v")
plotted
ggsave("pics/arcsinhPCmeanComparison.pdf", plotted, dpi = 250, height = 14, width = 22)


screeplot(pcaCount)
#log.expr <- log2(expr+1)

plot(pcaCount, type = "l")
summary(pcaCount)
screeplot(pcaCount)
cexprVals<-data.frame(pcaCount$x)
csampleVals<-data.frame(pcaCount$rotation)
cpcaSum <- as.data.frame(summary(pcaCount)$importance)
treat <- c("ps", "pS", "Ps", "PS")
# singeTreatment <- factor(singeTreatment, levels = singeTreatment)
ccoords <- data.frame(X=rep(0,4), Y=rep(0,4), csampleVals, Treatment=treat)

cexprVals <- merge(cexprVals, lipidGenes[,c(1,8:10)], by.x = 0, by.y = "Eutr", all.x = T)

cexprVals <- cexprVals[rev(order(cexprVals$family)),]
{meanplot <- ggplot(cexprVals, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = family), shape=21, color = "None") +
    geom_segment(data=ccoords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  geom_label_repel(data = cexprVals[which(abs(cexprVals$PC2) >= 0.1 | abs(cexprVals$PC4) >= 0.1),], aes(label =name, fill = family), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
    xlab(paste0("PC2 ", cpcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",cpcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1,1))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}
ggsave("~/R/Eutrema/PS/pics/asinhPSMeanPC24lipids.pdf", meanplot, dpi = 150, height = 12, width = 14)  # 


exprGS <- cexprVals
exprGS$Type <- ifelse(exprGS$Row.names %in% GSL$Eutr, "GSL", ifelse(exprGS$Row.names %in% Sassim$Eutr, "S assimilation", "other"))
exprGS$Type <- factor(exprGS$Type, levels = c("other", "GSL", "S assimilation"))
rbindGS <- rbind(GSL, Sassim)
exprGS <- merge(exprGS, rbindGS, all.x = T, by.x = "Row.names", by.y = "Eutr")
# rownames(exprGS) <- exprGS$Row.names

exprGS <- exprGS[order(exprGS$Type),]

# lncRNA stuff
{meanplot <- ggplot(exprGS, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = Type), shape=21, color = "None") +
    scale_fill_manual(values = c( "gray", "coral3", "chartreuse4")) +
    geom_segment(data=ccoords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    geom_label_repel(data = exprGS[which(abs(exprGS$PC2) >= 0.1 | abs(exprGS$PC4) >= 0.1),], aes(label = gene_name, fill = Type), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
    xlab(paste0("PC2 ", cpcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ", cpcaSum$PC4[2]*100,"%")) +
    #     scale_fill_manual(values = c( "gray", "coral3", "chartreuse4")) +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1,1))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}

ggsave("~/R/Eutrema/PS/pics/asinhPSMeanPC24GSL.pdf", meanplot, dpi = 150, height = 12, width = 14)  # 

theGenes <- read.table("~/R/Eutrema/PS/Genes", sep = "\t")
colnames(theGenes) <- c("family", "name", "gene_id", "Function")

cexprPlusGenes <- merge(cexprVals, theGenes, by.x = "Row.names", by.y = "gene_id", all.x = T)
cexprPlusGenes <- cexprPlusGenes[rev(order(cexprPlusGenes$family.y)),]
{
    meanplot <- ggplot(cexprPlusGenes, aes(x = PC2, y = PC4)) +
    geom_point(aes(fill = Function), shape=21, color = "None") +
    geom_segment(data=ccoords, aes(x=X, y=Y, xend=PC2, yend=PC4, color=Treatment), 
                 arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
  scale_color_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  geom_label_repel(data = cexprPlusGenes[which(abs(cexprPlusGenes$PC2) >= 0.05 | abs(cexprPlusGenes$PC4) >= 0.05),], aes(label = name.y, fill = Function), segment.alpha = 0.7, color = "white", alpha = 0.8, segment.color = "black") +
  xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",pcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-2,2), ylim=c(-1.5,1.5))+
theme(axis.title=element_text(size=15), legend.text=element_text(size=12),
axis.text=element_text(size=12))
meanplot}
ggsave("~/R/Eutrema/PS/pics/asinhPSMeanPC24phossulf.pdf", meanplot, dpi = 150, height = 12, width = 14)

cps <- (cps=rowMeans(countsNorm[,1:3]))
cPs <- (cPs=rowMeans(countsNorm[,4:6]))
cpS <- (cpS=rowMeans(countsNorm[,7:9]))
cPS <- (cPS=rowMeans(countsNorm[,10:12]))

cps <- sort(cps)
cPs <- sort(cPs)
cpS <- sort(cpS)
cPS <- sort(cPS)
ps <- sort(ps)
Ps <- sort(Ps)
pS <- sort(pS)
PS <- sort(PS)
aps <- sort(aps)
aPs <- sort(aPs)
apS <- sort(apS)
aPS <- sort(aPS)

expData <- data.frame(x=1:length(PS),ps=ps, Ps=Ps, pS=pS, PS=PS)
aexpData <- data.frame(x=1:length(aPS), aps=aps, aPs=aPs, apS=apS, aPS=aPS)
cexpData <- data.frame(x=1:length(cPS), cps=cps, cPs=cPs, cpS=cpS, cPS=cPS)

library(reshape2)
expMelt <- melt(expData, id.vars = "x")
aexpMelt <- melt(aexpData, id.vars = "x")
cexpMelt <- melt(cexpData, id.vars = "x")
{Tnorm <- cexpMelt %>% ggplot(., aes(x = x, y = value, group = variable, color = variable)) + 
    geom_line() + theme_bw() +
    labs(x = "sort order", y = "Normalized read counts", title = "Normalization")
}
{Tfpkm <- expMelt %>% ggplot(., aes(x = x, y = value, group = variable, color = variable)) + 
    geom_line() + theme_bw() +
    labs(x = "sort order", y = "log2(fpkm + 1)", title = "Logarithmic FPKM")
}
{Tasin <- aexpMelt %>% ggplot(., aes(x = x, y = value, group = variable, color = variable)) + 
    geom_line() + theme_bw() +
    labs(x = "sort order", y = "asinh(CountNorm)", title = "Hyperbolic arcsine")
}
transformPlots <- plot_grid(plotlist = list(Tnorm, Tfpkm, Tasin), nrow = 1, align = "h")
ggsave("~/R/Eutrema/PS/pics/transformPlots.pdf", transformPlots, dpi = 250, height = 5, width = 18)

