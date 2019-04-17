library(DESeq2) #DEG analysis
library(plyr) #data manipulation
library(gplots) #graphing and visualization
library(RColorBrewer) #graphing and visualization
# library(stringr) #idk
library(genefilter) #GFF maker
library(dplyr) #graphing
library(geneplotter) #graphing
library(ggplot2) #graphing
library(GenomicFeatures) #GFF functions makeTxdbFromGff
library(reshape2) #graphing
library(ashr) #DESeq2
# library(gridExtra) #used for Venn and graphs
# library(cowplot)  #used for log2 expression graphs
# library(venn) #used for venn
# library(VennDiagram) #used for venn
# library(grid) #used for venn and graphs
# library(topGO)  #used for GO enrichment (results not informative)
# library(Rgraphviz)  #used for GO enrichment
# library(pheatmap) #used for visualization of DESeq outputs
setwd("~/scratch/PS") #stuff in here

#####################SETUP###################
# condition <- c("F", "F", "F", "Sh", "Sh", rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3))
condition <- c("F", "F", "F", rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3))
sampleNO <- c("Y-F2003-1",
              "Y-F2003-2",
              "cracker",
              # "SWW-1",
              # "SWW-2",
              "ps10_S1",
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

#specifying the reference level
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPHS")
# DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "Sh")

dds <- DESeq(DESeq2Table, test = "Wald")
res <- results(dds, pAdjustMethod = "BH")
res
rld <- rlogTransformation(dds, blind = F)
# png("PCA_samplesSh.png", height = 4, width = 6, units = "in", res = 300)
DESeq2::plotPCA(rld, intgroup=c("condition"))
# dev.off()
#how many genes are significantly expressed
table(res$padj < 0.05)

resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_HPHS"    "condition_FCC_vs_HPHS"  "condition_HPLS_vs_HPHS" "condition_LPHS_vs_HPHS" "condition_LPLS_vs_HPHS"

# res1_2 <- lfcShrink(dds, coef = "condition_HPLS_vs_Sh", type = "ashr")
# res1_3 <- lfcShrink(dds, coef = "condition_LPHS_vs_Sh", type = "ashr")
# res1_4 <- lfcShrink(dds, coef = "condition_LPLS_vs_Sh", type = "ashr")
# res1_I <- lfcShrink(dds, coef = "Intercept", type = "ashr")
# res1_F <- lfcShrink(dds, coef = "condition_F_vs_Sh", type = "ashr")
# res1_1 <- lfcShrink(dds, coef = "condition_HPHS_vs_Sh", type = "ashr")

res1_2 <- lfcShrink(dds, coef = "condition_HPLS_vs_HPHS", type = "ashr")
res1_3 <- lfcShrink(dds, coef = "condition_LPHS_vs_HPHS", type = "ashr")
res1_4 <- lfcShrink(dds, coef = "condition_LPLS_vs_HPHS", type = "ashr")
res1_I <- lfcShrink(dds, coef = "Intercept", type = "ashr")
res1_F <- lfcShrink(dds, coef = "condition_F_vs_HPHS", type = "ashr")

#########################PLOTTING###############################################
#r1_1 <- data.frame(cbind(row.names(res1_1), rep(1, nrow(res1_1)), res1_1$lfcSE, rep("HPHS", nrow(res1_1))), stringsAsFactors = FALSE)
r1_I <- cbind.data.frame(row.names(res1_I), rep(0, nrow(res1_I)), res1_I$lfcSE, rep("HPHS", nrow(res1_I)))
colnames(r1_I) <- c("gene", "lfc", "se", "sample")
# r1_I <- cbind.data.frame(row.names(res1_I), rep(0, nrow(res1_1)), res1_1$lfcSE, rep("Sh", nrow(res1_I)))
# colnames(r1_I) <- c("gene", "lfc", "se", "sample")
# r1_1 <- cbind.data.frame(row.names(res1_1), res1_2$log2FoldChange, res1_2$lfcSE, rep("HPHS", nrow(res1_1)))
# colnames(r1_1) <- c("gene", "lfc", "se", "sample")
r1_2 <- cbind.data.frame(row.names(res1_2), res1_2$log2FoldChange, res1_2$lfcSE, rep("HPLS", nrow(res1_2)))
colnames(r1_2) <- c("gene", "lfc", "se", "sample")
r1_3 <- cbind.data.frame(row.names(res1_3), res1_3$log2FoldChange, res1_3$lfcSE, rep("LPHS", nrow(res1_3)))
colnames(r1_3) <- c("gene", "lfc", "se", "sample")
r1_4 <- cbind.data.frame(row.names(res1_4), res1_4$log2FoldChange, res1_4$lfcSE, rep("LPLS", nrow(res1_4)))
colnames(r1_4) <- c("gene", "lfc", "se", "sample")
r1_F <- cbind.data.frame(row.names(res1_F), res1_F$log2FoldChange, res1_F$lfcSE, rep("F", nrow(res1_F)))
colnames(r1_F) <- c("gene", "lfc", "se", "sample")
# r1_FCC <- cbind.data.frame(row.names(res1_FCC), res1_FCC$log2FoldChange, res1_FCC$lfcSE, rep("FCC", nrow(res1_FCC)))
# colnames(r1_FCC) <- c("gene", "lfc", "se", "sample")
lists <- rbind(r1_2, r1_3, r1_4, r1_F, r1_I)
# lists <- rbind(r1_1, r1_2, r1_3, r1_4, r1_F, r1_I)
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

 graph <- function(ID) {
  dat <- lists[lists$gene %in% ID, ]
  ggplot(dat, aes(x = sample, y = lfc)) +
    geom_point() +
    xlab(paste0(as.character(dat$gene[1]), " expression")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(x = sample, ymax = lfc + se, ymin = lfc - se), width = 0.1) +
    ylim(-5,5) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank())
}


################################TREES###################################

png("fiveVennvsS.png", width = 1000, height = 1000, res = 200)
venn(list(HPLS=rownames(res1_2[which(res1_2[,"padj"] < 0.05),]),
LPHS=rownames(res1_3[which(res1_3[,"padj"] < 0.05),]),
LPLS=rownames(res1_4[which(res1_4[,"padj"] < 0.05),]),
F=rownames(res1_F[which(res1_F[,"padj"] < 0.05),]),
HPHS=rownames(res1_1[which(res1_1[,"padj"] < 0.05),])), zcolor = "style")
dev.off()


png("smolvenn.png", width = 800, height = 800, res = 250)
venn(list(HPLS=rownames(res1_2[which(res1_2[,"padj"] < 0.05),]),
                       LPHS=rownames(res1_3[which(res1_3[,"padj"] < 0.05),]),
                       LPLS=rownames(res1_4[which(res1_4[,"padj"] < 0.05),])), zcolor = "style")
dev.off()


png("fieldvenn.png", width = 800, height = 800, res = 250)
venn(list(F2003=rownames(res1_F[which(res1_F[,"padj"] < 0.05),]),
FCC=rownames(res1_FCC[which(res1_FCC[,"padj"] < 0.05),])), zcolor = "style", family = "mono")
dev.off()

print(fieldvenn)
print(venn)
print(smallvenn)


upPS <- venn.diagram(list(HPLS=rownames(res1_2[which(res1_2[,"padj"] < 0.05 & res1_2[,"log2FoldChange"] > 0),]),
                  LPHS=rownames(res1_3[which(res1_3[,"padj"] < 0.05 & res1_3[,"log2FoldChange"] > 0),]),
                  LPLS=rownames(res1_4[which(res1_4[,"padj"] < 0.05 & res1_4[,"log2FoldChange"] > 0),])), filename = NULL, col = "red", rotation.degree = 100)
downPS <-venn.diagram(list(HPLS=rownames(res1_2[which(res1_2[,"padj"] < 0.05 & res1_2[,"log2FoldChange"] < 0),]),
                   LPHS=rownames(res1_3[which(res1_3[,"padj"] < 0.05 & res1_3[,"log2FoldChange"] < 0),]),
                   LPLS=rownames(res1_4[which(res1_4[,"padj"] < 0.05 & res1_4[,"log2FoldChange"] < 0),])), filename = NULL, col = "navyblue")
upF <-venn.diagram(list(F2003=rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] > 0),]),
                FCC=rownames(res1_FCC[which(res1_FCC[,"padj"] < 0.05 & res1_FCC[,"log2FoldChange"] > 0),])), filename = NULL, col = "red")

downF <- venn.diagram(list(F2003=rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] < 0),]),
                     FCC=rownames(res1_FCC[which(res1_FCC[,"padj"] < 0.05 & res1_FCC[,"log2FoldChange"] < 0),])), filename = NULL, col = "navyblue")



table <- rbind(HPLS=c(52,23), LPHS=c(1,0), LPLS=c(10, 9), F2003=c(3933,3477), FCC=c(2463,3289))
colnames(table) <- c("Upreg", "Downreg")
table <- as_tibble(table)
table <- mutate(table, Total=(Upreg + Downreg), Sample =c( "HPLS", "LPHS","LPLS", "F", "FCC") )

table

write.table(table, file = "expCounts.tab", quote = F, row.names = F)

png("vennDiagram.png", height = 800, width = 800, res = 125)
grid.arrange(gTree(children = gList(textGrob(".", gp=gpar(col = "white", fontfamily = "serif")))),
  gTree(children = gList(textGrob("Upregulated", gp = gpar(fontfamily = "serif")))),
  gTree(children = gList(textGrob("Downregulated", gp = gpar(fontfamily = "serif")))),
gTree(children = gList(textGrob("PS_LIB", rot = 90, gp = gpar(fontfamily = "serif")))),
gTree(children = upPS),
gTree(children = downPS),
gTree(children = gList(textGrob("Field", rot = 90, gp = gpar(fontfamily = "serif")))),
gTree(children = upF),
gTree(children = downF),  ncol = 3, widths = c(0.1, 1, 1), heights = c(0.1, 1, 1))
dev.off()



names <- rownames(res[which(res[,"padj"] < 0.05),])
names
thing <- sapply(names, graph, USE.NAMES = F, simplify = F)
plotted <- plot_grid(plotlist = thing,  ncol = 4, align =  "v")
#plotted <- plot_grid(plotted, legend, ncol = 2, rel_widths = c(1, 0.06))
plotted
plotted <- grid.arrange(arrangeGrob(plotted, left = textGrob("log2-fold Change", rot = 90, vjust = 1)))

ggsave("Thing.png", plotted, height = 9, width = 12, dpi = 250)

genes <- read.table("Genes")
genes2 <- genes[genes[,2] %in% lists$gene, ]
oldNames <- genes2[,2]

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
    ylim(min(dat$lfc-dat$se-1),max(dat$lfc+dat$se+1)) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank())
}

thing <- sapply(oldNames, graph, USE.NAMES = F, simplify = F)
plotted <- plot_grid(plotlist = thing,  ncol = 4, align =  "v")
plotted
plotted <- grid.arrange(arrangeGrob(plotted, left = textGrob("log2-fold Change", rot = 90, vjust = 1)))
ggsave("BeforeDot.png", plotted, height = 12, width = 10, dpi = 150)

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
GOres1_FCC <- as.data.frame(res1_FCC@listData)
rownames(GOres1_FCC) <- res1_FCC@rownames

dat <- list(HPLS=GOres1_2, LPHS=GOres1_3, LPLS=GOres1_4, F2003=GOres1_F, FCC=GOres1_FCC)
names <- c("HPLS", "LPHS", "LPLS", "F2003", "FCC")

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
                   description = paste0("HPHS vs. ", name),
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
  # png(paste0(name, "GOPos.png"), height = 1500, width = 1500, res = 200)
  # showSigOfNodes(pos_genes, score(resFtest), firstSigNodes = 10, useInfo = "all")
  # dev.off()

  #Down regulated genes
  Fdown <- data[data$log2FoldChange < 0,]$padj
  names(Fdown) <- rownames(data[data$log2FoldChange < 0, ])
  Fdown <- Fdown[complete.cases(Fdown)]

  neg_genes <- new("topGOdata",
                   description = paste0("HPHS vs. ", name),
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

  # png(paste0(name, "GONeg.png"), height = 1500, width = 1500, res = 200)
  # showSigOfNodes(neg_genes, score(resFtestN), firstSigNodes = 10, useInfo = "all")
  # dev.off()
  i <- i + 1
}

name <- "LPLS"
data <- GOres1_4
  Fup <- data[data$log2FoldChange > 0,]$padj
  names(Fup) <- rownames(data[data$log2FoldChange > 0, ])
  Fup <- Fup[complete.cases(Fup)]

  pos_genes <- new("topGOdata",
                   description = paste0("HPHS vs. ", name),
                   ontology = "BP",
                   allGenes = Fup,
                   geneSel = gene_filter,
                   nodeSize = 10,
                   annotationFun = annFUN.gene2GO,
                   gene2GO = GOs)

  pos <- topGO::genes(pos_genes)

  sig <- sigGenes(pos_genes)

  # test <- new("classicCount", testStatistic = GOFisherTest, name = "test",
  #             allMembers = pos)

  testFisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher")

  resFtest <- getSigGroups(pos_genes, testFisher)

  posTable <- GenTable(pos_genes, Fisher = resFtest, topNodes = 200)
  write.table(posTable, file = paste0(name, "_POSlong.tab"), sep = "\t", quote = F)
  png(paste0(name, "GOPos.png"), height = 1500, width = 1500, res = 200)
  showSigOfNodes(pos_genes, score(resFtest), firstSigNodes = 10, useInfo = "all")
  dev.off()

###########################GET LIST OF GENE NAMES #################


rownames(res1_2[which(res1_2[,"padj"] < 0.05 & res1_2[,"log2FoldChange"] > 0),])
rownames(res1_2[which(res1_2[,"padj"] < 0.05 & res1_2[,"log2FoldChange"] < 0),])
rownames(res1_3[which(res1_3[,"padj"] < 0.05 & res1_3[,"log2FoldChange"] > 0),])
rownames(res1_3[which(res1_3[,"padj"] < 0.05 & res1_3[,"log2FoldChange"] < 0),])
rownames(res1_4[which(res1_4[,"padj"] < 0.05 & res1_4[,"log2FoldChange"] > 0),])
rownames(res1_4[which(res1_4[,"padj"] < 0.05 & res1_4[,"log2FoldChange"] < 0),])
rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] > 0),])
rownames(res1_F[which(res1_F[,"padj"] < 0.05 & res1_F[,"log2FoldChange"] < 0),])
# rownames(res1_FCC[which(res1_FCC[,"padj"] < 0.05 & res1_FCC[,"log2FoldChange"] > 0),])
# rownames(res1_FCC[which(res1_FCC[,"padj"] < 0.05 & res1_FCC[,"log2FoldChange"] < 0),])

liHPLS <- (res1_2[which(res1_2[,"padj"] < 0.05),])
liLPHS <- (res1_3[which(res1_3[,"padj"] < 0.05),])
liLPLS <- (res1_4[which(res1_4[,"padj"] < 0.05),])
liF <- (res1_F[which(res1_F[,"padj"] < 0.05),])
# liFCC <- (res1_FCC[which(res1_FCC[,"padj"] < 0.05),])
# liHPLS
# liLPHS
# liLPLS
# liF
# liFCC
# writeLines(liHPLS, con = "HPLSnames")

#############################FPKM#######################


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
  inter2 <- intersect(rownames(lengths), rownames(dds))
  length(inter2)
  lens <- lengths[row.names(lengths) %in% inter2,]
  mcols(dds)$basepairs <- lens
    fpkms2 <- fpkm(dds, robust = T)
#GTF
# inter <- intersect(myGTF$geneID, rownames(DESeq2Table))
# head(inter)
# smolGTF <- myGTF %>% filter( geneID %in% inter )

#GFF

lists <- list(HPLS=res1_2, LPHS=res1_3, LPLS=res1_4, Field=res1_F)

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
###############PHEATMAP#############
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:6]
df <- as.data.frame(colData(dds)[,c("condition")])
#must have this or else have error msg: Error in check.length("fill") :
#                                         'gpar' element 'fill' must not be length 0
rownames(df) <- colnames(assay(ntd)[select,])
pheatmap::pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



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

bigmat <- assay(rld)
head(bigmat)

#Extracting only positive -> Bray analysis of negative values are trivial
bigmat[bigmat < 0] <- 0
logcount <- bigmat
head(logcount)
keep2 <- rowSums(logcount) >= 10
logcount <- logcount[keep2,]
nrow(logcount)

subsample <- logcount[sample(nrow(logcount), 2000),]
head(as.matrix(subsample))
bray1 <- vegdist(as.matrix(subsample), "bray")


head(bray1)
biplot(ape::pcoa(bray1))
# countdist <- phyloseq::distance(logcount, methods = "bray") doesn't work
# countdist <- vegdist(logcount, "bray") don't do this

bigmat2 <- assay(vst(dds))
head(bigmat2)
# matdist <- dist(bigmat, method = "minkowski") <- never do this
matdist <- vegdist(bigmat, "bray")
matdist2 <- vegdist(bigmat2, "bray")
# poca(matdist) <- never do this either
max(matdist)
max(matdist2)
min(matdist2)
hist(matdist2)

raw <- counts(dds)
head(raw)

cas <- CA(raw, graph = T)i

vst <- varianceStabilizingTransformation(dds)
vst <- assay(vst)
###############################################

#try a subset of data (2000 genes that don't have 0 or negative expression, randomly sampled...)
submat <- otu_table(subsample, taxa_are_rows = T)
taxmat <- matrix(rep(condition, each = nrow(submat)), nrow = nrow(submat))
colnames(taxmat) <- colnames(submat)
rownames(taxmat) <- rownames(submat)
taxmat <- tax_table(taxmat)
physeq <- phyloseq(submat)
sampledat <- sample_data(data.frame(row.names=sample_names(physeq), condition = condition, stringsAsFactors = F,
                                    sampleID = sampleNO))
physeq <- merge_phyloseq(physeq, sampledat)
# plot_bar(physeq)
physeq
ord <- ordinate(physeq, method = "PCoA", distance = "bray")
plot_ordination(physeq, ord, type = "samples", color = "condition")


#########big###
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


###Subsetteds dds

vst <- varianceStabilizingTransformation(dds)
subsetvst <- vst[, vst$condition %in% c("HPLS", "LPHS", "LPLS", "HPHS")]
subsetvst <- assay(subsetvst)
smallOtu <- otu_table(subsetvst, taxa_are_rows = T)
sampleSmall <- c("ps10_S1","ps11_S2", "ps12_S3",  "ps37_S4",   "ps38_S5",  "ps40_S6",
              "ps41_S7",   "ps44_S8",  "ps46_S9",  "ps48_S10","ps49_S11", "ps50_S12")

sampleDatS <- sample_data(data.frame(row.names=sample_names(smallOtu), condition = rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3),
                                     sampleID =sampleSmall, stringsAsFactors = F))
phySmall <- phyloseq(smallOtu, sampleDatS)
ordS <- ordinate(phySmall, method = "PCoA")
pocaS <- plot_ordination(phySmall, ordS, type = "samples", color = "condition")
pocaS
ggsave("vst-pcoaS.png", pocaS, height = 4, width = 6, dpi = 125)
