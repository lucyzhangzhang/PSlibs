setwd("~/R/Eutrema/PS/")
library(WGCNA)
library(reshape2)
library(ggdendro)
library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork) # do I really need this it's a ggplot api wrapper apparently
# yes I do need this in fact it allows me to combine two plots into one
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

fpkm.raw <- read.table("FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)

gsg <- goodSamplesGenes(fpkm.raw, verbose = 3)
gsg$allOK

fpkm <- t(fpkm.raw)
sampleTree <- hclust(dist(fpkm), method = "average")

# choose the power to use for clustering
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# already did it (chose 13)
# this takes a long time so don't do it unless you're 
# doing it for the first time (and don't forget to record the number)
sft <- pickSoftThreshold(fpkm, powerVector = powers, verbose = 3, 
                         networkType = "signed hybrid")
sizeGrWindow(9, 5)
par(mfrow = c(1, 2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               xlab = "Soft Threshold (power)",
                               ylab = "Sclae Free Topology Model Fit, signed R^2",
                               type = "n", main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red");
abline(h = 0.90, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

# making the co-expression network
network <- blockwiseModules(fpkm, power = 13,
                       TOMType = "none",
                       networkType = "signed hybrid",
                       randomSeed = 319, corType = "pearson", 
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       minCoreKME = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "droughtTOM",
                       verbose = 3, maxBlockSize = 10000)

clustcolours <- labels2colors(network$colors)
MEs <- moduleEigengenes(fpkm, clustcolours)$eigengenes

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(fpkm, clustcolours)$eigengenes
MEs <- orderMEs(MEs0)

# group by treatment
treatment <- c(rep("ps", 3), 
               rep("Ps", 3),
               rep("pS", 3),
               rep("PS", 3))
                   
treatment <- data.frame(treatment = as.factor(treatment))
rownames(treatment) <- colnames(fpkm.raw)
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))
treat <- treat[,-1]
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat));

# clustering the eigengenes
dissim <- (1-t(cor(MEs0, method = "p")))/2
hclustME <- hclust(as.dist(dissim), method = "average")

# visualizing clusters
par(mfrow = c(1, 1))
plot(hclustME, main = "Clustering tree based off the module eigengenes") 
sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise";
plotMat(t(scale(fpkm[,clustcolours==which.module ]) ),nrgcols=30,rlabels=T,
clabels=T,rcols=which.module,
title=which.module )

# for the second (blue) module we use
which.module="blue";
plotMat(t(scale(fpkm[,clustcolours==which.module ]) ),nrgcols=30,rlabels=T,
clabels=T,rcols=which.module,
title=which.module )
which.module="brown";
plotMat(t(scale(fpkm[,clustcolours==which.module ]) ),nrgcols=30,rlabels=T,
clabels=T,rcols=which.module,
title=which.module )

# relationship between clust module and module eigengenes
sizeGrWindow(8,7);
which.module="green"
ME=MEs0[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(fpkm[,clustcolours==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="array sample")

#tabularize cluster and order by Frequency
geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames

#plotting
library(gplots)
library(RColorBrewer)

colnames(moduleTraitCor) <- c("ps", "Ps", "pS", "PS")

sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))
moduleTraitCor[moduleTraitPvalue > 0.05] <- NA
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
head(moduleTraitCor)
corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
correlation_ggplot <- corr_clust%>% melt()

dx <- dendro_data(sample_clust)
dy <- dendro_data(module_clust)
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()

{heatmap <- ggplot(data = correlation_ggplot, 
                   aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#16316A", mid = "white", 
        high = "#61143F", na.value = "white", name = "Correlation", 
        limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  labs(x="Treatment", y="Clusters") +  theme_bw() +     
  theme(axis.text.y=element_text(angle=20,vjust=0.5), 
        plot.margin = margin(-0.75, 0, 0, 0, "cm"), legend.position = "bottom")
heatmap <- px + heatmap + plot_layout(ncol = 1, heights = c(1, 10))
 heatmap}
ggsave("wgcnaHeatmap.pdf", heatmap, dpi = 300, height = 12, width = 6)

##########SEPARATE TREATMENTS##################
                   
treat2 <- data.frame(as.factor(colnames(fpkm.raw)))
rownames(treat2) <- colnames(fpkm.raw)
treat22 <- model.matrix(~., data = treat2,
                      contrasts.arg = lapply(treat2, contrasts, contrasts = F))
treat22 <- treat22[,-1]
# testing for noise
moduleTraitCor <- WGCNA::cor(MEs, treat22, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat22));

moduleTraitPvalueFisher <- corPvalueFisher(moduleTraitCor,length(treat22));


#tabularize cluster and order by Frequency
geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames

viewcolors <- data.frame(table(clustcolours))
viewcolors <- viewcolors[order(-viewcolors$Freq),]

colnames(moduleTraitCor) <- colnames(fpkm.raw)

sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))
moduleTraitCor[moduleTraitPvalue> 0.05] <- NA
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
head(moduleTraitCor)
corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
correlation_ggplot <- corr_clust%>% melt()

dx <- dendro_data(sample_clust)
dy <- dendro_data(module_clust)

# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()

{heatmapT <- ggplot(data = correlation_ggplot, 
                   aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#16316A", mid = "white", 
        high = "#61143F", na.value = "white", name = "Correlation", 
        limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  labs(x="Treatment", y="Clusters") +  theme_bw() +     
  theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(-0.75, 0, 0, 0, "cm"), legend.position = "bottom")
heatmapT <- px + heatmapT + plot_layout(ncol = 1, heights = c(1, 10))
          heatmapT}
ggsave("wgcnaHeatmapall.pdf", heatmapT, dpi = 300, height = 12, width = 8)
heatmapT
heatmap2

# remove PS1 from the analysis
ps1Remove.raw <- fpkm.raw[,-1]

fpkm <- t(ps1Remove.raw)
sampleTree <- hclust(dist(fpkm), method = "average")

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

sft <- pickSoftThreshold(fpkm, powerVector = powers, verbose = 5, 
                         networkType = "signed hybrid")
sizeGrWindow(9, 5)
par(mfrow = c(1, 2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               xlab = "Soft Threshold (power)",
                               ylab = "Sclae Free Topology Model Fit, signed R^2",
                               type = "n", main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red");
abline(h = 0.90, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

# making the co-expression network
network <- blockwiseModules(fpkm, power = 16,
                       TOMType = "none",
                       networkType = "signed hybrid",
                       randomSeed = 319, corType = "pearson", 
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "droughtTOM",
                       verbose = 3, maxBlockSize = 10000)
clustcolours <- labels2colors(network$colors)
MEs <- moduleEigengenes(fpkm, clustcolours)$eigengenes
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(fpkm, clustcolours)$eigengenes
MEs <- orderMEs(MEs0)


#removal of the ps1 sample
treatment <- c(rep("ps", 2), 
               rep("Ps", 3),
               rep("pS", 3),
               rep("PS", 3))
                   
treatment <- data.frame(treatment = as.factor(treatment))
rownames(treatment) <- colnames(fpkm.raw)
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))
treat <- treat[,-1]
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat));

#tabularize cluster and order by Frequency
geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames

viewcolors <- data.frame(table(clustcolours))
viewcolors <- viewcolors[order(-viewcolors$Freq),]
#plotting

colnames(moduleTraitCor) <- c("ps", "Ps", "pS", "PS")

sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))
moduleTraitCor[moduleTraitPvalue > 0.05] <- NA
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
#rownames(moduleTraitCor) <- rep("", nrow(moduleTraitCor))
head(moduleTraitCor)

corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
# corr_clust <- corr_clust[ordr(rownames(viewcolors)),]
correlation_ggplot <- corr_clust%>% melt()

dx <- dendro_data(sample_clust)
dy <- dendro_data(module_clust)
# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()

{heatmap3 <- ggplot(data = correlation_ggplot, 
                   aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#16316A", mid = "white", 
        high = "#61143F", na.value = "white", name = "Correlation", 
        limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  labs(x="Treatment", y="Clusters") +  theme_bw() +     theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(-0.75, 0, 0, 0, "cm"), legend.position = "bottom")
heatmap3 <- px + heatmap3 + plot_layout(ncol = 1, heights = c(1, 10))
 heatmap3}
heatmap3
ggsave("wgcnaHeatmapsPS1removed.pdf", heatmap3, dpi = 300, height = 12, width = 6)

###########P and S independently###########
# S only
treatment <- c(rep("s", 3), 
               rep("s", 3),
               rep("S", 3),
               rep("S", 3))
                   
treatment <- data.frame(treatment = as.factor(treatment))
rownames(treatment) <- colnames(fpkm.raw)
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))
treat <- treat[,-1]
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat));

#tabularize cluster and order by Frequency
geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames

viewcolors <- data.frame(table(clustcolours))
viewcolors <- viewcolors[order(-viewcolors$Freq),]

colnames(moduleTraitCor) <- c("s","S")

moduleS <- moduleTraitCor

# P only
treatment <- c(rep("p", 3), 
               rep("P", 3),
               rep("p", 3),
               rep("P", 3))

treatment <- data.frame(treatment = as.factor(treatment))
rownames(treatment) <- colnames(fpkm.raw)
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))
treat <- treat[,-1]
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat));

#tabularize cluster and order by Frequency
geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames

colnames(moduleTraitCor) <- c("p","P")

moduleP <- moduleTraitCor

# moduleTraitCor <- cbind(moduleP, moduleS)

sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))
moduleTraitCor[moduleTraitPvalue > 0.05] <- NA
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
#rownames(moduleTraitCor) <- rep("", nrow(moduleTraitCor))
head(moduleTraitCor)
#colnames(moduleTraitCor) <- colnames(moduleTraitCor) %>% substr(., 8, nchar(colnames(moduleTraitCor)))
corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
# corr_clust <- corr_clust[ordr(rownames(viewcolors)),]
correlation_ggplot <- corr_clust%>% melt()

dx <- dendro_data(sample_clust)
dy <- dendro_data(module_clust)
# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()

{heatmap4 <- ggplot(data = correlation_ggplot, 
                   aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#16316A", mid = "white", 
        high = "#61143F", na.value = "white", name = "Correlation", 
        limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  labs(x="Treatment", y="Clusters") +  theme_bw() +     theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "bottom")
 heatmap4}
{heatmap5 <- ggplot(data = correlation_ggplot, 
aes(x = Var2, y = Var1)) +
geom_tile(aes(fill =value), colour="darkgrey") +
scale_fill_gradient2(midpoint = 0, low = "#16316A", mid = "white", 
high = "#61143F", na.value = "white", name = "Correlation", 
limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
labs(x="Treatment") +  theme_bw() + theme(axis.text.y=element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "none", 
axis.title.y = element_blank())
heatmap5}

heatmap6 <- heatmap4 + heatmap5 + plot_layout(nrow = 1)
heatmap6
ggsave("pics/wgcnaHeatmapTreatOnly.pdf", heatmap6, dpi = 300, height = 7, width = 5)

#extracting gene names
library(xlsx)
library(openxlsx)

mediumpurple3 <- names(clustcolours[grep("mediumpurple3", clustcolours)])
greenyellow <- names(clustcolours[grep("greenyellow", clustcolours)])

# annotation file
FpAnnot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", header = T, row.names = 1)

MP3 <- FpAnnot[rownames(FpAnnot) %in% mediumpurple3,]
GY <- FpAnnot[rownames(FpAnnot) %in% greenyellow,]

write.xlsx2(MP3, file = "PSLib_wgcna_clusters.xlsx", sheetName = "mediumpurple3")
write.xlsx2(GY, file = "PSLib_wgcna_clusters.xlsx", sheetName = "greenyellow", append = T)

# GO enrichment
library(topGO)

geneID2GO <- readMappings(file = "GO/Drought.go", sep = "\t", IDsep = ";") 
allGenes <- rownames(fpkm.raw)

fisherStat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

geneFilter <- function(clust) {
    return(clust ==1 )
}

#############################3
# turq <- factor(as.integer(allGenes %in% turquoise))
# names(turq) <- allGenes
# turq_GO <- new("topGOdata", ontology = "BP", allGenes = turq,
#                geneSel = geneFilter, nodeSize = 10,
#                annotationFun = annFUN.gene2GO, gene2GO = geneID2GO)
# 
# turq_genes <- topGO::genes(turq_GO)
# sig <- sigGenes(turq_GO)
# turqTest <- getSigGroups(turq_GO, fisherStat)
# pvals <- score(turqTest)
# adjustPval <- p.adjust(pvals, method = "fdr", n = length(pvals))
# sigGenes <- adjustPval[which(adjustPval <= 0.05)]
# turqTable <- GenTable(turq_GO, classic = turqTest, topNodes = length(sigGenes))
# res <- merge(turqTable, sigGenes, by.x = "GO.ID", by.y = 0)                      
# 
# par(cex = 0.7, ps = 2)  
# showSigOfNodes(turq_GO, score(turqTest), firstSigNodes = 10, useInfo = "all")
# 
# 
# turquoise <- names(clustcolours)[clustcolours == "turquoise"]
# turq_genes <- factor(as.integer(allGenes %in% turquoise))
# names(turq_genes) <- allGenes
# turq_GOdata <- new("topGOdata", ontology = "BP", allGenes = turq_genes,
#                   annot = annFUN.gene2GO, gene2GO = geneID2GO)
# fisherStat <- new("classicCount", testStatistic = GOFisherTest, name =
#                      "Fisher test")
# turq_resFisher <- getSigGroups(turq_GOdata, fisherStat)
# turq_pvals <- score(turq_resFisher)
# turq_adjustPval <- p.adjust(turq_pvals, method = "fdr", n=length(turq_pvals))
# turq_sigGenes <- turq_adjustPval[which(turq_adjustPval <= 0.05)]
# turq_results <- GenTable(turq_GOdata, classic = turq_resFisher,
#                     topNodes = length(turq_sigGenes))
# turq_GO <- topGO::genes(turq_GOdata)
# res <- merge(turq_results, sigGenes, by.x = "GO.ID", by.y = 0)
# showSigOfNodes(turq_GOdata, topGO::score(turq_pvals), firstSigNodes = 10, useInfo = "all")
# colnames(res)[7] <- "p.adj"
# write.table(res[,c(1,7)], file = "GO/turquoise.tab", 
#             row.names = F, col.names = F,
#             quote = F, sep = "\t")

# streamline getting GO terms from each cluster 
goenrichment <- function(colour, clust = clustcolours, 
                         allG = allGenes, g2GO = geneID2GO) {
    ID <- names(clust)[clust == colour]
    genes <- factor(as.integer(allGenes %in% ID))
    names(genes) <- allG
print(genes)
    # TopGO objects
    GOdata <- new("topGOdata", ontology = "BP", allGenes = genes,
                  annot = annFUN.gene2GO, gene2GO = g2GO)
    fisherStat <- new("classicCount", testStatistic = GOFisherTest,
                      name = "Fisher test")

    # Fisher test
    resFisher <- getSigGroups(GOdata, fisherStat)
    pvals <- score(resFisher)
    adjustPval <- p.adjust(pvals, method = "fdr", n = length(pvals))
    sigGenes <- adjustPval[which(adjustPval <= 0.05)]

    res <- GenTable(GOdata, classic = resFisher,
                        topNodes = length(sigGenes))
    res <- merge(res, sigGenes, by.x = "GO.ID", by.y = 0)
    colnames(res)[7] <- "p.adj"
    
    # write GO and pval for ReviGO
    write.table(res[,c(1,7)], file = paste0("GO/", colour, ".tab", sep = ""),
                row.names = F, col.names = F,
                quote = F, sep = "\t")
    return(res)
}
############################

# extract GO terms
library(zoo)
## Revigo
revigo_summarize <- function(revigo, sig_go){
  revigo <- revigo %>% mutate(newdesc = ifelse(revigo$eliminated == 1, NA, revigo$description))
  revigo$newdesc <- na.locf(revigo$newdesc)
  rev_go <- cbind(revigo$term_ID, revigo$newdesc)
  colnames(rev_go) <- c("GO.ID", "Revigo")
  rev_merge <- merge(sig_go, rev_go,  by= "GO.ID") #%>% group_by(., newdesc) %>% 
#    summarize(DEG=sum(Significant), Expected=sum(Expected))
#  rev_merge <- rev_merge %>% melt(., id.vars="newdesc")
  return(rev_merge)
}

# clusters of interest
clust <- c("turquoise", "mediumpurple3", "greenyellow")

turq_res <- goenrichment("turquoise")
mp3 <- goenrichment("mediumpurple3")
gy <- goenrichment("greenyellow")

# merging revigo result with TopGO result
turq_rev <- read.csv("GO/turquoise_rev.csv")
turq_rev <- revigo_summarize(turq_rev, turq_res)
turq_rev <- turq_rev %>% arrange(desc(Annotated))

gy_rev <- read.csv("GO/greenyellow_rev.csv")
gy_rev <- revigo_summarize(gy_rev, gy)
gy_rev <- gy_rev %>% arrange(desc(Annotated))

# write to excel using xlsx2 method
# write.xlsx2(turq_rev, file = "GO/PS-GO.xlsx", sheetName = "turquoise",
#             append = F)
# write.xlsx2(gy_rev, file = "GO/PS-GO.xlsx", sheetName = "greenyellow", 
#             append = T)

# write to excel using openxlsx method
# Note: zip::zip() is deprecated, please use zip::zipr() instead
wb <- createWorkbook("GOenrichment")

addWorksheet(wb, "turquoise")
writeData(wb, sheet = "turquoise", turq_rev)
freezePane(wb, sheet = "turquoise",
           firstRow = T, firstCol = T)

addWorksheet(wb, "greenyellow")
writeData(wb, sheet = "greenyellow", gy_rev)
freezePane(wb, sheet = "greenyellow",
           firstRow = T, firstCol = T)

saveWorkbook(wb, "GO/PS-GO.xlsx", overwrite = T)

# k-means clustering step improves biological features of WGCNA gene co-expression
library(km2gcn)
