setwd("~/R/Eutrema/PS/")
library(WGCNA)
library(reshape2)
library(ggdendro)
library(ggplot2)
library(tidyverse)
library(viridis)
# library(patchwork) # do I really need this it's a ggplot api wrapper apparently
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

fpkm.raw <- read.table("FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)

gsg <- goodSamplesGenes(fpkm.raw, verbose = 3)
gsg$allOK

fpkm <- t(fpkm.raw)
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
network <- blockwiseModules(fpkm, power = 8,
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
# nGenes <- nrow(allgenes)
# nSamples <- ncol(allgenes);
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(fpkm, clustcolours)$eigengenes
MEs <- orderMEs(MEs0)
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
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,31);

#plotting
library(gplots)
library(RColorBrewer)

colnames(moduleTraitCor) <- c("ps", "Ps", "pS", "PS")

sample_clust <- hclust(dist(t(moduleTraitCor), method="euclidean"))
module_clust <- hclust(as.dist(1-cor(t(moduleTraitCor), method="pearson")))
moduleTraitCor[moduleTraitPvalue > 0.05] <- NA
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
#rownames(moduleTraitCor) <- rep("", nrow(moduleTraitCor))
head(moduleTraitCor)
#colnames(moduleTraitCor) <- colnames(moduleTraitCor) %>% substr(., 8, nchar(colnames(moduleTraitCor)))
corr_clust <- moduleTraitCor[module_clust$order, sample_clust$order]
correlation_ggplot <- corr_clust%>% melt()

pal <- colorRampPalette(c("#E7B800","#E7B800", "white", "#67A7C1", "#67A7C1"),
                        alpha=TRUE, interp="spline")(100)
droughtcol <- c("#C0C8C8","#8D9B9B",  "#6A90D1", "#5776AC",
                "#C0C8C8","#8D9B9B",  "#6A90D1", "#5776AC")
heatmap.2(moduleTraitCor,
          #cellnote = mat_data,  # same data set for cell labels
         # main = "Correlation", # heat map title
          #notecol="black",      # change font color of cell labels to black
          Rowv=as.dendrogram(module_clust), 
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # wid1Gens margins around plot
          #col=rev(brewer.pal(11,"PuOr")),       # use on color palette defined earlier
          col = pal,
          #breaks=col_breaks,    # enable color transition at specified limits
         Colv=as.dendrogram(sample_clust), 
         dendrogram="col",
          #ColSideColors= droughtcol,
         na.color="white",
         key.xlab="Correlation", key.title=NA, keysize=0.75,key.par = list(cex=0.5))

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

#correlation_ggplot <- correlation_ggplot %>% 
#  mutate(ecotype = ifelse(substr(Var2, 1,1) == "S", "Shandong", "Yukon")) %>%
#  mutate(Condition = substr(as.character(Var2), 2,nchar(as.character(Var2)))) %>%
#  mutate(Interaction = interaction(Condition, ecotype))

(heatmap <- ggplot(data = correlation_ggplot, 
                   aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill =value), colour="darkgrey") +
    scale_fill_gradient2(midpoint = 0, low = "#5ca4a9", mid = "white", 
        high = "#f9a03f", na.value = "white", name = "Correlation", 
        limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  labs(x="Ecotype and Condition", y="Clusters") +  theme_bw() + #theme(axis.text.y=element_blank(), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
    theme(axis.text.y=element_text(angle=20,vjust=0.5), plot.margin = margin(-0.75, 0, 0,0 , "cm")))
#plot.margin = margin(2, 2, 2, 2, "cm")
 (px + heatmap + plot_layout(ncol = 1, heights = c(1, 5))) # + plot_layout(ncol=2, widths=(1,5))

geneClusters <- network$colors
genenames <- rownames(fpkm.raw)
names(clustcolours) <- genenames
clustcolours[names(clustcolours) %>% toupper() %in% pattern3genes]
clustcolours[names(clustcolours) %>% toupper() %in% pattern4genes]

