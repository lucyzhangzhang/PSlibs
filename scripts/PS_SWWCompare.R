library(DESeq2) #DEG analysis
library(dplyr) #data manipulation
setwd("~/R/Eutrema/PS") #stuff in here
condition <- as.factor(c(rep("SWW", 2), rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3)))

# lncRNA predictions using CREMA
# lmao okies
predictions <- read.csv("~/R/Eutrema/PS/crema/allG/final_ensemble_predictions.csv", header = T)
predictions$X <- toupper(predictions$X)
posiPred <- predictions[which(predictions$prediction == 1),]

# just a list of all the names that are used in this script
sampleNO <- c("SWW-1",
              "SWW-2",
              "ps12_S3",
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

# ps vs SWW
condition <- as.factor(c(rep("LPLS", 3), rep("SWW", 2)))

sampleNO <- c("ps12_S3",
              "ps48_S10",
              "ps49_S11",
              "SWW-1",
              "SWW-2")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(c(1, 1, 1, 2, 2)))

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
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "SWW")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
# [1] "Intercept"             "condition_LPLS_vs_SWW"

res.LPLSSWW <- results(dds, contrast = c("condition", "LPLS", "SWW"))
forMMLPLSSWW <- data.frame(names=rownames(res.LPLSSWW), PsPS=res.LPLSSWW$log2FoldChange)
#numbers
LL.up <- res.LPLSSWW[which(res.LPLSSWW[,"padj"] < 0.05 & res.LPLSSWW[,"log2FoldChange"] > 0),]
LL.dn <- res.LPLSSWW[which(res.LPLSSWW[,"padj"] < 0.05 & res.LPLSSWW[,"log2FoldChange"] < 0),]
#names
res.LPLSSWW.up <- res.LPLSSWW[which(res.LPLSSWW[,"padj"] < 0.05 & res.LPLSSWW[,"log2FoldChange"] > 0),]
res.LPLSSWW.dn <- res.LPLSSWW[which(res.LPLSSWW[,"padj"] < 0.05 & res.LPLSSWW[,"log2FoldChange"] < 0),]
# LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLSSWW.up),]
# LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLSSWW.dn),]

#group into matrix
LL.upp <- sum(rownames(LL.up) %in% posiPred$X)
LL.dnp <- sum(rownames(LL.dn) %in% posiPred$X)
LPLScompare <- matrix(c(nrow(LL.up), nrow(LL.dn), LL.upp, LL.dnp), ncol = 4, byrow = T)
LPLScompare

# Ps vs SWW
condition <- as.factor(c(rep("HPLS", 3), rep("SWW", 2)))

sampleNO <- c("ps37_S4",
              "ps41_S7",
              "ps50_S12",
              "SWW-1",
              "SWW-2")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(c(1, 1, 1, 2, 2)))

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
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "SWW")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
# [1] "Intercept"             "condition_HPLS_vs_SWW"

res.HPLSSWW <- results(dds, contrast = c("condition", "HPLS", "SWW"))
forMMHPLSSWW <- data.frame(names=rownames(res.HPLSSWW), PsPS=res.HPLSSWW$log2FoldChange)
#numbers
HL.up <- res.HPLSSWW[which(res.HPLSSWW[,"padj"] < 0.05 & res.HPLSSWW[,"log2FoldChange"] > 0),]
HL.dn <- res.HPLSSWW[which(res.HPLSSWW[,"padj"] < 0.05 & res.HPLSSWW[,"log2FoldChange"] < 0),]
#names
res.HPLSSWW.up <- res.HPLSSWW[which(res.HPLSSWW[,"padj"] < 0.05 & res.HPLSSWW[,"log2FoldChange"] > 0),]
res.HPLSSWW.dn <- res.HPLSSWW[which(res.HPLSSWW[,"padj"] < 0.05 & res.HPLSSWW[,"log2FoldChange"] < 0),]
# HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLSSWW.up),]
# HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLSSWW.dn),]

#group into matrix
HL.upp <- sum(rownames(HL.up) %in% posiPred$X)
HL.dnp <- sum(rownames(HL.dn) %in% posiPred$X)
HPLScompare <- matrix(c(nrow(HL.up), nrow(HL.dn), HL.upp, HL.dnp), ncol = 4, byrow = T)
HPLScompare

# pS vs SWW
condition <- as.factor(c(rep("LPHS", 3), rep("SWW", 2)))

sampleNO <- c("ps10_S1",
              "ps38_S5",
              "ps44_S8",
              "SWW-1",
              "SWW-2")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(c(1, 1, 1, 2, 2)))

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
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "SWW")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
# [1] "Intercept"             "condition_LPHS_vs_SWW"

res.LPHSSWW <- results(dds, contrast = c("condition", "LPHS", "SWW"))
forMMLPHSSWW <- data.frame(names=rownames(res.LPHSSWW), PsPS=res.LPHSSWW$log2FoldChange)
#numbers
LH.up <- res.LPHSSWW[which(res.LPHSSWW[,"padj"] < 0.05 & res.LPHSSWW[,"log2FoldChange"] > 0),]
LH.dn <- res.LPHSSWW[which(res.LPHSSWW[,"padj"] < 0.05 & res.LPHSSWW[,"log2FoldChange"] < 0),]
#names
res.LPHSSWW.up <- res.LPHSSWW[which(res.LPHSSWW[,"padj"] < 0.05 & res.LPHSSWW[,"log2FoldChange"] > 0),]
res.LPHSSWW.dn <- res.LPHSSWW[which(res.LPHSSWW[,"padj"] < 0.05 & res.LPHSSWW[,"log2FoldChange"] < 0),]
# LH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHSSWW.up),]
# LH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHSSWW.dn),]

#group into matrix
LH.upp <- sum(rownames(LH.up) %in% posiPred$X)
LH.dnp <- sum(rownames(LH.dn) %in% posiPred$X)
LPHScompare <- matrix(c(nrow(LH.up), nrow(LH.dn), LH.upp, LH.dnp), ncol = 4, byrow = T)
LPHScompare

# PS vs SWW
condition <- as.factor(c(rep("HPHS", 3), rep("SWW", 2)))

sampleNO <- c("ps37_S4",
              "ps41_S7",
              "ps50_S12",
              "SWW-1",
              "SWW-2")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(c(1, 1, 1, 2, 2)))

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
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "SWW")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
# [1] "Intercept"             "condition_HPHS_vs_SWW"

res.HPHSSWW <- results(dds, contrast = c("condition", "HPHS", "SWW"))
forMMHPHSSWW <- data.frame(names=rownames(res.HPHSSWW), PsPS=res.HPHSSWW$log2FoldChange)
#numbers
HH.up <- res.HPHSSWW[which(res.HPHSSWW[,"padj"] < 0.05 & res.HPHSSWW[,"log2FoldChange"] > 0),]
HH.dn <- res.HPHSSWW[which(res.HPHSSWW[,"padj"] < 0.05 & res.HPHSSWW[,"log2FoldChange"] < 0),]
#names
res.HPHSSWW.up <- res.HPHSSWW[which(res.HPHSSWW[,"padj"] < 0.05 & res.HPHSSWW[,"log2FoldChange"] > 0),]
res.HPHSSWW.dn <- res.HPHSSWW[which(res.HPHSSWW[,"padj"] < 0.05 & res.HPHSSWW[,"log2FoldChange"] < 0),]
# HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHSSWW.up),]
# HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHSSWW.dn),]

#group into matrix
HH.upp <- sum(rownames(HH.up) %in% posiPred$X)
HH.dnp <- sum(rownames(HH.dn) %in% posiPred$X)
HPHScompare <- matrix(c(nrow(HH.up), nrow(HH.dn), HH.upp, HH.dnp), ncol = 4, byrow = T)
HPHScompare

library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
DEGtable <- data.frame(rbind(LPLScompare, HPLScompare, LPHScompare, HPHScompare))
colnames(DEGtable) <- c("Upregulated", "Downregulated", "Pred_Up", "Pred_Down")
DEGtable$Upregulated <- DEGtable$Upregulated - DEGtable$Pred_Up
DEGtable$Downregulated <- DEGtable$Downregulated - DEGtable$Pred_Down
DEGtable <- data.frame(Treat=c("ps", "Ps", "pS", "PS"), DEGtable)
DEGtable$Downregulated <- -1 * DEGtable$Downregulated
DEGtable$"Pred_Down" <- -1 * DEGtable$"Pred_Down"

# "upp"= "#fcb3a7", "dnp"="#b5d1b4"
# "up"="#f54f31","dn"= "#1a8f10"
{DEG.melt <-  melt(DEGtable, id.vars = "Treat")
DEG.melt$label <- c(DEGtable$Upregulated + DEGtable$Pred_Up,
                    DEGtable$Downregulated + DEGtable$Pred_Down,
                    DEGtable$Pred_Up,
                    DEGtable$Pred_Down) 
 DEG.plot <-    ggplot(DEG.melt, aes(x = Treat, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + labs(x = "Treatment", y = "Number of DEGs") + theme_bw()+
    # scale_alpha_discrete(guide = F) +
geom_text_repel(data = DEG.melt[which(DEG.melt$variable == "Upregulated" | DEG.melt$variable == "Downregulated"),], 
                aes(x = Treat, y = label, label = abs(label), color = (variable)),
           alpha = 1, show.legend = F, seed = 7, point.padding = NA, direction = "y", nudge_y = DEG.melt$label[1:8]/35) +
scale_color_manual(values =c("#f54f31", "#1a8f10") ) +
scale_y_continuous(labels = abs) +
geom_text(data = DEG.melt[(DEG.melt$variable == "Pred_Up"), ], 
          aes(x = Treat, y = value, label = abs(value)), 
          alpha = 1, show.legend = F, color = "black", nudge_y = 80) +
geom_text(data = DEG.melt[(DEG.melt$variable == "Pred_Down"), ], 
          aes(x = Treat, y = value, label = abs(value)), 
          alpha = 1, show.legend = F, color = "black", nudge_y = -90) +
        scale_fill_manual(values = c("Pred_Up"="#f54f31","Pred_Down"= "#1a8f10","Upregulated"= "#fcb3a7", "Downregulated"="#b5d1b4"), 
                          name = "Expression", labels = str_wrap(c("Upregulated transcript", "Downregulated transcript", "Upregulated predicted lncRNA", "Downregulated predicted lncRNA" ), 13),
                          breaks = c("Upregulated", "Downregulated", "Pred_Up", "Pred_Down"), limits = c("Upregulated", "Downregulated", "Pred_Up", "Pred_Down"))
DEG.plot}
ggsave("pics/SWWcompare.pdf", DEG.plot, dpi = 350, height = 9, width = 4)

#############VENN############
library(venn)
# don't specify size when saving in pdf,  or else it won't print the text
# also don't specify the res value that is stup
pdf("pics/SWWupregVenn.pdf")
venn(list(ps=rownames(LL.up), Ps=rownames(HL.up), pS=rownames(LH.up), PS=rownames(HH.up)),
zcolor="style")
dev.off()
pdf("pics/SWWdownregVenn.pdf")
venn(list(ps=rownames(LL.dn), Ps=rownames(HL.dn), pS=rownames(LH.dn),  PS=rownames(HH.dn)),
zcolor="style")
dev.off()
venn(list(ps=union(rownames(LL.up), rownames(LL.dn)),
          Ps=union(rownames(HL.up), rownames(HL.dn)),
          pS=union(rownames(LH.up), rownames(LH.dn)),
          PS=union(rownames(HH.up), rownames(HH.dn))), 
zcolor="style")
#############END#############
