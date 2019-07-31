library(DESeq2) #DEG analysis
library(plyr) #data manipulation
library(gplots) #graphing and visualization
library(genefilter) #GFF maker
library(dplyr) #graphing
library(geneplotter) #graphing
setwd("~/R/Eutrema/PS") #stuff in here

#####################SETUP###################
condition <- as.factor(c(rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3)))
sampleNO <- c("ps12_S3",
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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPHS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_HPHS"  "condition_HPLS_vs_HPHS" "condition_LPHS_vs_HPHS" "condition_LPLS_vs_HPHS"

res.HPLS <- results(dds, contrast = c("condition", "HPLS", "HPHS"))
res.LPHS <- results(dds, contrast = c("condition", "LPHS", "HPHS"))
res.LPLS <- results(dds, contrast = c("condition", "LPLS", "HPHS"))

#numbers
HL.up <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),])
HL.dn <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),])
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),])
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),])
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])
#names
res.HPLS.up <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),]
res.HPLS.dn <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),]
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

HL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.up),]
HL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.dn),]
LH.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]
LL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

#group into matrix
HPHScompare <- matrix(c(HL.up, HL.dn, LH.up, LH.dn, LL.up, LL.dn), ncol = 2, byrow = T)
HPHScompare
##########################HPLS########################
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPLS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)

res.HPHS <- results(dds, contrast = c("condition", "HPHS", "HPLS"))
res.LPHS <- results(dds, contrast = c("condition", "LPHS", "HPLS"))
res.LPLS <- results(dds, contrast = c("condition", "LPLS", "HPLS"))


HH.up <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),])
HH.dn <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),] )
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),] )
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),] )
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])

res.HPHS.up <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),]
res.HPHS.dn <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),]
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

HH.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.up),]
HH.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.dn),]
LH.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]
LL.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LL.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

HPLScompare <- matrix(c(HH.up, HH.dn, LH.up, LH.dn, LL.up, LL.dn), ncol = 2, byrow = T)
HPLScompare

############################LPHS###############################
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "LPHS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)

res.HPHS <- results(dds, contrast = c("condition", "HPHS", "LPHS"))
res.HPLS <- results(dds, contrast = c("condition", "HPLS", "LPHS"))
res.LPLS <- results(dds, contrast = c("condition", "LPLS", "LPHS"))


HH.up <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),])
HH.dn <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),])
HL.up <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),])
HL.dn <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),])
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])

res.HPHS.up <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),]
res.HPHS.dn <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),]
res.HPLS.up <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),]
res.HPLS.dn <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),]
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

HH.LH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.up),]
HH.LH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.dn),]
HL.LH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.up),]
HL.LH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.dn),]
LL.LH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LL.LH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

LPHScompare <- matrix(c(HH.up, HH.dn, HL.up, HL.dn, LL.up, LL.dn), ncol = 2, byrow = T)
LPHScompare
##########################LPLS########################
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "LPLS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)

res.HPHS <- results(dds, contrast = c("condition", "HPHS", "LPLS"))
res.HPLS <- results(dds, contrast = c("condition", "HPLS", "LPLS"))
res.LPHS <- results(dds, contrast = c("condition", "LPHS", "LPLS"))


HH.up <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),])
HH.dn <- nrow(res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),] )
HL.up <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),])
HL.dn <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),])
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),] )
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),] )

res.HPHS.up <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] > 0),]
res.HPHS.dn <- res.HPHS[which(res.HPHS[,"padj"] < 0.05 & res.HPHS[,"log2FoldChange"] < 0),]
res.HPLS.up <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),]
res.HPLS.dn <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),]
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]

HH.LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.up),]
HH.LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPHS.dn),]
HL.LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.up),]
HL.LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.dn),]
LL.LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LL.LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]

LPLScompare <- matrix(c(HH.up, HH.dn, HL.up, HL.dn, LH.up, LH.dn), ncol = 2, byrow = T)
LPLScompare

data <- rbind(HPHScompare, HPLScompare, LPHScompare, LPLScompare)
data <- cbind(c(rep("HPHS", 3), rep("HPLS", 3), rep("LPHS", 3), rep("LPLS", 3)), c("HPLS", "LPHS", "LPLS", "HPHS", "LPHS", "LPLS", "HPHS", "HPLS", "LPLS", "HPHS", "HPLS", "LPHS"), data)
colnames(data) <- c("denom", "numerator", "up", "down")
data <- as.data.frame(data, stringsAsfactors = F)
rownames(data) <- NULL
data

#latex tabular form
write.table(data, "~/committee/DEGcount", sep = " & ", quote = F, eol = " \\\\\n", row.names = F)

##################Individual effects of P/S i.e. combine the other treatment##########
# P
condition <- as.factor(c(rep(c("LP", "HP", "LP", "HP"), each = 3)))
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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HP")

dds <- DESeq(DESeq2Table)

resultsNames(dds)
res.HP <- results(dds, contrast = c("condition", "LP", "HP"))

HP.up <- nrow(res.HP[which(res.HP[,"padj"] < 0.05 & res.HP[,"log2FoldChange"] > 0),] )
HP.dn <- nrow(res.HP[which(res.HP[,"padj"] < 0.05 & res.HP[,"log2FoldChange"] < 0),] )
HPcompare <- c(HP.up, HP.dn)
HPcompare

###
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "LP")

dds <- DESeq(DESeq2Table)

resultsNames(dds)
res.LP <- results(dds, contrast = c("condition", "HP", "LP"))

LP.up <- nrow(res.LP[which(res.LP[,"padj"] < 0.05 & res.LP[,"log2FoldChange"] > 0),] )
LP.dn <- nrow(res.LP[which(res.LP[,"padj"] < 0.05 & res.LP[,"log2FoldChange"] < 0),] )
LPcompare <- c(LP.up, LP.dn)
LPcompare


res.LP[which(res.LP[,"padj"] < 0.05),]
res.HP[which(res.HP[,"padj"] < 0.05),]

# S
condition <- as.factor(c(rep(c("LS", "LS", "HS", "HS"), each = 3)))
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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)
res.HS <- results(dds, contrast = c("condition", "LS", "HS"))

HS.up <- nrow(res.HS[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] > 0),] )
HS.dn <- nrow(res.HS[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] < 0),] )
HScompare <- c(HS.up, HS.dn)
HScompare

###
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "LS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)
res.LS <- results(dds, contrast = c("condition", "HS", "LS"))

LS.up <- nrow(res.LS[which(res.LS[,"padj"] < 0.05 & res.LS[,"log2FoldChange"] > 0),] )
LS.dn <- nrow(res.LS[which(res.LS[,"padj"] < 0.05 & res.LS[,"log2FoldChange"] < 0),] )
LScompare <- c(LS.up, LS.dn)
LScompare


res.HS.up <- res.HP[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] > 0),]
res.HS.dn <- res.HP[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] < 0),]

#load annotation and FPKM files, generated previously with DESeq2, annotation file generated by Caitlin
expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)

pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.up),]

neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.dn),]

library(xlsx)
library(openxlsx)

write.xlsx2(pos, file = "Sulfurcomparison.xlsx", sheetName = "plusS", append = F, row.names = T)
write.xlsx2(neg, file = "Sulfurcomparison.xlsx", sheetName = "minusS", append = T, row.names = T)


