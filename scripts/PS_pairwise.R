library(DESeq2) #DEG analysis
library(dplyr) #data manipulation
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

expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)
#specifying the reference level
keep <- rowSums(counts(DESeq2Table)) >= 1
DESeq2Table <- DESeq2Table[keep,]
#correct for batch effect
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2, 3))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPHS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)

{res.HPLS <- results(dds, contrast = c("condition", "HPLS", "HPHS"))
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

# HL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.up),]
# HL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.dn),]
# LH.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
# LH.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]
# LL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
# LL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]
# 
# write.table(rownames(HL.HH.pos), file = "names/HH/HL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(HL.HH.neg), file = "names/HH/HL/dn/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LH.HH.pos), file = "names/HH/LH/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LH.HH.neg), file = "names/HH/LH/dn/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LL.HH.pos), file = "names/HH/LL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LL.HH.neg), file = "names/HH/LL/dn/names", row.names = F, quote = F, col.names = F)
#group into matrix
HPHScompare <- matrix(c(HL.up, HL.dn, LH.up, LH.dn, LL.up, LL.dn), ncol = 2, byrow = T)
HPHScompare}
##########################HPLS########################
DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPLS")

dds <- DESeq(DESeq2Table)

resultsNames(dds)
{
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
}
write.table(rownames(LH.HL.pos), file = "names/HL/LH/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.HL.neg), file = "names/HL/LH/dn/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HL.pos), file = "names/HL/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HL.neg), file = "names/HL/LL/dn/names", row.names = F, quote = F, col.names = F)

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

# write.table(rownames(LL.LH.pos), file = "names/LH/LL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LL.LH.neg), file = "names/LH/LL/dn/names", row.names = F, quote = F, col.names = F)
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
LH.LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]

LPLScompare <- matrix(c(HH.up, HH.dn, HL.up, HL.dn, LH.up, LH.dn), ncol = 2, byrow = T)
LPLScompare

data <- rbind(HPHScompare, HPLScompare, LPHScompare, LPLScompare)
data <- cbind(c(rep("HPHS", 3), rep("HPLS", 3), rep("LPHS", 3), rep("LPLS", 3)), c("HPLS", "LPHS", "LPLS", "HPHS", "LPHS", "LPLS", "HPHS", "HPLS", "LPLS", "HPHS", "HPLS", "LPHS"), data)
colnames(data) <- c("denom", "numerator", "up", "down")
data <- as.data.frame(data, stringsAsfactors = F)
rownames(data) <- NULL
data

#latex tabular form
# write.table(data, "~/committee/DEGcount", sep = " & ", quote = F, eol = " \\\\\n", row.names = F)

##################Individual effects of P/S i.e. combine the other treatment##########
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

res.HS.up <- res.HS[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] > 0),]
res.HS.dn <- res.HS[which(res.HS[,"padj"] < 0.05 & res.HS[,"log2FoldChange"] < 0),]

# write.table(rownames(res.HS.up), file = "names/HS/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(res.HS.dn), file = "names/HS/dn/names", row.names = F, quote = F, col.names = F)
#load annotation and FPKM files, generated previously with DESeq2, annotation file generated by Caitlin
expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)

#pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.up),]

#neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.dn),]

library(xlsx)
library(openxlsx)

# write.xlsx2(pos, file = "Sulfurcomparison.xlsx", sheetName = "plusS", append = F, row.names = T)
# write.xlsx2(neg, file = "Sulfurcomparison.xlsx", sheetName = "minusS", append = T, row.names = T)

################Separate LPLS because it fucks up distribution calc#######
condition <- as.factor(c(rep(c("HPLS", "HPHS"), each = 3)))
sampleNO <- c("ps11_S2",
              "ps40_S6",
              "ps46_S9",
              "ps37_S4",
              "ps41_S7",
              "ps50_S12")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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
forMMHPLS <- data.frame(names=rownames(res.HPLS), PsPS=res.HPLS$log2FoldChange)
#numbers
HL.up <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),])
HL.dn <- nrow(res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),])
#names
res.HPLS.up <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] > 0),]
res.HPLS.dn <- res.HPLS[which(res.HPLS[,"padj"] < 0.05 & res.HPLS[,"log2FoldChange"] < 0),]


HL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.up),]
HL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HPLS.dn),]

#group into matrix
HPHScompare <- matrix(c(HL.up, HL.dn), ncol = 2, byrow = T)
HPHScompare

# write.table(rownames(HL.HH.pos), file = "names/HH/HL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(HL.HH.neg), file = "names/HH/HL/dn/names", row.names = F, quote = F, col.names = F)
######################################################################
condition <- as.factor(c(rep(c("LPHS", "HPHS"), each = 3)))
sampleNO <- c(              "ps10_S1",
              "ps38_S5",
              "ps44_S8",
              "ps37_S4",
              "ps41_S7",
              "ps50_S12")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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

res.LPHS <- results(dds, contrast = c("condition", "LPHS", "HPHS"))
forMMLPHSS <- data.frame(names=rownames(res.LPHS), pSPS=res.LPHS$log2FoldChange)

#numbers
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),])
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),])
#names
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]

LH.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]
# write.table(rownames(LH.HH.pos), file = "names/HH/LH/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LH.HH.neg), file = "names/HH/LH/dn/names", row.names = F, quote = F, col.names = F)

#group into matrix
HPHScompare <- matrix(c(LH.up, LH.dn), ncol = 2, byrow = T)
HPHScompare

##############################LPLS vs HPHS##################################
condition <- as.factor(c(rep(c("LPLS", "HPHS"), each = 3)))
sampleNO <- c("ps12_S3",
              "ps48_S10",
              "ps49_S11",
              "ps37_S4",
              "ps41_S7",
              "ps50_S12")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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

res.LPLS2 <- results(dds, contrast = c("condition", "LPLS", "HPHS"))
forMMLPLS2S <- data.frame(names=rownames(res.LPLS2), psPS=res.LPLS2$log2FoldChange)

#numbers
LL.up <- nrow(res.LPLS2[which(res.LPLS2[,"padj"] < 0.05 & res.LPLS2[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS2[which(res.LPLS2[,"padj"] < 0.05 & res.LPLS2[,"log2FoldChange"] < 0),])
#names
res.LPLS2.up <- res.LPLS2[which(res.LPLS2[,"padj"] < 0.05 & res.LPLS2[,"log2FoldChange"] > 0),]
res.LPLS2.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

LL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS2.up),]
LL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS2.dn),]

# write.table(rownames(LL.HH.pos), file = "names/HH/LL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LL.HH.neg), file = "names/HH/LL/dn/names", row.names = F, quote = F, col.names = F)
#group into matrix
HPHScompare <- matrix(c(LL.up, LL.dn), ncol = 2, byrow = T)
HPHScompare

condition <- as.factor(c(rep(c("HPLS", "LPHS"), each = 3)))
sampleNO <- c(              "ps11_S2",
              "ps40_S6",
              "ps46_S9",
              "ps10_S1",
              "ps38_S5",
              "ps44_S8"
             )

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPLS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_HPHS"  "condition_HPLS_vs_HPHS" "condition_LPHS_vs_HPHS" "condition_LPLS_vs_HPHS"

res.LPHS <- results(dds, contrast = c("condition", "LPHS", "HPLS"))
forMMLPHS <- data.frame(names=rownames(res.LPHS), PspS=res.LPHS$log2FoldChange)

#numbers
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),])
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),])
#names
res.LPHS2.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS2.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]

LH.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS2.up),]
LH.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS2.dn),]

# write.table(rownames(LH.HL.pos), file = "names/HL/LH/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LH.HL.neg), file = "names/HL/LH/dn/names", row.names = F, quote = F, col.names = F)
#group into matrix
HPHScompare <- matrix(c(LH.up, LH.dn), ncol = 2, byrow = T)
HPHScompare
condition <- as.factor(c(rep(c("LPLS", "HPLS"), each = 3))) 
sampleNO <- c("ps12_S3",
              "ps48_S10",
              "ps49_S11",
              "ps11_S2",
              "ps40_S6",
              "ps46_S9")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPLS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_HPHS"  "condition_HPLS_vs_HPHS" "condition_LPHS_vs_HPHS" "condition_LPLS_vs_HPHS"

res.LPLS <- results(dds, contrast = c("condition", "LPLS", "HPLS"))
forMMLPLSS <- data.frame(names=rownames(res.LPLS), psPs=res.LPLS$log2FoldChange) 

#numbers
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])
#names
res.LPLS3.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS3.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

LL.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS3.up),]
LL.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS3.dn),]

#group into matrix
HPHScompare <- matrix(c(LL.up, LL.dn), ncol = 2, byrow = T)
# HPHScompare
# write.table(rownames(LL.HL.pos), file = "names/HL/LL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LL.HL.neg), file = "names/HL/LL/dn/names", row.names = F, quote = F, col.names = F)

condition <- as.factor(c(rep(c("LPLS", "LPHS"), each = 3)))
sampleNO <- c("ps12_S3",
              "ps48_S10",
              "ps49_S11",
                            "ps10_S1",
              "ps38_S5",
              "ps44_S8")

metadata <- data.frame(sampleNO = sampleNO,
                       condition = condition,
                       libraryName = sampleNO,
                        reps = as.factor(rep(c(1, 2, 3), 2)))

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

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "LPHS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)
#[1] "Intercept"              "condition_F_vs_HPHS"  "condition_HPLS_vs_HPHS" "condition_LPHS_vs_HPHS" "condition_LPLS_vs_HPHS"

res.LPLS <- results(dds, contrast = c("condition", "LPLS", "LPHS"))
forMMLPLS <- data.frame(names = rownames(res.LPLS), pspS = res.LPLS$log2FoldChange)

#numbers
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])
#names
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

LH.LL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LH.LL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

#group into matrix
LPLScompare <- matrix(c(LL.up, LL.dn), ncol = 2, byrow = T)
LPLScompare

# write.table(rownames(LH.LL.pos), file = "names/LH/LL/up/names", row.names = F, quote = F, col.names = F)
# write.table(rownames(LH.LL.neg), file = "names/LH/LL/dn/names", row.names = F, quote = F, col.names = F)

# Venn diagrams
library(VennDiagram)
#S
#runs out of order, just run the required groups above
Sup <- venn.diagram(x=list(PSxPs=rownames(HL.HH.pos),pSxps=rownames(LH.LL.pos)), filename = NULL)
Sdn <- venn.diagram(x=list(PSxPs=rownames(HL.HH.neg),pSxps=rownames(LH.LL.neg)), filename = NULL)
#draw plot
dev.off()
grid.draw(Sup)
plot.new()
grid.draw(Sdn)
dev.off()

Sinter <- intersect(rownames(HL.HH.pos),rownames(LH.LL.pos))
Sinter.descs <- expr_annot[rownames(expr_annot) %in% inter,]

#P
Pup <- venn.diagram(x=list(PSxpS=rownames(LH.HH.pos),Psxps=rownames(LL.HL.pos)), filename = NULL)
Pdn <- venn.diagram(x=list(PSxpS=rownames(LH.HH.neg),Psxps=rownames(LL.HL.neg)), filename = NULL)


plot.new()
grid.draw(Pup)
plot.new()
grid.draw(Pdn)

#############################################################################################

# I don't fuckn know man, wtf is up with this shit, so damn inconsistent

fpkm.raw <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", row.names=1, header = T)
fpkm <- fpkm.raw[,c(1:12)]

ps <- log2(fpkm[,1:3] + 1)
Ps <- log2(fpkm[,4:6] + 1)
pS <- log2(fpkm[,7:9] + 1)
PS <- log2(fpkm[,10:12] + 1)

pS <- rowMeans(pS)
Ps <- rowMeans(Ps)
ps <- rowMeans(ps)
PS <- rowMeans(PS)

ps.mean <- cbind(ps, pS, Ps, PS)
colnames(ps.mean) <- c("ps", "pS", "Ps", "PS")


## all genes
pca<- prcomp(ps.mean, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)

dim(exprVals)
dim(sampleVals)


pos.PC <- exprVals[rownames(exprVals) %in% rownames(pos),] 
neg.PC <- exprVals[rownames(exprVals) %in% rownames(neg),]


pos.combined <- as.data.frame(append(pos, 
                                     list(PC2=pos.PC$PC2, PC4=pos.PC$PC4), 
                                     after = 17))
neg.combined <- as.data.frame(append(neg, 
                                     list(PC2=neg.PC$PC2, PC4=neg.PC$PC4), 
                                     after = 17))

rownames(pos.combined) <- rownames(pos)
rownames(neg.combined) <- rownames(neg)

pos.combined <- pos.combined[order(pos.combined$PC2, decreasing = T),]
neg.combined <- neg.combined[order(neg.combined$PC2, decreasing = T),]

#write to excel file lmao

# write.xlsx2(pos.combined, file = "ScomparisonPCs.xlsx", 
#             sheetName = "plusS", append = F, row.names = T)
# write.xlsx2(neg.combined, file = "ScomparisonPCs.xlsx", 
#             sheetName = "minusS", append = T, row.names = T)


######################################################################################
# GO term semantic redundancy

library(topGO)

# load DEGs
 wd <- "~/R/Eutrema/PS/GO/"
setwd(wd)

allGenes <- rownames(fpkm.raw)
geneID2GO <- readMappings(file = "Drought.go", sep = "\t", IDsep = ";") 
Sup <- scan("~/R/Eutrema/PS/crema/DEGs/HS/up/names", what = "character") 
Sdn <- scan("~/R/Eutrema/PS/crema/DEGs/HS/dn/names", what = "character") 

S <- c(Sup, Sdn)

Sgenes <- factor(as.integer(allGenes %in% S)) 
names(Sgenes) <- allGenes

pipeLine <- function(genes, name = "topGO") {
	fisherStat <- new("classicCount", testStatistic = GOFisherTest, 
                          name = "Fisher test")
	genefilter <- function(score){
	    return(score == 1)
	}
	GOdat <- new("topGOdata",
	             description = "DEGs",
	             ontology = "BP",
	             allGenes = genes,
	             geneSel = genefilter,
	             annot = annFUN.gene2GO, gene2GO = geneID2GO
	)
	
	GOgenes <- topGO::genes(GOdat)
	
	resFtest <- getSigGroups(GOdat, fisherStat)
	pvals <- score(resFtest)
	pAdjust <- p.adjust(pvals, method = "fdr", n = length(pvals))
	sigGenes <- pAdjust[which(pAdjust <= 0.05 )]
	
	res <- GenTable(GOdat, classic = resFtest, 
                        topNodes = length(sigGenes))
	
	res <- merge(res, sigGenes, by.x = "GO.ID", by.y = 0)
	    colnames(res)[7] <- "p.adj"
	    write.table(res[,c(1,7)], file = paste0(name, ".tab", sep = ""),
	                row.names = F, col.names = F,
	                quote = F, sep = "\t")
	
	return(res)
}

pipeLine(Sgenes, name = "S")



# Top GO all DEGs

DEG.names <- scan("~/R/Eutrema/PS/crema/DEGs/allNames", what = "character")
DEGs <- factor(as.integer(allGenes %in% DEG.names))
names(DEGs) <- allGenes

DEG.result <- pipeLine(DEGs, name = "DEG_GO")

DEG.rev <- read.csv("DEG_rev.csv", header = T)

DEG.rev <- as.data.frame(DEG.rev, StringsAsFactors = F)

library(ggplot2)
library(ggrepel)

DEG.rev <- mutate(DEG.rev, semSim = 1-uniqueness)
DEG.rev <- mutate(DEG.rev, Freq = as.numeric(gsub("%$", "", frequency))) 

{GO.plot <- ggplot(DEG.rev, aes(x = log10.p.value, y = semSim)) +
    labs(x = "log10 p-value", y = "GO Term semantic similarity") +
    geom_vline(xintercept = log10(0.05), linetype = "dashed",
               color = "grey") + 
    geom_text(data = DEG.rev, mapping = aes(x = log10(0.05), y = 0, 
                                            label = "log10(0.05)"), color = "red",
              size = 2) + 
    geom_point(color = "black", size = 1.8) +
    geom_point(aes(color = semSim), size = 1.5) +
    xlim(-30, 6) +
    scale_color_gradientn(colors = rev(rainbow(8))) +
    geom_label_repel(data = subset(DEG.rev,grepl("(.*stress(.*|$)| .*glu.*|.*nutr.*|.*sulf.*|.*starv.*|.*auxin.*|(.*|$)S-gly.*|.*phosp.*)", DEG.rev$description)),
                     aes(x = log10.p.value, y = semSim,
                         label = description), min.segment.length = unit(0, "lines"), 
                     size = 2, segment.size = 0.2, alpha = 0.5, seed = 12) + 
    geom_label_repel(data = subset(DEG.rev,grepl("(.*stress(.*|$)| .*glu.*|.*nutr.*|.*sulf.*|.*starv.*|.*auxin.*|(.*|$)S-gly.*|.*phosp.*)", DEG.rev$description)),
                     aes(x = log10.p.value, y = semSim,
                         label = description), min.segment.length = unit(0, "lines"), 
                     size = 2, segment.size = 0.2, fill = NA, seed = 12) + 
theme_bw()+
theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=15))
GO.plot}
ggsave("../pics/goSemSim.pdf", GO.plot, dpi = 350, height = 4, width = 5.5)
pdf("../pics/goSemSim.pdf", width = 6, height = 6)
GO.plot
dev.off()
# all labels
{GO.plot.all <- ggplot(DEG.rev, aes(x = log10.p.value, y = semSim)) +
    labs(x = "log10 p-value", y = "GO Term semantic similarity") +
    geom_vline(xintercept = log10(0.05), linetype = "dashed",
               color = "grey") + 
    geom_text(data = DEG.rev, mapping = aes(x = log10(0.05), y = 0, 
                label = "log10(0.05)"), color = "red", size = 2) + 
    xlim(-30, 5) + ylim(0, 0.5) +
    geom_point(color = "black", size = 1.8) +
    geom_point(aes(color = semSim), size = 1.5) +
    scale_color_gradientn(colors = rev(rainbow(8))) +
    geom_label_repel(data = DEG.rev,
                     aes(x = log10.p.value, y = semSim,
                         label = description), min.segment.length = unit(0, "lines"), 
                     size = 2, segment.size = 0.2, alpha = 0.5, seed = 12) + 
    geom_label_repel(data = DEG.rev,
                     aes(x = log10.p.value, y = semSim,
                         label = description), min.segment.length = unit(0, "lines"), 
                     size = 2, segment.size = 0.2, fill = NA, seed = 12) + 
    theme_bw()+ theme(legend.position = "none")
GO.plot.all}
ggsave("../pics/goSemSimAll.pdf", GO.plot.all, dpi = 350, height = 9, width = 14)

annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", header = T, row.names = 1)

auxin1 <-topGO::genesInTerm(GOdat, "GO:0009850")
auxin2 <-topGO::genesInTerm(GOdat, "GO:0009851")
auxin3 <- union(unlist(auxin1), unlist(auxin2))
auxin3 <- auxin3[which(auxin3 %in% DEG.names)]
auxins <- annot[which(rownames(annot) %in% auxin3),c(1:18, 29)]
colnames(auxins)[6:17] <- paste0(c("ps", "Ps", "pS", "PS"), 1:3)

write.xlsx2(auxins, file = "differentialAuxins.xlsx", sheetName = "auxins", append = F)
sulfurs <- subset(DEG.rev, grepl(".*sulf.*", DEG.rev$description))
sulf1 <- genesInTerm(GOdat, "GO:0006790")
sulf2 <- genesInTerm(GOdat, "GO:0044272")
S <- union(unlist(sulf1), unlist(sulf2))
length(S[which(S %in% DEG.names)])

##########################STRING############################
library(tidyr)
library(igraph)
library(httr)
res.HS.up
res.HS.dn

# Arabidopsis conversion list
araConvert <- read.table("~/Eutrema/FPKM/eutremaToArabidopsis.names")
names(araConvert) <- c("Eutr", "Ara")
# length(unique(araConvert$Eutr))
# Sup <- data.frame(Eutr=rownames(res.HS.up))
# Sup <- merge(Sup, araConvert, by = "Eutr")
# 
# Sdn <- data.frame(Eutr=rownames(res.HS.dn))
# Sdn <- merge(Sdn, araConvert, by = "Eutr")
#     better than string I guess??
# Sgenes <- rbind(Sup, Sdn)
#testing
# ConvertList <- genes
# Sgenes$Expr <- ifelse(Sgenes$Eutr %in% Sup$Eutr, "Sup", "Sdn")
making_a_network <- function(ConvertList) {
	interactions <- data.frame()
        araGenes <- list()
        araComplete <- ConvertList[complete.cases(ConvertList),]
        araComplete <- araComplete[grepl("AT.*", araComplete$Ara),]
        isolates <- ConvertList[!grepl("AT.*", ConvertList$Ara),]
        pb <- txtProgressBar(min = 0, max = nrow(araComplete), initial = 0, style = 3)
	for (i in 1:nrow(araComplete)) {
            # freakking smacc
	    httpReq <- GET(paste0("http://string-db.org/api/tsv/network?identifier=",
	                             araComplete[i,2], "&species=3702&required_score=700")) 
	    httpReq <- httr::content(httpReq, as = "text", encoding = "ISO-8859-1")
            entry <- araComplete[i,]
            araGeneNames <- data.frame()
	    # read tabel from string value if the data set is smol enough
	    # can't specify header when reading from a string
            if (!grepl("Error.*", httpReq)) {
                # replace the first row, which is a header
                httpReq <- read.table(text = httpReq, stringsAsFactors = F)
                if (nrow(httpReq) >= 2) {
            	names(httpReq) <- httpReq[1,] 
                httpReq <- httpReq[2:nrow(httpReq),]

                # what the common name is corresponding to the Gene ID
                #                 localInteraction <- data.frame()
                #                 for (j in 1:nrow(httpReq)) {
                # ID_A ID_B COMMON_NAME_A COMMON_NAME_B
                #filter the scores
                # pscore = 9, ascore = 10, escore = 11, dscore = 12
                    #                 Arow <- httpReq[j,]
                    #                 if (any(Arow[9:12] > 0)) {
                    #                     localInteraction <- rbind(localInteraction, 
                    #                                               data.frame(as.character(Arow[3]), 
                    #                                                          as.character(Arow[4])))
                araGeneNames <- cbind(c(httpReq[,1], httpReq[,2]), c(httpReq[,3], httpReq[,4]))
                interactions <- rbind(interactions, data.frame(as.character(httpReq[,3]), 
                                                               as.character(httpReq[,4]))) 
                #                 }
                #         }
                #             if ((nrow(localInteraction) > 0)) {
                #                 interactions <- rbind(interactions, localInteraction)
                araGenes <- rbind(araGenes, araGeneNames)
                #                     }
            }  else {
                    isolates <- rbind(isolates, entry)
                }
            } else {
                    isolates <- rbind(isolates, entry)
            }

            setTxtProgressBar(pb, i)
                     
	}	    

        return(list(Int=interactions, Name=araGenes, Iso=isolates))
}

library(tidyverse)
library(qgraph)
library(semnet)

# HPLS vs. LPHS 127D 11 U

# pS vs. ps 73U 18D

# PS vs. ps 17U 3D

# PS vs. pS 2U 0D
# Ps vs. ps 0U 11D

# pS vs. PS 2U 0D
# res.LPHS2.up
# res.LPHS2.dn
# LPHS2up <- data.frame(Eutr=rownames(res.LPHS2.up))
# LPHS2up <- merge(LPHS2up, araConvert, by = "Eutr")
# 
# LPHS2dn <- data.frame(Eutr=rownames(res.LPHS2.dn))
# LPHS2dn <- merge(LPHS2dn, araConvert, by = "Eutr")
#     better than string I guess??
# LPHS2genes <- rbind(LPHS2up, LPHS2dn)
# LPHS2genes$Expr <- ifelse(LPHS2genes$Eutr %in% LPHS2up$Eutr, "Sup", "Sdn")


# DEGlist <- list(
# A_S_vs_s=Sgenes,
# B_Ps_vs_PS=HPHSgenes,
# C_pS_vs_Ps=LPHSgenes,
# D_ps_vs_pS=LPLSgenes,
# E_ps_vs_PS=LPLS2genes,
# F_ps_vs_Ps=LPLS3genes,
# G_pS_vs_PS=LPHS2genes)
# 
# library(openxlsx)
# library(xlsx)
# 
# for (i in 1:length(DEGlist)) {
#     if (i == 1) {
#         App = F
#     } else {
#         App = T
#     }
#     write.xlsx2(DEGlist[[i]], file = "DEGsList.xlsx", sheetName = names(DEGlist[i]), append = App)
# }


# S and HPHS
Snames <- Sgenes$Eutr
HPHSnames <- HPHSgenes$Eutr
SHnames <- data.frame(intersect(Snames,HPHSnames))
SHnames$Expr <- ifelse(SHnames[,1] %in% union(Sgenes[which(Sgenes$Expr == "Sup"),]$Eutr,
                                              HPHSgenes[which(HPHSgenes$Expr == "Sup"), ]$Eutr), "Up", "Down")
names(SHnames) <- c("Eutr", "Expr")

FPKM$Eutr <- rownames(FPKM)
SHmerge <- merge(FPKM, SHnames, by = "Eutr")
write.xlsx2(SHmerge, file = "S_DEG_Overlap.xlsx", append = F)

draw_network <- function(uplist, dnlist,  name) {
    # Reorganize the data
    up <- data.frame(Eutr=rownames(uplist))
    up <- merge(up, araConvert, by = "Eutr")
    
    dn <- data.frame(Eutr=rownames(dnlist))
    dn <- merge(dn, araConvert, by = "Eutr")
    # better than string I guess??
    genes <- rbind(up, dn)
    genes$Expr <- ifelse(genes$Eutr %in% up$Eutr, "Sup", "Sdn")

    # STRING database query
    StringDat <- making_a_network(genes)

    interactions <- StringDat[[1]]
    Conversions <- StringDat[[2]]
    isolates <- StringDat[[3]]
    Conversions <- unique(Conversions)
    isolates <- unique(isolates)
    colnames(Conversions) <- c("Ara", "Name")
    
    ConMerge <- merge(genes[complete.cases(genes),], Conversions, by = "Ara")
    ConMerge$Name <- unlist(ConMerge$Name)
    ConMerge <- ConMerge %>% distinct()
    EutrGenes <- merge(araConvert, Conversions, by = "Ara")

    # calculating the graph object
    qnode <- interactions %>% distinct()
    qnode <- as.matrix(qnode)
    
    # Finding unique nodes
    UniqueNodes <- unique(c(as.character(interactions[,1]), as.character(interactions[,2])))    
    names(UniqueNodes) <- 1:length(UniqueNodes)
    SuniqueNodes <- as.character(ConMerge$Name)
    
    # Defining groups for colors
    Up <- UniqueNodes[which(UniqueNodes %in% ConMerge[which(ConMerge$Expr == "Sup"),]$Name)]
    Down <- UniqueNodes[which(UniqueNodes %in% ConMerge[which(ConMerge$Expr == "Sdn"),]$Name)]
    Eutr <- UniqueNodes[which(UniqueNodes %in% setdiff(EutrGenes$Name,  union(Up, Down)))]
    Ara <- UniqueNodes[-which(UniqueNodes %in% union(Eutr, SuniqueNodes))]
            
    # Groups
    qgroups <- list(Up=as.numeric(names(Up)),
                    Down=as.numeric(names(Down)),
                    Eutr=as.numeric(names(Eutr)),
                    Ara=as.numeric(names(Ara)))

    # Draw graph
    #     qgraph(qnode, directed = F, groups = qgroups, 
    #            color = c("darksalmon", "lightskyblue", "azure3", "white"),
    #     edge.width = 0.3, filetype = "pdf", filename = paste0("pics/S"),
    #     border.width = 0.4, repulsion = 0.92, vsize=0.8, normalize = F)
    # 
    qgraph(qnode, directed = F, groups = qgroups, 
           color = c("darksalmon", "lightskyblue", "azure3", "white"),
    edge.width = 0.3, filetype = "pdf", filename = paste0("pics/", name),
    border.width = 0.4, repulsion = 0.92, vsize=0.8, normalize = F)

    # Write isolates in a latex table
    write.table(isolates, paste0(name, "DEG_isolates.tab"), sep = " & ", 
                row.names = F, col.names = T, quote = F, eol = " \\\\\n")
    return(list(DEGs=ConMerge, other=Eutr, Convert=Conversions, Interactions = interactions, groups = qgroups))
}


Snetwork <- draw_network(res.HS.up, res.HS.dn, "S")
HPLSnetwork <- draw_network(res.HPLS.up, res.HPLS.dn, "PsPS")
LPHSnetwork <- draw_network(res.LPHS.up, res.LPHS.dn, "pSPS")
LPLSnetwork <- draw_network(res.LPLS.up, res.LPLS.dn, "pspS")
LPLS2network <- draw_network(res.LPLS2.up, res.LPLS2.dn, "psPS")
LPHS2network <- draw_network(res.LPHS2.up, res.LPHS2.dn, "pSPs")
LPLS3network <- draw_network(res.LPLS3.up, res.LPLS3.dn, "psPs")

#####################WGCNA#####################
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = F)
# for draw_network outputs
# [[1]] Conversions between Arabidopsis and Eutrema merged with DEGS
# [[2]] Eutr genes that popped up in network that are not DEGs
# [[3]] The conversion between gene name and gene ID for Arabidopsis

# What's happening here:
# run WGCNA on the DEGs and plot them on the protein interaction network
# i.e. plotted as colour code


gsg <- goodSamplesGenes(FPKM, verbose = 3)
gsg$allOK
fpkm <- t(FPKM)
sampleTree <- hclust(dist(fpkm), method = "average") # 

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2)) # 

# already did it (chose 8)
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
network <- blockwiseModules(fpkm, power = 8,
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

#############################MAPMAN#######################
# Static P
forMM1 <- merge(forMMHPLS, forMMLPHS, all = T)
forMM <- merge(forMM1, forMMLPLS, all = T)

# Static S
forMMLPHSS
forMMLPLSS
forMMLPLS2S
forMM2 <- merge(forMMLPHSS, forMMLPLSS, all = T)
forMM3 <- merge(forMM2, forMMLPLS2S, all = T)

colnames(forMM3)[1] <- "IDENTIFIER"

write.table(forMM3, file = "~/R/Eutrema/MapMan/data/staticS.txt", quote = F, row.names = F, col.names = T, na = "0") 
