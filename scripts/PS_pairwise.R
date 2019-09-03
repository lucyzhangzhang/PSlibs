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
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]
#correct for batch effect
DESeq2Table$batch <- factor(DESeq2Table$reps, levels = c(1, 2, 3))

DESeq2Table$condition <- relevel(DESeq2Table$condition, ref = "HPHS")

dds <- DESeq(DESeq2Table)
resultsNames(dds)

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

write.table(rownames(HL.HH.pos), file = "names/HH/HL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(HL.HH.neg), file = "names/HH/HL/dn/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.HH.pos), file = "names/HH/LH/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.HH.neg), file = "names/HH/LH/dn/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HH.pos), file = "names/HH/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HH.neg), file = "names/HH/LL/dn/names", row.names = F, quote = F, col.names = F)
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

write.table(rownames(LL.LH.pos), file = "names/LH/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.LH.neg), file = "names/LH/LL/dn/names", row.names = F, quote = F, col.names = F)
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
write.table(data, "~/committee/DEGcount", sep = " & ", quote = F, eol = " \\\\\n", row.names = F)

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

write.table(rownames(res.HS.up), file = "names/HS/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(res.HS.dn), file = "names/HS/dn/names", row.names = F, quote = F, col.names = F)
#load annotation and FPKM files, generated previously with DESeq2, annotation file generated by Caitlin
expr_annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", row.names=1)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)

pos <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.up),]

neg <- expr_annot[rownames(expr_annot) %in% rownames(res.HS.dn),]

library(xlsx)
library(openxlsx)

write.xlsx2(pos, file = "Sulfurcomparison.xlsx", sheetName = "plusS", append = F, row.names = T)
write.xlsx2(neg, file = "Sulfurcomparison.xlsx", sheetName = "minusS", append = T, row.names = T)

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

write.table(rownames(HL.HH.pos), file = "names/HH/HL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(HL.HH.neg), file = "names/HH/HL/dn/names", row.names = F, quote = F, col.names = F)
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

#numbers
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),])
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),])
#names
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]

LH.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]
write.table(rownames(LH.HH.pos), file = "names/HH/LH/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.HH.neg), file = "names/HH/LH/dn/names", row.names = F, quote = F, col.names = F)

#group into matrix
HPHScompare <- matrix(c(LH.up, LH.dn), ncol = 2, byrow = T)
HPHScompare
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

res.LPLS <- results(dds, contrast = c("condition", "LPLS", "HPHS"))

#numbers
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])
#names
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

LL.HH.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LL.HH.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

write.table(rownames(LL.HH.pos), file = "names/HH/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HH.neg), file = "names/HH/LL/dn/names", row.names = F, quote = F, col.names = F)
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

#numbers
LH.up <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),])
LH.dn <- nrow(res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),])
#names
res.LPHS.up <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] > 0),]
res.LPHS.dn <- res.LPHS[which(res.LPHS[,"padj"] < 0.05 & res.LPHS[,"log2FoldChange"] < 0),]

LH.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.up),]
LH.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPHS.dn),]

write.table(rownames(LH.HL.pos), file = "names/HL/LH/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.HL.neg), file = "names/HL/LH/dn/names", row.names = F, quote = F, col.names = F)
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

#numbers
LL.up <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),])
LL.dn <- nrow(res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),])
#names
res.LPLS.up <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] > 0),]
res.LPLS.dn <- res.LPLS[which(res.LPLS[,"padj"] < 0.05 & res.LPLS[,"log2FoldChange"] < 0),]

LL.HL.pos <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.up),]
LL.HL.neg <- expr_annot[rownames(expr_annot) %in% rownames(res.LPLS.dn),]

#group into matrix
HPHScompare <- matrix(c(LL.up, LL.dn), ncol = 2, byrow = T)
HPHScompare
write.table(rownames(LL.HL.pos), file = "names/HL/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LL.HL.neg), file = "names/HL/LL/dn/names", row.names = F, quote = F, col.names = F)

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

write.table(rownames(LH.LL.pos), file = "names/LH/LL/up/names", row.names = F, quote = F, col.names = F)
write.table(rownames(LH.LL.neg), file = "names/LH/LL/dn/names", row.names = F, quote = F, col.names = F)

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


pos.combined <- as.data.frame(append(pos, list(PC2=pos.PC$PC2, PC4=pos.PC$PC4), after = 17))
neg.combined <- as.data.frame(append(neg, list(PC2=neg.PC$PC2, PC4=neg.PC$PC4), after = 17))

rownames(pos.combined) <- rownames(pos)
rownames(neg.combined) <- rownames(neg)

pos.combined <- pos.combined[order(pos.combined$PC2, decreasing = T),]
neg.combined <- neg.combined[order(neg.combined$PC2, decreasing = T),]

#write to excel file lmao

write.xlsx2(pos.combined, file = "ScomparisonPCs.xlsx", sheetName = "plusS", append = F, row.names = T)
write.xlsx2(neg.combined, file = "ScomparisonPCs.xlsx", sheetName = "minusS", append = T, row.names = T)
