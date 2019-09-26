library(dplyr)
library(tibble)  # for `rownames <- to <- column` and `column <- to <- rownames`
library(methods) # to make the S4 object

wd  <- "~/R/Eutrema/PS/names"
setwd(wd)

# being compared to
Denom <- c(rep("HH", 6), rep("HL", 4), rep("LH", 2), rep("HS", 2))
Dlabel <- c(rep("PS", 6), rep("Ps", 4), rep("pS", 2), rep("S", 2))
Expr <- rep(c("up", "dn"), length(Denom)/2)
# the comparison
Num <- rep(c("HL", "LH", "LL", "LH", "LL", "LL", ""), each = 2)
Nlabel <- rep(c("Ps", "pS", "ps", "pS", "ps", "ps", "s"), each = 2)

Meta <- data.frame(Dlabel = Dlabel, Expr = Expr, Nlabel = Nlabel)

# add uniqe IDs
Meta <- Meta %>% mutate(ID = paste(Nlabel, Dlabel, Expr, sep = "_"))

# move last column to the beginning of the table
Meta <- Meta[, c(ncol(Meta), 1:(ncol(Meta) - 1))]

# add file names
Meta <- Meta %>% mutate(file = paste(wd, Denom, Num, Expr, "final_ensemble_predictions.csv", sep = "/"))
Meta <- Meta %>% mutate(names = paste(wd, Denom, Num, Expr, "names", sep = "/"))

# count the number of genes in each group
# print 0 if the RNA prediction file doesn't exist
Meta <- Meta %>% mutate(Genes = sapply(names, function(x) {
                                           tryCatch(nrow(read.table(x)), error = function(e) {return(0)})})) 

# how many are predicted as lncRNA by crema?
# print 0 if the RNA prediction file doesn't exist
Meta <- Meta %>% mutate(pred = sapply(file, function(x) {
                                          ifelse(!file.exists(x), 0, 
                                                 length(which(read.csv(x, header = T)$prediction == 1)))}))

# annotation file
expr.annot <- read.csv("~/R/Eutrema/PS/FPKMS_LengthAdjusted.csv", header = T, row.names=1)

# new prediction class
predictions <- setClass("predictions", slots = c(ID = "character", Dlabel = "factor",
                                                 Nlabel = "factor", Expr = "factor",
                                                 file = "character", names = "character",
                                                 Genes = "numeric", pred = "numeric",
                                                 predAnnot = "data.frame", annot = "data.frame"))

# function to generated list of new class "predictions"
newPred <- function(Meta, annot) {
    # initialize output list
    out <- list()

    # iterate every row in metadata
    for (i in 1:nrow(Meta)) {
        metaRow <- Meta[i,]

        # find annotations of the predicted lncRNAs
        predNames <- tryCatch(read.csv(metaRow$file, header = T, row.names = 1), 
                              error = function(e) {return(NA)})
        # don't use ifelse(), that shit gay
        if(class(predNames) == "data.frame"){
            predNames <- predNames %>% 
                rownames_to_column("genes") %>% 
                filter(prediction == 1) %>%
                column_to_rownames("genes")
            predAnnot <- annot[tolower(rownames(annot)) %in% rownames(predNames),]
        }else{
            # force data frame because class restrictions
            predNames <- data.frame(NA)
            predAnnot <- data.frame(NA)
        }
        
        # find annotations of the DEGs
        names <- tryCatch(scan(metaRow$names, what = "character", quiet = T), 
                          error = function(e) {return(NA)})
        Annot <- annot[rownames(annot) %in% names,]
        
        # create new predictions S4 object
        entry <- predictions(ID = metaRow$ID, Dlabel = metaRow$Dlabel,
                             Nlabel = metaRow$Nlabel, file = metaRow$file,
                             names = metaRow$names, Genes = metaRow$Genes,
                             pred = metaRow$pred, annot = Annot,
                             predAnnot = predAnnot)

        # add to output list
        out <- append(out, entry)
    }

    return(out)
}

# generate a new prediction class object useing precious function
pred <- suppressWarnings(newPred(Meta, expr.annot))

# view annotations of all the predicted lncRNAs
predAnnots <- c()
for (i in 1:length(pred)) {
    add <- list(pred[[i]]@predAnnot)
    names(add) <- pred[[i]]@ID
    predAnnots <- append(predAnnots, add)
}

head(predAnnots)
test <- data.frame(Genes="No predicted genes")
library(xlsx)
library(openxlsx)
{lapply(names(predAnnots), function(name) if (is.na(predAnnots[[name]]) || (nrow(predAnnots[[name]]) == 0)) {

           write.xlsx2(test, file = "predictions.xlsx",
                       sheetName = name, append = T)

} else if (which(names(predAnnots) == name) == 1) {
           write.xlsx2(data.frame(predAnnots[[name]]), file = "predictions.xlsx",  
                       sheetName = name, append = F)

} else {
           write.xlsx2(data.frame(predAnnots[[name]]), file = "predictions.xlsx", 
                       sheetName = name, append = T)}
)}

{lapply(names(predAnnots), function(name) write.table(predAnnots[[name]], file = paste0(name, ".tab"), sep = ",", quote = T))}

# for (name in names(predAnnots)) {
#     write.xlsx2(predAnnots[[name]], file = "predictions.xlsx", sheetName = name, append = T)
# }

################################VISUALIZATION##########################
library(reshape2)
library(ggplot2)
library(ggrepel)
# predAnnots[[11]]$Ath.desc
data <- Meta[, c("Dlabel", "Nlabel", "Expr", "Genes", "pred")]
data <- data %>% mutate(ID = paste(Nlabel, Dlabel, sep = ".vs."))
data <- data[, !(names(data) %in% c("Dlabel", "Nlabel"))]   
data$Genes <- data$Genes - data$pred

label  <- c(Meta[, "Genes"], Meta[, "pred"])
names(label)  <- NULL
data.melt <- melt(data, id.vars = c("ID", "Expr"), measure.vars = c("Genes", "pred"),
value.name = "Genes")


data.melt$Genes[data.melt$Expr == "dn"]  <- data.melt$Genes[data.melt$Expr == "dn"] * -1

data.melt <- cbind(data.melt, label)

data.melt$Expr  <- factor(data.melt$Expr, levels = c("up", "dn"))
data.melt$label[data.melt$Expr == "dn"]  <- data.melt$label[data.melt$Expr == "dn"] * -1
data.melt <- mutate(data.melt, grouping = rep(rep(c(rep("PS", 3), rep("Ps", 2), "pS", "S"), each = 2), 2))
data.melt <- mutate(data.melt, Nlabel = rep(Nlabel, 2))

#hack in P data
pdat <- data.frame(ID="p.vs.P", Expr="up", variable="pred", 
                                Genes=as.numeric(1), label=as.numeric(1),grouping="P", Nlabel="p", stringsAsFactors = F) 
data.melt <- rbind.data.frame(data.melt, pdat, stringsAsFactors = F) 

# plotting
plot <- {ggplot(data.melt, aes(ID, Genes, alpha = variable, fill = Expr)) + # 
    geom_bar(stat = "identity") + theme_bw() +
    labs(x = "Comparison", y = "Number of transcripts") +
    scale_fill_discrete(name = "Expression", labels = c("Upregulated", "Downregulated"),
                        breaks = c("up", "dn")) +
scale_x_discrete(breaks = data.melt$ID, labels = data.melt$Nlabel) + 
facet_grid(. ~grouping, scales = "free_x", space = "free_x") + 
scale_alpha_discrete(name = "Category", labels = c("Transcripts", "Predicted lncRNA"),
                     range = c(0.3,1)) +
geom_text_repel(data = data.melt[data.melt$variable == "Genes", ], 
                aes(x = ID, y = label, label = abs(Genes), color = Expr), 
          nudge_y = c(5, -5), alpha = 1, show.legend = F, seed = 7, point.padding = NA, direction = "y") +
geom_text(data = data.melt[(data.melt$variable == "pred") & (data.melt$label != 0) & (data.melt$Genes > 0), ], 
          aes(x = ID, y = label, label = abs(Genes)), 
                 alpha = 1, show.legend = F, color = "black", nudge_y = 5) +
geom_text(data = data.melt[(data.melt$variable == "pred") & (data.melt$label != 0) & (data.melt$Genes < 0), ], 
          aes(x = ID, y = label, label = abs(Genes)), 
                 alpha = 1, show.legend = F, color = "black", nudge_y = -5) +
    theme(panel.spacing = unit(0, "lines"))}
    plot

ggsave("predict.pdf", plot, dpi = 300, height = 7, width = 8)

up <- data %>% filter(Expr == "up")
dn <- data %>% filter(Expr == "dn") %>% mutate(Genes = Genes * -1, pred = pred * -1)
write.table(cbind(num = rownames(up), up), file = "up", row.names = F, quote = F, col.names = T, sep = "\t")
write.table(cbind(num = rownames(dn), dn), file = "dn", row.names = F, quote = F, col.names = T, sep = "\t")

#############################PCA loading predictions###############################
 setwd("~/R/Eutrema/PS/crema/pcaLoading")

library(ggplot2)
FPKM <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", header = T, row.names = 1)
sampNames <- c("ps1", "ps2", "ps3", "Ps1", "Ps2", "Ps3", "pS1", "pS2", "pS3", "PS1", "PS2", "PS3")
colnames(FPKM) <- sampNames

p2.raw <- read.csv("~/R/Eutrema/PS/crema/pcaLoading/p2/final_ensemble_predictions.csv", header = T, row.names = 1)
n2.raw <- read.csv("~/R/Eutrema/PS/crema/pcaLoading/n2/final_ensemble_predictions.csv", header = T, row.names = 1)
p4.raw <- read.csv("~/R/Eutrema/PS/crema/pcaLoading/p4/final_ensemble_predictions.csv", header = T, row.names = 1)
n4.raw <- read.csv("~/R/Eutrema/PS/crema/pcaLoading/n4/final_ensemble_predictions.csv", header = T, row.names = 1)

p2.predname <- rownames(p2.raw[which(p2.raw$prediction == 1),])
n2.predname <- rownames(n2.raw[which(n2.raw$prediction == 1),])
p4.predname <- rownames(p4.raw[which(p4.raw$prediction == 1),])
n4.predname <- rownames(n4.raw[which(n4.raw$prediction == 1),])

pred.names <- union(union(p2.predname, n2.predname),union(p4.predname, n4.predname))
write(pred.names, "~/R/Eutrema/PS/crema/pcaLoading/pred.names", sep = "\n")
fpkm <- read.table("~/R/Eutrema/PS/FPKMS_LengthAdjusted2.tab", row.names=1, header = T) #scp from McMaster cluster

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

## all genes
pca<- prcomp(ps_mean, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here

sdev <- pca$sdev

screeplot(pca)
#log.expr <- log2(expr+1)

plot(pca, type = "l")
summary(pca)
exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)

plot.points <- function(points, data){
    subset <- data[which(toupper(rownames(data)) %in% toupper(points)),]
    plot <- ggplot(subset, aes_string("PC2", "PC4"))+
        geom_point(shape = 19, size = 2)
}

test <- plot.points(pred.names, exprVals)
samples <- c("ps", "pS", "Ps", "PS")
coords <- data.frame(X=rep(0, 4), Y=rep(0, 4),sampleVals, Samples = samples)
coords$Treatment <- factor(coords$Samples, c("ps", "pS", "Ps", "PS"))
pcaSum <- as.data.frame(summary(pca)$importance)
(mean <- ggplot(exprVals, aes_string("PC2", "PC4")) +
    geom_point(shape=19, alpha=0.3) +
    geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Treatment), arrow=arrow(length = unit(0.3, "cm"), angle = 45), size = 1.5) +
    geom_point(data = exprVals[which(toupper(rownames(exprVals)) %in% toupper(pred.names)) ,], colour="#2851b5", size=2) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), nudge_y = 0.07) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10018244m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10018244m.g" ,], aes(PC2,PC4, label = "SULTR1;2"), nudge_y = 0.07) +
    #     geom_point(data = exprVals[rownames(exprVals) == "Thhalv10005068m.g" ,], colour="#2851b5", size=2) +
    #     geom_shadowtext(data = exprVals[rownames(exprVals) == "Thhalv10005068m.g" ,], aes(PC2,PC4, label = "Rhodanese"), nudge_y = 0.07) +
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
    xlab(paste0("PC2 ", pcaSum$PC2[2]*100,"%")) + ylab(paste0("PC4 ",pcaSum$PC4[2]*100,"%")) +
    coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) )
mean
ggsave("~/R/Eutrema/PS/PSMeanPred1Up.pdf", mean, dpi = 250, height = 6, width = 8)

#############ALL GENES##################


ps_mean <- data.frame(ps, pS, Ps, PS)
ps_filt <- ps_mean[(ps_mean$ps >= 1) | (ps_mean$Ps >= 1) | (ps_mean$pS >= 1) | (ps_mean$PS >= 1),]
colnames(ps_filt) <- colnames(ps_mean)
psn <- rownames(ps.filt[which(ps.filt$ps != 0),])
Psn <- rownames(ps.filt[which(ps.filt$Ps != 0),])
pSn <- rownames(ps.filt[which(ps.filt$pS != 0),])
PSn <- rownames(ps.filt[which(ps.filt$PS != 0),])

allGN <- union(union(psn, PSn), union(pSn, Psn))
write(allGN, "~/R/Eutrema/PS/crema/allG/allG.names", sep = "\n")
