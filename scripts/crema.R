library(dplyr)
library(tibble)  # for `rownames <- to <- column` and `column <- to <- rownames`
library(methods)

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

