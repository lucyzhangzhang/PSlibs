### PCA gene averages
library(openxlsx)
setwd("~/scratch/PS/caitlin/") #working directory goes here
#setwd("~/Documents/mcmaster/phd/rscripts/pca/pca_gene_avg/") # hard drive failure...getting scripts from Dropbox currently...please save this i guess??!?!?!
fpkm <- read.csv("~/scratch/PS/caitlin/20161230-fpkm.csv", row.names=1) #scp from McMaster cluster

#drought <- fpkm[,c(13:28,44:45)]
ps <- fpkm[,c(48:59)]

LPHS <- ps[,1:3]
HPLS <- ps[,4:6]
LPLS <- ps[,7:9]
HPHS <- ps[,10:12]

LPHS <- apply(LPHS, 1, median)
HPLS <- apply(HPLS, 1, median)
LPLS <- apply(LPLS, 1, median)
HPHS <- apply(HPHS, 1, median)

ps_med <- cbind(LPHS, HPLS, LPLS, HPHS)

ps.log <- log2(ps_med+1)

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

# ## Getting gene lists of "intereting" components
# # PC2 (neg and pos)
#
# # PC2
# posPC2 <- as.matrix(exprVals[order(exprVals$PC2, decreasing=TRUE)[1:50],2])
# posPC2Names <- rownames(exprVals[order(exprVals$PC2, decreasing=TRUE)[1:50],])
# rownames(posPC2) <- posPC2Names
#
# negPC2 <- as.matrix(exprVals[order(exprVals$PC2, decreasing=FALSE)[1:50],2])
# negPC2Names <- rownames(exprVals[order(exprVals$PC2, decreasing=FALSE)[1:50],])
# rownames(negPC2) <- negPC2Names
# # PC4
#
# posPC4 <- as.matrix(exprVals[order(exprVals$PC4, decreasing=TRUE)[1:50],4])
# posPC4Names <- rownames(exprVals[order(exprVals$PC4, decreasing=TRUE)[1:50],])
# rownames(posPC4) <- posPC4Names
#
# negPC4 <- as.matrix(exprVals[order(exprVals$PC4, decreasing=FALSE)[1:50],4])
# negPC4Names <- rownames(exprVals[order(exprVals$PC4, decreasing=FALSE)[1:50],])
# rownames(negPC4) <- negPC4Names
#
expr_annot <- read.csv("~/Eutrema/FPKM/2016-12-30_fpkm_annot.csv", row.names=1)
#
# neg2genes <- match(rownames(negPC2), rownames(expr_annot))
# pos2genes <- match(rownames(posPC2), rownames(expr_annot))
# neg4genes <- match(rownames(negPC4), rownames(expr_annot))
# pos4genes <- match(rownames(posPC4), rownames(expr_annot))
#
#
# neg2 <- cbind(negPC2, expr_annot[neg2genes,7:89])
# pos2 <- cbind(posPC2, expr_annot[pos2genes,7:89])
# neg4 <- cbind(negPC4, expr_annot[neg4genes,7:89])
# pos4 <- cbind(posPC4, expr_annot[pos4genes,7:89])
#
# #write.csv(neg2, file="PC2_negative.csv")
# #write.csv(pos2, file="PC2_positive.csv")
# #write.csv(neg4, file="PC4_negative.csv")
# #write.csv(pos4, file="PC4_positive.csv")
#
# #####################
# # Want Pos/Neg PC2, and around 0 PC4
# #####################
#
# # < 0 PC2, <0.1 PC4
# # < 0 PC2, >0.1 PC4
#
# #positive
# pos2_specific <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] < 0.1 & exprVals[,4] > -0.1),c(2,4)]
# pos2genes_specific <- expr_annot[match(rownames(pos2_specific), rownames(expr_annot)),]
# pos2genes_expr <- cbind(pos2_specific, pos2genes_specific[,-c(1:6)])
# pos2genes_top <- pos2genes_expr[order(pos2genes_expr$PC2, decreasing=TRUE)[1:50],]
# #negative
# neg2_specific <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] < 0.1 & exprVals[,4] > -0.1),c(2,4)]
# neg2genes_specific <- expr_annot[match(rownames(neg2_specific), rownames(expr_annot)),]
# neg2genes_expr <- cbind(neg2_specific, neg2genes_specific[,-c(1:6)])
# neg2genes_top <- neg2genes_expr[order(neg2genes_expr$PC2, decreasing=TRUE)[1:50],]
#
# #write.csv(neg2genes_top, file="PC2_negative-PC4_origin.csv")
# #write.csv(pos2genes_top, file="PC2_positive-PC4_origin.csv")
#
# exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,]

####################
# Pos PC2, Neg PC4 #
####################
# top 50?

pos2neg4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] < 0), c(2,4)]

pos2_200 <-pos2neg4[order(pos2neg4$PC2, decreasing=TRUE)[1:200],]
neg4_200 <- pos2neg4[order(pos2neg4$PC4, decreasing = FALSE)[1:200],]

pos2neg4_58 <- pos2_200[rownames(pos2_200) %in% rownames(neg4_200),]
#pos2neg4_58 <- neg4_200[na.omit(match(rownames(pos2_200), rownames(neg4_200))),] # locations in pos2_200
#na.omit(match(rownames(pos2neg4_58), rownames(expr_annot)))
#pos2neg4_genes <- expr_annot[na.omit(match(rownames(pos2neg4_58), rownames(expr_annot))),]
pos2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2neg4_58),]
#pos2neg4_genes <- cbind(pos2neg4_58, pos2neg4_genes[,7:ncol(pos2neg4_genes)])
pos2neg4_genes <- merge(pos2neg4_58, pos2neg4_genes[,7:ncol(pos2neg4_genes)], by=0)
rownames(pos2neg4_genes) <- pos2neg4_genes$Row.names
pos2neg4_genes <- subset(pos2neg4_genes, select = -c(Row.names))
#write.csv(pos2neg4_genes, file="pos2neg4_genes.csv")

####################
# Pos PC2, Pos PC4 #
####################
# top 50?

pos2pos4 <- exprVals[which(exprVals[,2] > 0 & exprVals[,4] > 0), c(2,4)]

pos2_200 <-pos2pos4[order(pos2pos4$PC2, decreasing=TRUE)[1:300],]
pos4_200 <- pos2pos4[order(pos2pos4$PC4, decreasing = TRUE)[1:300],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#pos2pos4_top <- pos2_200[na.omit(match(rownames(pos4_200), rownames(pos2_200))),] # 33
pos2pos4_top <- pos2_200[rownames(pos2_200) %in% rownames(pos4_200),]
pos2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(pos2pos4_top),]
#pos2pos4_genes <- cbind(pos2pos4_top, pos2pos4_genes[,7:ncol(pos2pos4_genes)])
pos2pos4_genes <- merge(pos2pos4_top, pos2pos4_genes[,7:ncol(pos2pos4_genes)], by=0)
rownames(pos2pos4_genes) <- pos2pos4_genes$Row.names
pos2pos4_genes <- subset(pos2pos4_genes, select = -c(Row.names))

write.csv(pos2pos4_genes, file="pos2pos4_genes.csv")
####################
# Neg PC2, Pos PC4 #
####################
# top 50?

neg2pos4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] > 0), c(2,4)]

neg2_200 <-neg2pos4[order(neg2pos4$PC2, decreasing=FALSE)[1:200],]
pos4_200 <- neg2pos4[order(neg2pos4$PC4, decreasing = TRUE)[1:200],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#neg2pos4_top <- neg2_200[na.omit(match(rownames(pos4_200), rownames(neg2_200))),] # 78
neg2pos4_top <- neg2_200[rownames(neg2_200) %in% rownames(pos4_200),]

neg2pos4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2pos4_top),]
#neg2pos4_genes <- cbind(neg2pos4_top, neg2pos4_genes[,7:ncol(neg2pos4_genes)])
neg2pos4_genes <- merge(neg2pos4_top, neg2pos4_genes[,7:ncol(neg2pos4_genes)], by=0)
row.names(neg2pos4_genes) <- neg2pos4_genes$Row.names
neg2pos4_genes <- subset(neg2pos4_genes, select = -c(Row.names))

#write.csv(neg2pos4_genes, file="neg2pos4_genes.csv")
####################
# Neg PC2, Neg PC4 #
####################
# top 50?

neg2neg4 <- exprVals[which(exprVals[,2] < 0 & exprVals[,4] < 0), c(2,4)]

neg2_300 <-neg2neg4[order(neg2neg4$PC2, decreasing=FALSE)[1:300],]
neg4_300 <- neg2neg4[order(neg2neg4$PC4, decreasing = FALSE)[1:300],]
#na.omit(match(rownames(pos2_200), rownames(neg4_200)))
#neg2neg4_top <- neg2_200[na.omit(match(rownames(neg4_200), rownames(neg2_200))),] # 78
neg2neg4_top <- neg2_300[rownames(neg2_300) %in% rownames(neg4_300),]

# neg2_200[rownames(neg2_200) %in% rownames(pos4_200),]
neg2neg4_genes <- expr_annot[rownames(expr_annot) %in% rownames(neg2neg4_top),]
#neg2neg4_genes <- cbind(neg2neg4_top, neg2neg4_genes[,7:ncol(neg2neg4_genes)])
neg2neg4_genes <- merge(neg2neg4_top, neg2neg4_genes[,7:ncol(neg2neg4_genes)], by=0)
row.names(neg2neg4_genes) <- neg2neg4_genes$Row.names
neg2neg4_genes <- subset(neg2neg4_genes, select = -c(Row.names))


#write.csv(neg2neg4_genes, file="neg2neg4_genes.csv")

##################
# Write to excel #
##################

pca_list <- list(neg2neg4_genes,  neg2pos4_genes, pos2neg4_genes, pos2pos4_genes)
names(pca_list) = c("neg2neg4",  "neg2pos4", "pos2neg4", "pos2pos4")
#write.xlsx(pca_list, file = "pca_genelists_updated2.xlsx", rowNames=TRUE)
#####
## Visualizing PCAs...

#samples <-colnames(ps_med)
samples <- c("pS", "Ps", "ps", "PS")
coords <- data.frame(X=rep(0, 4), Y=rep(0, 4),sampleVals, Samples = samples)
coords$Samples2 <- factor(coords$Samples, c("PS", "Ps", "pS", "ps"))

library(ggplot2)
theme_set(theme_bw())

labels_pos2neg4 <- ifelse(rownames(exprVals) %in% rownames(pos2neg4_58), "Yes", "No")
labels_pos2pos4 <- ifelse(rownames(exprVals) %in% rownames(pos2pos4_top), "Yes", "No")
labels_neg2neg4 <- ifelse(rownames(exprVals) %in% rownames(neg2neg4_top), "Yes", "No")
exprVals_lab <- cbind(exprVals, labels_pos2neg4, labels_pos2pos4, labels_neg2neg4)

### Labelling plot!
(pc24plot_labelled <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
  geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples2), arrow=arrow(), size=1) +
  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
  geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)))
#png("figs/pc24_labelled_recoloured.png", bg="transparent")
#print(pc24plot_labelled)
#dev.off()

(pc24plot_pos2neg4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
 # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
 # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
 # geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
  geom_point(data = exprVals_lab[exprVals_lab$labels_pos2neg4 == "Yes",], colour = "mediumpurple1", alpha = 0.5) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size =1) +
 # geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
 # geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
 # geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16)))

(pc24plot_pos2pos4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
#  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
#  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
#  geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
  geom_point(data = exprVals_lab[exprVals_lab$labels_pos2pos4 == "Yes",], colour = "mediumpurple1", alpha = 0.7) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
#  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
#  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
#  geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16)))

(pc24plot_neg2neg4lab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
  #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="red") +
  #  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="red") +
  #  geom_point(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], colour="red") +
  geom_point(data = exprVals_lab[exprVals_lab$labels_neg2neg4 == "Yes",], colour = "mediumpurple1", alpha = 0.4) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
  #  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "At4"), vjust=1.5, colour="white") +
  #  geom_text(data = exprVals[rownames(exprVals) == "XLOC_003912" ,], aes(PC2,PC4, label = "XLOC"), vjust=-0.5, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16)))

#ggsave("pos2pos4_labelled.pdf", pc24plot_pos2pos4lab)
#ggsave("pos2neg4_labelled.pdf", pc24plot_pos2neg4lab)
#ggsave("neg2neg4_labelled.pdf", pc24plot_neg2neg4lab)



(pc24plot_recol <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="#3366FF", size=4) +
 # geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=2) +
  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#3366FF", size=4) +
  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), vjust=-1, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
#scale_colour_manual(values=c("#fd8a19", "#fcfc35", "#f91301", "#febfcb"))
#scale_colour_manual(values=c( "#febfcb", "#fcfc35", "#f91301", "#fd8a19" ))
  #scale_colour_manual(values=c( "coral2", "#fcfc35", "#f91301", "#fd8a19" ))
scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16)))
#ggsave(pc24plot_recol, file="pc24_recoloured.pdf")

pc24plot_nolab <- ggplot(exprVals, aes_string("PC2", "PC4")) +
  geom_point(shape=19, alpha=0.3) +
 # geom_point(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], colour="#3366FF", size=4) +
 geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=2) +
#  geom_point(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], colour="#3366FF", size=4) +
#  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10004656m.g" ,], aes(PC2,PC4, label = "SDI1"), vjust=-1, colour="black") +
#  geom_text(data = exprVals[rownames(exprVals) == "Thhalv10015137m.g" ,], aes(PC2,PC4, label = "IPS2"), vjust=-1, colour="white") +
  xlab("PC2 (0.05%)") + ylab("PC4 (0.02%)") +
  coord_cartesian(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) +
  #scale_colour_manual(values=c("#fd8a19", "#fcfc35", "#f91301", "#febfcb"))
  #scale_colour_manual(values=c( "#febfcb", "#fcfc35", "#f91301", "#fd8a19" ))
  #scale_colour_manual(values=c( "coral2", "#fcfc35", "#f91301", "#fd8a19" ))
  scale_colour_manual(values=c( "#febfcb", "gold", "#f91301", "#fd8a19" ), name = "") +
  theme_set(theme_bw(base_size=16))
ggsave(pc24plot_nolab, file="pc24_recoloured_noLabels.pdf")





## loop for all PC1 combinations
#for (tryThis in 2:4){
 # toUse<-noquote(paste("PC", tryThis, sep=""))
  #fileName<-paste("figs/PC1",toUse, ".pdf",sep="")
  #pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC1", toUse)) +
         # geom_point(aes(color = Samples), size=5) +
          geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=toUse, colour=Samples), arrow=arrow(), size=1))# +
         # coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)))
  #dev.off()
#}

pc12plot <- ggplot(coords, aes_string("PC1", "PC2")) +
          geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC2, colour=Samples), arrow=arrow(), size=1) +
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))


pc13plot <- ggplot(coords, aes_string("PC1", "PC3")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC3, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))

pc14plot <- ggplot(coords, aes_string("PC1", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC1, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))

pc23plot <- ggplot(coords, aes_string("PC2", "PC3")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC3, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))

pc24plot <- ggplot(coords, aes_string("PC2", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC2, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))






pc34plot <- ggplot(coords, aes_string("PC3", "PC4")) +
  geom_segment(data=coords, aes(x=X, y=Y, xend=PC3, yend=PC4, colour=Samples), arrow=arrow(), size=1) +
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8))

#pdf("figs/PC1PC2.pdf")
#print(pc12plot)
#dev.off()
#
#pdf("figs/PC1PC3.pdf")
#print(pc13plot)
#dev.off()
#
#pdf("figs/PC1PC4.pdf")
#print(pc14plot)
#dev.off()
#
#pdf("figs/PC2PC3.pdf")
#print(pc23plot)
#dev.off()
#
#pdf("figs/PC2PC4.pdf")
#print(pc24plot)
#dev.off()
#
#pdf("figs/PC3PC4.pdf")
#print(pc34plot)
#dev.off()
