library(ggplot2)
library(ggrepel)
library(viridis)

library(DESeq2)
library(readr)
library(dplyr)
#setwd("/Users/caitlinsimopoulos/Dropbox/McMaster/plant_biotech/ireland/talk/figs/") #working directory goes here
# setwd("/Users/caitlinsimopoulos/Dropbox/McMaster/PhD/Written/Manuscripts/drought/figs")
sedwd("~/scratch/misc")

fpkm <- read.csv("~/Eutrema/FPKM/2018-10-15-drought_geneFPKM.csv", sep=",") ## isoform level

unique_fpkm <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]

rownames(unique_fpkm) <- unique_fpkm$gene_id

unique_fpkm <- unique_fpkm[,-1]
y.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"y_lnc", drop=F] == 1),])
s.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"s_lnc", drop=F] == 1),])
rep_lib <- unique_fpkm[,c(9:15, 17:22, 23:25, 29:35, 37:43, 44:46)]
fpkmunique <- unique_fpkm[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
colnames(fpkmunique)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")


#fpkmunique <- fpkmunique[apply(fpkmunique > 1, 1, any),]
dim(fpkmunique)

fpkmlog <- log2(fpkmunique+1)

pca <- prcomp(t(fpkmlog), center=T, scale=FALSE)
plot(pca, type = "l")
x.var <- pca$sdev ^ 2
x.pvar <- x.var/sum(x.var)

plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')

summary(pca)
#summary(pca.top)

colnames(rep_lib)[c(1,2,8,9,17,18,24,25)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")
screeplot(pca)
pcarep <- prcomp(log2(rep_lib+1), scale=F)
exprVals<-data.frame(pcarep$x)
sampleVals<-data.frame(pcarep$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(rep_lib),1,nchar(colnames(rep_lib))-2)
ecotype <- c(rep("Yukon", 16), rep("Shandong", 17))


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 33), Y=rep(0, 33), sampleVals,
                   Samples = factor(samples), Name = colnames(rep_lib),
                   Condition = factor(substr(samples, 2,nchar(samples))),
                   Ecotype = ecotype)
coords$Condition[c(13,30)] <- c("D2","D2")
#coords$Condition[c(8, 18, 19, 20)] <- c("D1", "D2" , "WW1", "WW1")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (pcarep$sdev^2/sum(pcarep$sdev^2)) * 100
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

#coords$EcotypeLib<- c("Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B",
#                       "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B",
#                       "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B",
#                       "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-B")
(pc24plot_rep<-ggplot(exprVals, aes(PC2, PC4)) +
    geom_point(data = coords, aes(x=PC2, y=PC4, fill = Condition, shape = Ecotype),
               size=5) +  #, show.legend=F) +
    coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC2", PC2per)) +
    scale_y_continuous(name = paste("PC4", PC4per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5,size=4) +
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))
#ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/drought/figs/pc24plot_rep.pdf", plot=pc24plot_rep,
       # width=7, height=5, units="in")

########
# FIGS #
########

exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(fpkmlog),1,nchar(colnames(fpkmlog))-2)
ecotype <- c(rep("Yukon", 15), rep("Shandong", 16))


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 31), Y=rep(0, 31), exprVals,
                   Samples = factor(samples), Name = colnames(fpkmlog),
                   Condition = factor(substr(samples, 2,nchar(samples))),
                   Ecotype = ecotype)

#coords$Condition[c(8, 18, 19, 20)] <- c("D1", "D2" , "WW1", "WW1")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (pca$sdev^2/sum(pca$sdev^2)) * 100
PC1per <- paste0("(", round(PoV[1], 2), "%", ")")
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

coords$EcotypeLib <- c("Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B",
                           "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B",
                           "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B",
                           "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-B")
(pc12plot<-ggplot(exprVals, aes(PC1, PC2)) +
    geom_point(data = coords, aes(x=PC1, y=PC2, fill = Condition, shape = EcotypeLib),
               size=5) +  #, show.legend=F) +
    #coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC1", PC1per)) +
    scale_y_continuous(name = paste("PC2", PC2per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC1,PC2,label=Name), force=1, point.padding = 0.5,size=4) +
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,24, 22,25)) +
    theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))

(pc13plot_new<-ggplot(exprVals, aes(PC1, PC3)) +
    geom_point(data = coords, aes(x=PC1, y=PC3, fill = Condition, shape = Ecotype),
               size=5) +  #, show.legend=F) +
    # coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC1", PC1per)) +
    scale_y_continuous(name = paste("PC3", PC3per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC1,PC3,label=Name), force=1, point.padding = 0.5,size=4) +
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18),
          panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))

(pc23plot<-ggplot(exprVals, aes(PC2, PC3)) +
  geom_point(data = coords, aes(x=PC2, y=PC3, fill = Condition, shape = EcotypeLib),
             size=5) +  #, show.legend=F) +
  #coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC3", PC3per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC3,label=Name), force=1, point.padding = 0.5,size=4) +
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
  theme_classic() +
  scale_shape_manual(values=c(21,24, 22,25)) +
  theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21)))) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/drought/figs/pc23plot.pdf", plot=pc23plot, width=7, height=5, units="in")
#ggsave("yukpc23plot.pdf", plot=pc23plot, width=9, height=7, units="in")
#viridis::scale_color_viridis(discrete = TRUE, option = "C") +
#scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +

(pc24plot<-ggplot(exprVals, aes(PC2, PC4)) +
    geom_point(data = coords, aes(x=PC2, y=PC4, fill = Condition, shape = Ecotype),
               size=5) +  #, show.legend=F) +
   # coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC2", PC2per)) +
    scale_y_continuous(name = paste("PC4", PC4per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5,size=4) +
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18),
          panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))
ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/drought/figs/pc24plot.pdf", plot=pc24plot, width=9, height=7, units="in")
grey-green <- c("#88AB75") #if you need your old color

 PC34plot<-ggplot(exprVals, aes(PC3, PC4)) +
  #geom_point(data = coords, aes(x=PC3, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC4,label=Name), force=1, point.padding = 0.5, size=4) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15))






## lets use rlog, blind=T instrad of log2 for now
# yuk <- drought[,1:20]
# yuk_f <- drought[,c(1:20, 41:43)]
# yuk_f.log <- log2(yuk_f+1)
# shan <- drought[,21:40]
# yuk.log <- log2(yuk+1)
# shan.log <- log2(shan+1)
yuk <- rep_lib[,1:20]
yuk_f <- rep_lib[,c(1:20, 41:43)]
yuk_f.log <- log2(yuk_f+1)
shan <- rep_lib[,21:40]
yuk.log <- log2(yuk+1)
shan.log <- log2(shan+1)

pca.yuk<- prcomp(yuk.log, scale=FALSE)  # PCA with centering and scaling
pca.yuk$rotation  # The loadings are here


plot(pca.yuk, type = "l")
summary(pca.yuk)
#summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(pca.yuk$x)
sampleVals<-data.frame(pca.yuk$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(yuk),1,nchar(colnames(yuk))-2)
samples[c(8, 14,18)] <- c("YD1w", "YD2", "YD2w")
#samples <- c(samples, rep("Y-F2003", 2))
ecotype <- c(rep("Yukon", 20))
#cond <- substring(colnames(drought)[1:16],2,3)
#cond <- c(cond, rep("field", 2))
label <- c(rep("yes", 43))


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 20), Y=rep(0, 20), sampleVals,
                   Samples = samples, Name = colnames(yuk),
                   Condition = factor(substr(samples, 2,nchar(samples))),
                   Ecotype = ecotype)

coords$Condition[c(8, 18, 19, 20)] <- c("D1", "D2" , "WW1", "WW1")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (pca.yuk$sdev^2/sum(pca.yuk$sdev^2)) * 100
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")
#theme_set(theme_bw()) # make background white
library(ggplot2)
library(ggrepel)
pc23plot<-ggplot(exprVals, aes(PC2, PC3)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC3, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC3, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC3", PC3per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC3,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
pc23plot
ggsave("yukpc23plot.pdf", plot=pc23plot, width=9, height=7, units="in")
#viridis::scale_color_viridis(discrete = TRUE, option = "C") +
#scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +


pc24plot<-ggplot(exprVals, aes(PC2, PC4)) +
  #geom_point(data = coords, aes(x=PC2, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5, size=4) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15))
pc24plot
PC34plot<-ggplot(exprVals, aes(PC3, PC4)) +
  #geom_point(data = coords, aes(x=PC3, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC4,label=Name), force=1, point.padding = 0.5, size=4) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15))


pc2pc5 <- ggplot(coords, aes(PC2, PC5)) +
  geom_point(data=coords, aes(x=PC2, y=PC5, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC5, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC5", PC5per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC5,label=Name),force=1, point.padding = 0.2) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic()

pc2pc6 <- ggplot(coords, aes(PC2, PC6)) +
  geom_point(data=coords, aes(x=PC2, y=PC6, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC6, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC6", PC6per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC6,label=Name),force=1, point.padding = 0.2) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic()


### shandong

pca.shan<- prcomp(shan.log, scale=FALSE)  # PCA with centering and scaling
pca.shan$rotation  # The loadings are here

plot(pca.shan, type = "l")
summary(pca.shan)
#summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(pca.shan$x)
sampleVals<-data.frame(pca.shan$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(shan),1,nchar(colnames(shan))-2)
samples[c(8, 15)] <- c("SD1w", "SD2")
#samples <- c(samples, rep("Y-F2003", 2))
ecotype <- c(rep("Shandong", 20))
#cond <- substring(colnames(drought)[1:16],2,3)
#cond <- c(cond, rep("field", 2))
label <- c(rep("yes", 43))


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 20), Y=rep(0, 20), sampleVals,
                   Samples = samples, Name = colnames(shan),
                   Condition = factor(substr(samples, 2,nchar(samples))),
                   Ecotype = ecotype)

coords$Condition[c(8,19,20)] <- c("D1", "WW1", "WW1" )
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (pca.shan$sdev^2/sum(pca.shan$sdev^2)) * 100
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")
#theme_set(theme_bw()) # make background white
library(ggplot2)
library(ggrepel)
pc23plot<-ggplot(exprVals, aes(PC2, PC3)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC3, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC3, color = Condition, shape = Ecotype),
             size=4) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC3", PC3per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC3,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
#viridis::scale_color_viridis(discrete = TRUE, option = "C") +
#scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
ggsave("shanpc23plot.pdf", plot=pc23plot, width=9, height=7, units="in")

pc24plot<-ggplot(exprVals, aes(PC2, PC4)) +
  #geom_point(data = coords, aes(x=PC2, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5, size=4) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15))

PC34plot<-ggplot(exprVals, aes(PC3, PC4)) +
  #geom_point(data = coords, aes(x=PC3, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC4,label=Name), force=1, point.padding = 0.5, size=4) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15))


pc2pc5 <- ggplot(coords, aes(PC2, PC5)) +
  geom_point(data=coords, aes(x=PC2, y=PC5, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC5, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC5", PC5per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC5,label=Name),force=1, point.padding = 0.2) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic()

pc2pc6 <- ggplot(coords, aes(PC2, PC6)) +
  geom_point(data=coords, aes(x=PC2, y=PC6, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC6, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC6", PC6per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC6,label=Name),force=1, point.padding = 0.2) +
  scale_color_brewer(palette = "Set3") +
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic()

## yuk with field

pca.yuk_f<- prcomp(yuk_f.log, scale=FALSE)  # PCA with centering and scaling
pca.yuk_f$rotation  # The loadings are here


plot(pca.yuk_f, type = "l")
summary(pca.yuk_f)
#summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(pca.yuk_f$x)
sampleVals<-data.frame(pca.yuk_f$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(yuk_f),1,nchar(colnames(yuk_f))-2)
samples[c(8, 14, 18,23)] <- c("YD1w", "YD2", "YD2w", "YCC")
#samples <- c(samples, rep("Y-F2003", 2))
ecotype <- c(rep("Yukon", 23))
#cond <- substring(colnames(drought)[1:16],2,3)
#cond <- c(cond, rep("field", 2))
label <- c(rep("yes", 43))


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 23), Y=rep(0, 23), sampleVals,
                   Samples = samples, Name = colnames(yuk_f),
                   Condition = substr(samples, 2,nchar(samples)),
                   Ecotype = ecotype, stringsAsFactors=FALSE)

coords$Condition[c(8, 18, 19, 20, 23)] <- c("D1", "D2" , "WW1", "WW1", "YCC")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2", "F2003", "YCC"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (pca.yuk_f$sdev^2/sum(pca.yuk_f$sdev^2)) * 100
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")
#theme_set(theme_bw()) # make background white
library(ggplot2)
library(ggrepel)
pc23plot<-ggplot(exprVals, aes(PC2, PC3)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC3, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC3, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC3", PC3per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC3,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)

ggsave("yuk_fpc23plot.pdf", plot=pc23plot, width=9, height=7, units="in")
#viridis::scale_color_viridis(discrete = TRUE, option = "C") +
#scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +

## genes for Elizabeth's class
pos2 <- exprVals[which(exprVals[,2] > 0), c(2,3)]
neg2 <- exprVals[which(exprVals[,2] < 0), c(2,3)]

pos3 <- exprVals[which(exprVals[,3] > 0), c(2,3)]
#neg3 <- exprVals[which(exprVals[,3] < 0), c(2,3)]

pos2_50 <- pos2[order(pos2$PC2, decreasing=TRUE)[1:50],]
#pos2_50_fpkm <- merge(pos2_50, fpkm, by.x=0, by.y=0)
neg2_50 <- neg2[order(neg2$PC2, decreasing=FALSE)[1:50],]
#neg2_50_fpkm <- merge(neg2_50, fpkm, by.x=0, by.y=0)

pos3_50 <- pos3[order(pos3$PC3, decreasing=TRUE)[1:50],]


## pos2_50 (YWW1 cabinet)
pos2_fpkm <- merge(pos2_50, fpkm[,c(1:29,50:64)], by.x=0, by.y="gene_id")
neg2_fpkm <- merge(neg2_50, fpkm[,c(1:29,50:64)], by.x=0, by.y="gene_id")
pos3_fpkm <- merge(pos3_50, fpkm[,c(1:29,50:64)], by.x=0, by.y="gene_id")


write.table(pos2_fpkm, "pos2fpkm.txt", sep="\t", quote=F, col.names=NA)
write.table(neg2_fpkm, "neg2fpkm.txt", sep="\t", quote=F, col.names=NA)
write.table(pos3_fpkm, "pos3fpkm.txt", sep="\t", quote=F, col.names=NA)

### PCA of lncRNAs only...
lncrna <- fpkm %>%
  filter(y_lnc==1 | s_lnc==1)
lncrna <- unique(lncrna[,c(1,10:52)])

rownames(lncrna) <- lncrna$gene_id
lncrna <- lncrna[,-1]
colnames(lncrna)[c(1,2,9,10,21,22,29,30,41,42)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                      "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2",
                                                      "YF2003.1", "YF2003.2")
lncrna <- lncrna[,-c(8, 14, 18, 19, 20, 28, 35, 39:43)]

lncrna.log <- log2(lncrna+1)
lnc.pca <- prcomp(lncrna.log, scale=FALSE)  # PCA with centering and scaling
lnc.pca$rotation  # The loadings are here


plot(lnc.pca, type = "l")
summary(lnc.pca)
#summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(lnc.pca$x)
sampleVals<-data.frame(lnc.pca$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(lncrna.log),1,nchar(colnames(lncrna.log))-2)
ecotype <- c(rep("Yukon", 15), rep("Shandong", 16))



## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 31), Y=rep(0, 31), sampleVals,
                   Samples = samples, Name = colnames(lncrna.log),
                   Condition = substr(samples, 2,nchar(samples)),
                   Ecotype = ecotype, stringsAsFactors=FALSE)

#coords$Condition[c(8, 18, 19, 20, 23)] <- c("D1", "D2" , "WW1", "WW1", "YCC")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (lnc.pca$sdev^2/sum(lnc.pca$sdev^2)) * 100
PC1per <- paste0("(", round(PoV[1], 2), "%", ")")
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

pc12plot<-ggplot(exprVals, aes(PC1, PC2)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC1, y=PC2, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC1, y=PC2, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC1", PC1per)) +
  scale_y_continuous(name = paste("PC2", PC2per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC1,PC2,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)

pc32plot<-ggplot(exprVals, aes(PC3, PC2)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC2, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC2, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC2", PC2per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC2,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)


pc24plot<-ggplot(exprVals, aes(PC2, PC4)) +
  #geom_point(data = coords, aes(x=PC4, y=PC2,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=20)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
ggsave("~/Documents/jobapps/figeys/talk/figs/lnc_pc24plot.pdf", plot=pc24plot, width=11, height=6, units="in")

## PCA of everythign BUT lnc!

### PCA of lncRNAs only...
nolncrna <- fpkm %>%
  filter(y_lnc==0 | s_lnc==0)

nolncrna <- unique(nolncrna[,c(1,10:52)])

rownames(nolncrna) <- nolncrna$gene_id
nolncrna <- nolncrna[,-1]
colnames(nolncrna)[c(1,2,9,10,21,22,29,30,41,42)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                     "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2",
                                                     "YF2003.1", "YF2003.2")
nolncrna <- nolncrna[,-c(8, 14, 18, 19, 20, 28, 35, 39:43)]

nolncrna.log <- log2(nolncrna+1)
nolnc.pca <- prcomp(nolncrna.log, scale=FALSE)  # PCA with centering and scaling
nolnc.pca$rotation  # The loadings are here


plot(nolnc.pca, type = "l")
summary(nolnc.pca)
  #summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(nolnc.pca$x)
sampleVals<-data.frame(nolnc.pca$rotation)


dim(exprVals)
dim(sampleVals)

samples <- substring(colnames(nolncrna.log),1,nchar(colnames(nolncrna.log))-2)
ecotype <- c(rep("Yukon", 15), rep("Shandong", 16))



## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 31), Y=rep(0, 31), sampleVals,
                   Samples = samples, Name = colnames(nolncrna.log),
                   Condition = substr(samples, 2,nchar(samples)),
                   Ecotype = ecotype, stringsAsFactors=FALSE)

#coords$Condition[c(8, 18, 19, 20, 23)] <- c("D1", "D2" , "WW1", "WW1", "YCC")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))
#coords$Label <- c(rep("", 7), "YD1.4w", rep("", 4), "YD2.1a", "YD2.1b", rep("", 3), "YD2.5w", "YWWold.1", "YWWold.1",
#                  rep("", 7), "SD1.4w", rep("", 5), "SD2.2a", "SD2.2b", rep("", 3), "SWWold.1", "SWWold.2", "", "", "")

PoV <- (nolnc.pca$sdev^2/sum(nolnc.pca$sdev^2)) * 100
PC1per <- paste0("(", round(PoV[1], 2), "%", ")")
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

pc12plot<-ggplot(exprVals, aes(PC1, PC2)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC1, y=PC2, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC1, y=PC2, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC1", PC1per)) +
  scale_y_continuous(name = paste("PC2", PC2per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC1,PC2,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)

pc32plot<-ggplot(exprVals, aes(PC3, PC2)) +
  #geom_point(data = coords, aes(x=PC2, y=PC3,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC2, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC2, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC2", PC2per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC2,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=15)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)


pc24plot<-ggplot(exprVals, aes(PC2, PC4)) +
  #geom_point(data = coords, aes(x=PC4, y=PC2,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC2, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC2, y=PC4, color = Condition, shape = Ecotype),
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5,size=4) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  theme(text=element_text(size=20)) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
pc24plot
ggsave("~/Documents/jobapps/figeys/talk/figs/lnc_pc24plot.pdf", plot=pc24plot, width=11, height=6, units="in")

