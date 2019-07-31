library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
setwd("~/R/Eutrema/PS")



trimCounts <- read.table("counts")

samples <- c("ps10_S1",
            "ps11_S2",
            "ps12_S3",
            "ps37_S4",
            "ps38_S5",
            "ps40_S6",
            "ps41_S7",
            "ps44_S8",
            "ps46_S9",
            "ps48_S10",
            "ps49_S11",
            "ps50_S12")
trimCounts <- cbind(trimCounts, trimCounts[,3]/trimCounts[,2],
                    rep(c(1, 2), nrow(trimCounts)/2),
                    rep(c(1, 1, 2, 2), nrow(trimCounts)/4),
                    rep(samples, each = 4))
colnames(trimCounts) <- c("File", "Raw", "Trimmed", "SurvivedRatio", "Read", "Lane", "Sample")

tCmelt <- melt(trimCounts[,c("File", "Raw", "Trimmed")], id = "File")


Trim <- ggplot(tCmelt, aes(x = File, y = value/1000000, fill = factor(variable))) +
  geom_bar(stat = "identity")+
  labs(y = "Millions of reads") +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())
ggsave("TrimStack.png", Trim, height = 6, width = 8, dpi = 125)
trimCounts$Sample <- as.factor(trimCounts$Sample)

trimRatio <- ggplot(trimCounts, aes(x = Sample, y = SurvivedRatio)) +
  geom_boxplot(notch = F)+
  ylim(0.75, 0.9) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

ggsave("TrimRatio.png", trimRatio, height = 5, width = 6, dpi = 125)

#tSmelt <- melt(trimCounts[,c("Sample", "Raw", "Trimmed")], id = "Sample")
#
# ggplot(tSmelt, aes(x = Sample, y = value/1000000, fill = factor(variable))) +
#   geom_bar(stat = "identity", position = "dodge")+
#   labs(y = "Average of millions of reads") +
#   stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
#   theme(text = element_text(size=10,
#                             family="serif"),
#         axis.text.x = element_text(size=10,angle=75, hjust=1),
#         legend.position = "top",
#         legend.title = element_blank())

trimAve <- aggregate(. ~ Sample, trimCounts[,c(2, 3, 7)], mean)
trimAve
trimSD <- aggregate(. ~ Sample, trimCounts[,c(2, 3, 7)], sd)
trimSD
tAvem <- melt(trimAve, id="Sample")
tAvem
tSDm <- melt(trimSD, id = "Sample")
tAvem <- cbind(tAvem, tSDm[,3])
colnames(tAvem) <- c("Sample", "Treat", "Count", "SD")
tAvem <- transform(tAvem, Count = Count/1000000, SD = SD/1000000)


#to make the error bars center at the middle of the bar
#you have to determine the width of the bar, then the
#dodge position of the error bar would be the same as that
#of the bar and the width of the error bar would be half
#the value of its dodge position
tAve <- ggplot(tAvem, aes(x = Sample, y = Count, fill = Treat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  geom_errorbar(aes(ymin=Count - SD, ymax = Count + SD), size = 0.5, position = position_dodge(0.9), width = 0.45) +
  ylab("Millions of reads") + xlab("Sample (n = 4)") +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())
tAve

ggsave("TrimAve.png", tAve, height = 5, width = 6, dpi = 125)

merge_data <- read.table("~/scratch/PS/STAR/tot", header = T)

merge_data <- cbind(samples, merge_data)

mergMelt <- melt(merge_data, id = "samples")
mergMelt

merge_plot <- ggplot(mergMelt, aes(x = samples, y = value/1000000, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Samples", y = "Millions of mapped reads") +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

print(merge_plot)



######MAPQ######
library(dplyr)

unMerge <- read.table("unMerge_count")

Merge <- read.table("Merged_count")

counts <- cbind(rep(samples, each = 4), unMerge, Merge)

colnames(counts) <- c("Sample", "Unmerged", "Score", "Merged")

counts <- counts[,1:ncol(counts)-1]

UnMerged <- counts[,c(1, 3, 2)]
unmergeTot <- rep(aggregate(UnMerged$Unmerged, by = list(Category = UnMerged$Sample), FUN = sum)[,2], each = 4)
UnMerged <- cbind(UnMerged, unmergeTot)
UnMerged <- cbind(UnMerged, UnMerged$Unmerged/UnMerged$unmergeTot)
colnames(UnMerged)[ncol(UnMerged)] <- "Ratio"
UnMerged
uniq <- UnMerged %>% filter(UnMerged$Score == 255)
trim.sum <- tapply(trimCounts$Trimmed, trimCounts$Sample, FUN = sum)
uniq <- cbind(uniq, trim.sum)
uniq$TrimRation <- uniq$unmergeTot/uniq$trim.sum*100
uniq$unmergeTot <- as.double(uniq$unmergeTot)
formatted <- data.frame(format(uniq[, 4], digits = 3, scientific = T), format(uniq[,  7], digits = 3, scientific = F), stringsAsFactors = F)
rownames(formatted) <- rownames(uniq)
colnames(formatted) <- c("Counts", "mapRatio")
write.table(formatted, file = "mapRatio.tab", quote = F, row.names = , col.names = F, sep = " & ", eol = " \\\\\n")

M_erged <- counts[,c(1, 3, 4)]
mergeTot <- rep(aggregate(M_erged$Merged, by = list(Category = M_erged$Sample), FUN = sum)[,2], each = 4)
M_erged <- cbind(M_erged, mergeTot)
M_erged <- cbind(M_erged, M_erged$Merged/M_erged$mergeTot)
colnames(M_erged)[ncol(M_erged)] <- "Ratio"
M_erged

UnMerge_plot <- ggplot(UnMerged, aes(x = Sample, y = Unmerged/unmergeTot, fill = factor(Score))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

ggsave("UnMerged.png", UnMerge_plot, height = 4, width = 5, dpi = 125)

Merge_plot <- ggplot(M_erged, aes(x = Sample, y = Merged, fill = factor(Score))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10,angle=75, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

ggsave("Merged.png", Merge_plot, height = 4, width = 5, dpi = 125)

unMergeUniq <- UnMerged %>% filter( Score == 255)
unMergeAOV <- aov(unMergeUniq$Ratio ~ unMergeUniq$Sample)

MergeUniq <- M_erged %>% filter(Score == 255)
MergeAOV <- aov(MergeUniq$Ratio ~ MergeUniq$Sample)


res <- t.test(unMergeUniq$Ratio, MergeUniq$Ratio, val.equal = T)
res

allUniq <- cbind(samples, unMergeUniq$Ratio, MergeUniq$Ratio)
colnames(allUniq) <- c("Sample", "Unmerged", "Merged")

#find outliers oooo
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

allUniq <- melt(as.data.frame(allUniq), id = "Sample")


#get the outliers and put the sample name in the outliers column and other ones as NA
mergeBox <- allUniq %>%
  group_by(variable) %>% mutate( outlier = ifelse(is_outlier(as.numeric(value)), samples[Sample], as.numeric(NA))) %>%
  ggplot(aes(x = as.factor(variable), y = as.numeric(value), fill =as.factor(variable))) +
  geom_boxplot() +
  labs(x = "Unmerged/Merged reads (n = 12)", y = "Ratio of uniquely mapped reads") +
  geom_text(aes(label = outlier), na.rm = TRUE, vjust = -1, size=3, family="serif") +
  ylim(0.985, 0.999) +
  theme(text = element_text(size=10,
                            family="serif"),
        axis.text.x = element_text(size=10, hjust=1),
        legend.position = "top",
        legend.title = element_blank())

ggsave("mergeBox.png", mergeBox, height = 4, width = 3, dpi = 125)

