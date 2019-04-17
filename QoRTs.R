library(QoRTs)
library(DESeq2)

setwd("~/scratch/PS/QoRTs/smolMerged")
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
completeAndCheckDecoder(samples)
incompleteDecoder <- data.frame(unique.ID = samples, group.ID = rep(c("LPLS", "HPLS", "LPHS", "HPHS"), each = 3))

decoder.data <- completeAndCheckDecoder(incompleteDecoder)

res <- read.qc.results.data("./", decoder = decoder.data, calc.DESeq2 = T, calc.edgeR = T)

makeMultiPlot.all(res, outfile.dir = "./", plot.device.name = "pdf")
