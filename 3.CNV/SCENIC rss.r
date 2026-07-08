# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

loom <- open_loom("/cluster3/yflu/STS/scenic/loom/development_reference_40.loom")

STS_group <- read.csv("/cluster3/yflu/STS/scenic/auc/all_auc_40_group.csv")
rownames(STS_group) <- STS_group$barcodekey
AUC_mtx <- fread("/cluster3/yflu/STS/scenic/auc/all_auc_40_transposed.csv",sep = ",",header = T)
rownames(AUC_mtx) <- AUC_mtx$V1
AUC_mtx <- AUC_mtx[,-1]
STS_group <- STS_group[colnames(AUC_mtx),]

AUC_mtx <- as.matrix(as.data.frame(AUC_mtx))
rownames(STS_group) <- STS_group$barcodekey
group_df <- STS_group[, "STS_meta$\"2p25_1_Amp\"", drop = FALSE]
colnames(group_df) <- "Group"
group_df$Group <- as.character(group_df)

rss <- calcRSS(AUC_mtx, STS_group[,"STS_meta$\"2p25_1_Amp\"", drop = FALSE])
rownames(rss) <- rownames(AUC_mtx)
rss <- fread("/cluster3/yflu/STS/scenic/auc/all_auc_40_rss.csv",sep = ",",header = T)
rss <- as.data.frame(rss)
rownames(rss) <- rss$V1
rss <- rss[,-1]