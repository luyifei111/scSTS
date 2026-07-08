library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(BiocParallel)
library(SeuratDisk)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(oncoPredict)
library(SeuratDisk)
library(Seurat)
library(rhdf5)
library(anndata)
library(stringr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(sctransform)
library(ggsci)
library(umap)
library(limma)
library(viridis)
library(patchwork)
library(SingleR)
library(xlsx)
library(harmony)
library(NMF)
library(irGSEA)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(scales)
library(DoubletFinder)
library(caret)
library(hdf5r)
library(ggpubr)
library(openxlsx)
library(pheatmap)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(MetaNeighbor)
library(ggplot2)
library(ggpubr)
library(Hmisc)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

metadata <- read_h5ad("/cluster3/yflu/STS/microenvironmnet/data_normal_TNK.h5ad")
STS_pega_TNK <- LoadH5Seurat("/cluster3/yflu/STS/microenvironmnet/data_normal_TNK.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega_TNK@meta.data <- metadata$obs
Idents(STS_pega_TNK) <- STS_pega_TNK@meta.data$louvain_labels_4
celltype <- read.xlsx("/cluster3/yflu/STS/microenvironmnet/T/celltypes res4.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS_pega_TNK)
STS_pega_TNK <- RenameIdents(STS_pega_TNK, new.cluster.ids)
DimPlot(STS_pega_TNK, reduction = "umap",label = TRUE,raster = F)

STS_pega_TNK@meta.data$Celltype_new <- Idents(STS_pega_TNK)
celltypes <- unique(levels(STS_pega_TNK@meta.data$Celltype_new))[-15]
STS_pega_TNK_1 <- subset(STS_pega_TNK,Celltype_new %in% celltypes)
STS_pega_TNK_1@meta.data$Celltype_new <- droplevels(STS_pega_TNK_1@meta.data$Celltype_new)
DimPlot(STS_pega_TNK_1, reduction = "umap",label = TRUE,raster = F)

cell_order <- levels(STS_pega_TNK_1@meta.data$Celltype_new)[c(1,8,9,4,7,2,11,5,14,10,12,3,6,13,15)]
STS_pega_TNK_1$Celltype_new <- factor(STS_pega_TNK_1$Celltype_new,levels = cell_order)

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(12)[c(1:4)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(12)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(12)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(11,12)])(8)[c(1)],
         "#ADD8E6",
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)])
DimPlot(
  STS_pega_TNK_1,
  reduction = "umap",
  label = TRUE,
  raster = FALSE,
  group.by = "Celltype_new",
  cols = cols,
)
FeaturePlot(STS_pega_TNK_1,features = c("FOXP3"))
STS_pega_TNK_1$Channel <- as.character(STS_pega_TNK_1$Channel)
STS_pega_TNK_2 <- subset(STS_pega_TNK_1,Channel %in% rownames(anno_sample_cluster_extended))

anno <- anno_sample_cluster_extended

# 取出要添加的新向量，名称为样本 ID
louvain_vec <- anno$Louvain
names(louvain_vec) <- rownames(anno)

# 匹配到 Seurat 对象的 Channel
matched_louvain <- louvain_vec[STS_pega_TNK_2$Channel]
names(matched_louvain) <- rownames(STS_pega_TNK_2@meta.data)

# 添加到 metadata
STS_pega_TNK_2$Louvain <- matched_louvain
matched_louvain <- unname(matched_louvain)
STS_pega_TNK_2$louvain <- matched_louvain 
lv <- STS_pega_TNK_2$louvain
cats <- sort(unique(lv))

# 为每个类别创建一列
for (k in cats) {
  new_col <- ifelse(lv == k,
                    paste0("Louvain_", k),
                    paste0("Not Louvain_", k))
  STS_pega_TNK_2[[paste0("louvain_", k)]] <- new_col
}

DimPlot(STS_pega_TNK_2,split.by = "louvain_2")  
STS_pega_TNK_2_sce <- as.SingleCellExperiment(STS_pega_TNK_2)
STS_pega_TNK_2_milo <- Milo(STS_pega_TNK_2_sce)

STS_pega_TNK_2_milo <- buildGraph(STS_pega_TNK_2_milo, k = 10, d = 15, reduced.dim = "PCA_REGRESSED_HARMONY")
STS_pega_TNK_2_milo <- makeNhoods(STS_pega_TNK_2_milo, prop = 0.05, k = 10, d=15, refined = TRUE, reduced_dims = "PCA_REGRESSED_HARMONY")
plotNhoodSizeHist(STS_pega_TNK_2_milo)
STS_pega_TNK_2_milo <- countCells(STS_pega_TNK_2_milo, meta.data = as.data.frame(colData(STS_pega_TNK_2_milo)), sample="Channel")
Counts <- as.data.frame(nhoodCounts(STS_pega_TNK_2_milo))

TNK_design <- data.frame(colData(STS_pega_TNK_2_milo))[,c("Channel", "louvain","louvain_1","louvain_2","louvain_3","louvain_4","louvain_5")]
TNK_design$Channel <- as.factor(TNK_design$Channel) 

TNK_design <- distinct(TNK_design)
rownames(TNK_design) <- TNK_design$Channel

options(warn = -1)   # 不要把 warning 当成 error

# rd <- reducedDim(STS_pega_TNK_2_milo, "PCA_REGRESSED_HARMONY")
# rd_new <- as(rd, "CsparseMatrix")  # 推荐格式（不会产生 warning）
# 
# reducedDim(STS_pega_TNK_2_milo, "PCA_REGRESSED_HARMONY") <- rd_new

patched_calc_distance <- function(in.x) {
  dist.list <- lapply(seq_len(nrow(in.x)), FUN = function(i) {
    i.dist <- apply(in.x, 1, FUN = function(P) sqrt(sum((P - in.x[i, ])^2)))
    list(
      rowIndex = rep(i, nrow(in.x)),
      colIndex = seq_len(length(i.dist)),
      dist = i.dist
    )
  })
  
  dist.df <- do.call(rbind.data.frame, dist.list)
  
  # ★ 修复点：repr="C" 而不是 "T"
  out.dist <- sparseMatrix(
    i = dist.df$rowIndex,
    j = dist.df$colIndex,
    x = dist.df$dist,
    dimnames = list(rownames(in.x), rownames(in.x)),
    repr = "C"
  )
  
  return(out.dist)
}

# 注入 miloR 命名空间，覆盖原函数
assignInNamespace(".calc_distance", patched_calc_distance, ns = "miloR")

STS_pega_TNK_2_milo <- calcNhoodDistance(STS_pega_TNK_2_milo, d=15, reduced.dim = "PCA_REGRESSED_HARMONY")
saveRDS(STS_pega_TNK_2_milo,"/cluster3/yflu/STS/STS_pega_TNK_2_milo.rds")

da_results <- testNhoods(STS_pega_TNK_2_milo, design = ~ louvain_1, design.df = TNK_design)
saveRDS(da_results,"/cluster3/yflu/STS/cpdb/da_results_1.rds")
da_results <- testNhoods(STS_pega_TNK_2_milo, design = ~ louvain_2, design.df = TNK_design)
saveRDS(da_results,"/cluster3/yflu/STS/cpdb/da_results_2.rds")
da_results <- testNhoods(STS_pega_TNK_2_milo, design = ~ louvain_3, design.df = TNK_design)
saveRDS(da_results,"/cluster3/yflu/STS/cpdb/da_results_3.rds")
da_results <- testNhoods(STS_pega_TNK_2_milo, design = ~ louvain_4, design.df = TNK_design)
saveRDS(da_results,"/cluster3/yflu/STS/cpdb/da_results_4.rds")
da_results <- testNhoods(STS_pega_TNK_2_milo, design = ~ louvain_5, design.df = TNK_design)
saveRDS(da_results,"/cluster3/yflu/STS/cpdb/da_results_5.rds")

STS_pega_TNK_2_milo <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/STS_pega_TNK_2_milo.rds")

da_results_1 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_1.rds")
da_results_1 <- annotateNhoods(STS_pega_TNK_2_milo, da_results_1, coldata_col = "Celltype_new")
da_results_1$Celltype_new <- factor(da_results_1$Celltype_new, levels = rev(unique(da_results_1$Celltype_new)[order(unique(da_results_1$Celltype_new))]))
p1 <- plotDAbeeswarm(da_results_1, group.by = "Celltype_new")

da_results_2 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_2.rds")
da_results_2 <- annotateNhoods(STS_pega_TNK_2_milo, da_results_2, coldata_col = "Celltype_new")
da_results_2$Celltype_new <- factor(da_results_2$Celltype_new, levels = rev(unique(da_results_2$Celltype_new)[order(unique(da_results_2$Celltype_new))]))
p2 <- plotDAbeeswarm(da_results_2, group.by = "Celltype_new")

da_results_3 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_3.rds")
da_results_3 <- annotateNhoods(STS_pega_TNK_2_milo, da_results_3, coldata_col = "Celltype_new")
da_results_3$Celltype_new <- factor(da_results_3$Celltype_new, levels = rev(unique(da_results_3$Celltype_new)[order(unique(da_results_3$Celltype_new))]))
p3 <- plotDAbeeswarm(da_results_3, group.by = "Celltype_new")

da_results_4 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_4.rds")
da_results_4 <- annotateNhoods(STS_pega_TNK_2_milo, da_results_4, coldata_col = "Celltype_new")
da_results_4$Celltype_new <- factor(da_results_4$Celltype_new, levels = rev(unique(da_results_4$Celltype_new)[order(unique(da_results_4$Celltype_new))]))
p4 <- plotDAbeeswarm(da_results_4, group.by = "Celltype_new")

da_results_5 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_5.rds")
da_results_5 <- annotateNhoods(STS_pega_TNK_2_milo, da_results_5, coldata_col = "Celltype_new")
da_results_5$Celltype_new <- factor(da_results_5$Celltype_new, levels = rev(unique(da_results_5$Celltype_new)[order(unique(da_results_5$Celltype_new))]))
p5 <- plotDAbeeswarm(da_results_5, group.by = "Celltype_new")

(p1 | p2 | p3 | p4 | p5)

STS_pega_M_2_milo <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/STS_pega_M_2_milo.rds")

da_results_M_1 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_M_1.rds")
da_results_M_1 <- annotateNhoods(STS_pega_M_2_milo, da_results_M_1, coldata_col = "Celltype_new")
da_results_M_1$Celltype_new <- factor(da_results_M_1$Celltype_new, levels = rev(unique(da_results_M_1$Celltype_new)[order(unique(da_results_M_1$Celltype_new))]))
p1 <- plotDAbeeswarm(da_results_M_1, group.by = "Celltype_new")

da_results_M_2 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_M_2.rds")
da_results_M_2 <- annotateNhoods(STS_pega_M_2_milo, da_results_M_2, coldata_col = "Celltype_new")
da_results_M_2$Celltype_new <- factor(da_results_M_2$Celltype_new, levels = rev(unique(da_results_M_2$Celltype_new)[order(unique(da_results_M_2$Celltype_new))]))
p2 <- plotDAbeeswarm(da_results_M_2, group.by = "Celltype_new")

da_results_M_3 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_M_3.rds")
da_results_M_3 <- annotateNhoods(STS_pega_M_2_milo, da_results_M_3, coldata_col = "Celltype_new")
da_results_M_3$Celltype_new <- factor(da_results_M_3$Celltype_new, levels = rev(unique(da_results_M_3$Celltype_new)[order(unique(da_results_M_3$Celltype_new))]))
p3 <- plotDAbeeswarm(da_results_M_3, group.by = "Celltype_new")

da_results_M_4 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_M_4.rds")
da_results_M_4 <- annotateNhoods(STS_pega_M_2_milo, da_results_M_4, coldata_col = "Celltype_new")
da_results_M_4$Celltype_new <- factor(da_results_M_4$Celltype_new, levels = rev(unique(da_results_M_4$Celltype_new)[order(unique(da_results_M_4$Celltype_new))]))
p4 <- plotDAbeeswarm(da_results_M_4, group.by = "Celltype_new")

da_results_M_5 <- readRDS("/cluster3/yflu/STS/cpdb/milo_30/da_results_M_5.rds")
da_results_M_5 <- annotateNhoods(STS_pega_M_2_milo, da_results_M_5, coldata_col = "Celltype_new")
da_results_M_5$Celltype_new <- factor(da_results_M_5$Celltype_new, levels = rev(unique(da_results_M_5$Celltype_new)[order(unique(da_results_M_5$Celltype_new))]))
p5 <- plotDAbeeswarm(da_results_M_5, group.by = "Celltype_new")

(p1 | p2 | p3 | p4 | p5)