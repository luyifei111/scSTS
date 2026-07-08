library(numbat)
library(Seurat)
library(tidyverse)
library(dplyr)
library(anndata)
library(openxlsx)
library(ggplot2)
library(glmnet)
library(reshape2)
library(ComplexHeatmap)
library(MetaNeighbor)
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
library(hdf5r)
library(sceasy)
library(scCustomize)

files <- list.files(
  path = "/cluster3/yflu/STS/numbat/result/run_numbat_out_rds",
  pattern = "_all_out\\.rds$",
  full.names = TRUE
)
samplenames <- substr(files,53,nchar(files)-12)
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
anno_sample_cluster_extended <- anno_sample_cluster_extended[samplenames,]

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5ad")
STS.pega.tumor <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.tumor@meta.data <- metadata$obs

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/RMS_pegasus_250305.h5ad")
STS.pega.tumor_RMS <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/RMS_pegasus_250305.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.tumor_RMS@meta.data <- metadata$obs

fits = readRDS("/cluster3/yflu/STS/development/regression_sim/lrFoetalClustersV4_development.RDS")

preds = readRDS("/cluster3/yflu/STS/development/log_regression_250305.rds")
pp = do.call(cbind,preds)
colnames(pp) = names(fits)
pp = (1+exp(-pp))**-1

pp_RMS <- pp
pp_RMS <- as.data.frame(pp_RMS)
rownames(pp_RMS) <- colnames(STS.pega.tumor)
pp_RMS <- pp_RMS[colnames(STS.pega.tumor_RMS),]
STS.pega.tumor_RMS@meta.data <- cbind(STS.pega.tumor_RMS@meta.data,pp_RMS)

Plot_Density_Custom(STS.pega.tumor_RMS, features = c("Myoblasts-Myocytes","Neural progenitors","Muscle progenitors","Satellite cells"),
                    viridis_palette = "inferno",limit = c(0,0.025))

RMS_external <- readRDS("/cluster3/yflu/STS/public_data/RMS_atlas_final_20240130.rds")
#saveRDS(preds,"log_regression_RMS_external_250312.rds")
preds = readRDS("/cluster3/yflu/STS/development/log_regression_RMS_external_250312.rds")
pp = do.call(cbind,preds)
colnames(pp) = names(fits)
pp = (1+exp(-pp))**-1
rownames(pp) <- colnames(RMS_external)
RMS_external@meta.data <- cbind(RMS_external@meta.data,pp)

library(ggplot2)
library(viridis)
library(rlang)
library(dplyr)
library(tidyr)

ggplot_multi_meta_continuous <- function(sts, meta_features, reduction = "umap", point_size = 1) {
  
  # 检查降维
  if(!reduction %in% names(sts@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in Seurat object."))
  }
  
  # 提取降维坐标
  coords <- Embeddings(sts[[reduction]])
  coords_df <- as.data.frame(coords)
  coords_df$cell <- rownames(coords_df)
  
  # 检查 metadata 列是否存在
  missing_cols <- meta_features[!meta_features %in% colnames(sts@meta.data)]
  if(length(missing_cols) > 0) stop(paste("Metadata column(s) not found:", paste(missing_cols, collapse = ", ")))
  
  # 检查是否连续
  non_numeric <- meta_features[!sapply(meta_features, function(x) is.numeric(sts@meta.data[[x]]))]
  if(length(non_numeric) > 0) stop(paste("Metadata column(s) not numeric:", paste(non_numeric, collapse = ", ")))
  
  # 将数据整合为长表
  long_df <- coords_df %>%
    dplyr::select(1:2) %>%
    bind_cols(sts@meta.data[, meta_features, drop = FALSE]) %>%
    tidyr::pivot_longer(cols = all_of(meta_features), names_to = "feature", values_to = "value")
  
  # 按输入顺序设置 factor
  long_df$feature <- factor(long_df$feature, levels = meta_features)
  
  # 绘图
  p <- ggplot(long_df, aes(x = long_df[[1]], y = long_df[[2]], color = value)) +
    geom_point(size = point_size) +
    scale_color_viridis(option = "viridis", direction = 1) +
    facet_wrap(~feature, scales = "free") +
    theme_classic() +
    labs(color = "Value")
  
  return(p)
}

meta_features <- c("Myoblasts-Myocytes","Neural progenitors","Muscle progenitors","Satellite cells")

ggplot_multi_meta_continuous(STS.pega.tumor_RMS, meta_features = meta_features, reduction = "umap", point_size = 0.4)
ggplot_multi_meta_continuous(RMS_external, meta_features = meta_features, reduction = "umap_rpca", point_size = 0.1)

STS.pega.tumor_RMS$cell <- rownames(STS.pega.tumor_RMS@meta.data)

anno_RMS <- subset(anno_sample_cluster_extended,Disease == 'RMS')
samplenames_RMS <- rownames(anno_RMS)
files_RMS <- paste("/cluster3/yflu/STS/numbat/result/run_numbat_out_rds/",rownames(anno_RMS),"_all_out.rds",sep = "")
pdf_RMS <- paste("/cluster3/yflu/STS/numbat/dev/fig/",rownames(anno_RMS),".pdf",sep = "")

i=1
numbat_obj <- readRDS(files_RMS[i])
clone <- numbat_obj$clone_post
clone <- subset(clone,compartment_opt == "tumor")
clone <- subset(clone,clone_opt != 1)
clone <- as.data.frame(clone)
rownames(clone) <- clone$cell
clone <- subset(clone,cell %in% STS.pega.tumor_RMS$cell)
clone <- clone[,c(1,2)]
clone$clone_el <- ifelse(
  clone$clone == min(clone$clone),
  "early clone",
  "late clone"
)

for (i in 2:length(samplenames_RMS)) {
  numbat_obj <- readRDS(files_RMS[i])
  clone_1 <- numbat_obj$clone_post
  clone_1 <- subset(clone_1,compartment_opt == "tumor")
  clone_1 <- subset(clone_1,clone_opt != 1)
  clone_1 <- as.data.frame(clone_1)
  rownames(clone_1) <- clone_1$cell
  clone_1 <- subset(clone_1,cell %in% STS.pega.tumor_RMS$cell)
  clone_1 <- clone_1[,c(1,2)]
  clone_1$clone_el <- ifelse(
    clone_1$clone_opt == min(clone_1$clone_opt),
    "early clone",
    "late clone"
  )
  print(paste(i,"/21"))
  clone <- rbind(clone,clone_1)
}

STS.pega.tumor_RMS_1 <- subset(STS.pega.tumor_RMS,cell %in% clone$cell)
STS.pega.tumor_RMS_1$clone_el <- clone[colnames(STS.pega.tumor_RMS_1),]$clone_el

DimPlot(STS.pega.tumor_RMS_1,group.by ="clone_el",raster = F)

avg_scores_list <- list()
for (i in 1:length(samplenames_RMS)) {
  numbat_obj <- readRDS(files_RMS[i])
  clone <- numbat_obj$clone_post
  clone <- subset(clone,compartment_opt == "tumor")
  clone <- subset(clone,clone_opt != 1)
  STS.pega.tumor_RMS_sub <- subset(STS.pega.tumor_RMS,cell %in% clone$cell)
  clone <- as.data.frame(clone)
  rownames(clone) <- clone$cell
  clone <- clone[STS.pega.tumor_RMS_sub$cell,]
  STS.pega.tumor_RMS_sub$clone <- clone$clone_opt
  STS.pega.tumor_RMS_sub$clone_el <- ifelse(
    STS.pega.tumor_RMS_sub$clone == min(STS.pega.tumor_RMS_sub$clone),
    "early clone",
    "late clone"
  )
  
  # 要计算平均值的列
  meta_cols <- c("Myoblasts-Myocytes","Neural progenitors","Muscle progenitors","Satellite cells")
  # 按 clone_el 分组计算平均值，排除为0的细胞
  meta_df <- STS.pega.tumor_RMS_sub@meta.data
  
  # 1️⃣ baseline（所有细胞中 clone_el 的比例）
  baseline <- meta_df %>%
    group_by(clone_el) %>%
    summarise(total_n = n(), .groups = "drop") %>%
    mutate(total_prop = total_n / sum(total_n))
  
  # 2️⃣ 转 long 格式
  meta_long <- meta_df %>%
    pivot_longer(cols = all_of(meta_cols),
                 names_to = "feature",
                 values_to = "value")
  
  # 3️⃣ 对每个 feature 单独筛选 >0.25，并计算 enrichment
  enrich_df <- meta_long %>%
    filter(value > 0.25) %>%
    group_by(feature, clone_el) %>%
    summarise(high_n = n(), .groups = "drop") %>%
    group_by(feature) %>%
    mutate(high_prop = high_n / sum(high_n)) %>%
    left_join(baseline, by = "clone_el") %>%
    mutate(enrichment = high_prop / total_prop)
  
  avg_scores_list[[i]] <- enrich_df
  names(avg_scores_list)[i] <- samplenames_RMS[i]
  
  p=DimPlot(STS.pega.tumor_RMS_sub,group.by ="clone_el") +
    DimPlot(STS.pega.tumor_RMS_sub,group.by ="clone") +
    ggplot_multi_meta_continuous(STS.pega.tumor_RMS_sub, meta_features = meta_features, reduction = "umap", point_size = 0.4)
  p
  ggsave(filename = pdf_RMS[i], plot = p, width = 15, height = 5)  # 宽高可调
  print(paste(i,"/21"))
}

# 阈值向量
thresholds <- seq(0.1, 0.5, by = 0.05)

# 存储结果
avg_scores_thresh <- list()

for (i in 1:length(samplenames_RMS)) {
  numbat_obj <- readRDS(files_RMS[i])
  clone <- numbat_obj$clone_post
  clone <- subset(clone, compartment_opt == "tumor")
  clone <- subset(clone, clone_opt != 1)
  STS.pega.tumor_RMS_sub <- subset(STS.pega.tumor_RMS, cell %in% clone$cell)
  clone <- as.data.frame(clone)
  rownames(clone) <- clone$cell
  clone <- clone[STS.pega.tumor_RMS_sub$cell,]
  STS.pega.tumor_RMS_sub$clone <- clone$clone_opt
  STS.pega.tumor_RMS_sub$clone_el <- ifelse(
    STS.pega.tumor_RMS_sub$clone == min(STS.pega.tumor_RMS_sub$clone),
    "early clone",
    "late clone"
  )
  
  meta_cols <- c("Myoblasts-Myocytes","Neural progenitors","Muscle progenitors","Satellite cells")
  meta_df <- STS.pega.tumor_RMS_sub@meta.data
  
  # baseline
  baseline <- meta_df %>%
    group_by(clone_el) %>%
    summarise(total_n = n(), .groups = "drop") %>%
    mutate(total_prop = total_n / sum(total_n))
  
  # 转长格式
  meta_long <- meta_df %>%
    pivot_longer(cols = all_of(meta_cols),
                 names_to = "feature",
                 values_to = "value")
  
  # 在每个阈值下计算 enrichment
  for (thresh in thresholds) {
    enrich_df <- meta_long %>%
      filter(value > thresh) %>%
      group_by(feature, clone_el) %>%
      summarise(high_n = n(), .groups = "drop") %>%
      group_by(feature) %>%
      mutate(high_prop = high_n / sum(high_n)) %>%
      left_join(baseline, by = "clone_el") %>%
      mutate(enrichment = high_prop / total_prop,
             threshold = thresh,
             sample = samplenames_RMS[i])
    
    avg_scores_thresh[[paste(samplenames_RMS[i], thresh, sep = "_")]] <- enrich_df
  }
  
  print(i)
}

# 合并所有样本和阈值结果
enrich_by_threshold <- bind_rows(avg_scores_thresh)

# 保留同时有 early 和 late 的样本
enrich_by_threshold <- enrich_by_threshold %>%
  group_by(sample, feature, threshold) %>%
  filter(n_distinct(clone_el) == 2) %>%
  ungroup() %>%
  mutate(log2_enrichment = log2(enrichment + 1e-6))

# Wilcoxon 检验
stat_res_thresh <- enrich_by_threshold %>%
  group_by(threshold, feature) %>%
  summarise(
    p_value = wilcox.test(enrichment ~ clone_el, paired = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

colnames(enrich_by_threshold)[8] <- "hold"
df <- subset(as.data.frame(enrich_by_threshold),hold > 0.15&hold < 0.16)
write.xlsx(df[,c(1,2,9,10)],"/cluster3/yflu/STS/figure_data/2I.xlsx")
