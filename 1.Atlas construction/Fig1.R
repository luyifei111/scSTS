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

plot_cellnumber_compare <- function(df, group_col = "Disease", mal_col = "Malignancy",
                                    start_col = 1, end_col = NULL,
                                    hide_ns = TRUE, step = 0.05,
                                    cols = NULL, dodge_width = 0.8, box_width = 0.6) {
  if (is.null(end_col)) {
    end_col <- ncol(df) - 2  # 默认最后两列是 Disease / Malignancy
  }
  plot_list <- list()
  
  # 如果提供了颜色，保证 group_col 是 factor 且顺序一致
  if (!is.null(cols)) {
    df[[group_col]] <- factor(df[[group_col]], levels = names(cols))
  }
  
  for (i in start_col:end_col) {
    tmp <- df[, c(i, which(colnames(df) %in% c(group_col, mal_col)))]
    colnames(tmp)[1] <- "proportion"
    tmp[[mal_col]] <- as.character(tmp[[mal_col]])
    
    diseasenames <- sort(unique(tmp[[mal_col]]))
    my_comparisons <- combn(diseasenames, 2, simplify = FALSE)
    
    my_comparisons_sig <- list()
    for (j in 1:length(my_comparisons)) {
      por1 <- tmp[tmp[[mal_col]] == my_comparisons[[j]][1], "proportion"]
      por2 <- tmp[tmp[[mal_col]] == my_comparisons[[j]][2], "proportion"]
      if (length(por1) > 1 & length(por2) > 1) {
        test <- t.test(por1, por2)
        if (is.na(test$p.value)) test$p.value <- 1
        if (test$p.value < 0.05) {
          my_comparisons_sig <- append(my_comparisons_sig, list(my_comparisons[[j]]))
        }
      }
    }
    
    # 关键修改：interaction，确保 group 宽度一致
    P <- ggplot(tmp, aes(x = !!sym(mal_col), y = proportion, fill = !!sym(group_col))) +
      stat_boxplot(geom = "errorbar", width = 0.15, position = position_dodge(width = dodge_width)) +
      geom_boxplot(position = position_dodge(width = dodge_width), width = box_width, outlier.color = "white") +
      geom_jitter(shape = 21, size = 1.5, position = position_dodge(width = dodge_width)) +
      ggtitle(colnames(df)[i]) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )
    
    # 应用自定义配色
    if (!is.null(cols)) {
      P <- P + scale_fill_manual(values = cols, drop = FALSE)
    }
    
    # 显著性标注
    if (length(my_comparisons_sig) > 0) {
      ymax <- max(tmp$proportion, na.rm = TRUE)
      label_y <- ymax + (seq_along(my_comparisons_sig)) * (ymax * step)
      
      P <- P + stat_compare_means(
        comparisons = my_comparisons_sig,
        method = "t.test",
        label = "p.signif",
        hide.ns = hide_ns,
        label.y = label_y
      )
    }
    
    plot_list[[colnames(df)[i]]] <- P
  }
  return(plot_list)
}

anno_merged <- read.xlsx("/cluster/home/yflu/STS/pegasus/anno_merged.xlsx","Sheet2")

aurocs_Disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")

Disease_order <- as.data.frame(table(anno_merged$Disease))
Disease_order <- Disease_order[order(Disease_order$Freq,decreasing = T),]
Disease_order <- as.character(Disease_order$Var1)

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')

Diseasenames <- rownames(aurocs_Disease)[p$tree_row$order]
Diseasenames <- substr(Diseasenames,5,nchar(Diseasenames))

cols <- c(
  "Lipoblastoma"        = "#8DD3C7",
  "Synovial sarcoma"    = "#97D7C5",
  "Liposarcoma"         = "#A1DBC3",
  "Spindle cell tumor"  = "#ACDFC1",
  "Infantile fibrosarcoma" = "#B6E3BF",
  "Aggressive fibromatosis" = "#C0E7BD",
  "IMT"                 = "#CBEBBC",
  "ASPS"                = "#D5EFBA",
  "NF"                  = "#FFFFB3",
  "Schwannoma"          = "#F5F5B8",
  "Lymphangioma"        = "#FCCDE5",
  "Hemangioma"          = "#F0D1E1",
  "KHE"                 = "#E4D5DD",
  "MPNST"               = "#B3DE69",
  "Undifferentiated sarcoma" = "#80B1D3",
  "Angiosarcoma"        = "#FDB462",
  "RMS"                 = "#ADD8E6",
  "Pecoma"              = "#BC80BD",
  "EWS"            = "#BE8FBE",
  "MRT"                 = "#C09EBF",
  "Peritumor" = "#F18EAC"
)

cols <- c(
  "LPB"        = "#8DD3C7",
  "SS"    = "#97D7C5",
  "LPS"         = "#A1DBC3",
  "SCT"  = "#ACDFC1",
  "IFS" = "#B6E3BF",
  "AF" = "#C0E7BD",
  "IMT"                 = "#CBEBBC",
  "ASPS"                = "#D5EFBA",
  "NF"                  = "#FFFFB3",
  "SWN"          = "#F5F5B8",
  "LYM"        = "#FCCDE5",
  "HE"          = "#F0D1E1",
  "KHE"                 = "#E4D5DD",
  "MPNST"               = "#B3DE69",
  "US" = "#80B1D3",
  "AS"        = "#FDB462",
  "RMS"                 = "#ADD8E6",
  "PECOMA"              = "#BC80BD",
  "EWS"            = "#BE8FBE",
  "MRT"                 = "#C09EBF",
  "Peritumor" = "#F18EAC"
)

anno_merged$Disease <- factor(anno_merged$Disease,names(cols))
anno_merged$Malignancy <- factor(anno_merged$Malignancy,c("Benign","Malignant","Peritumor"))

anno_sorted <- anno_merged %>%
  arrange(Malignancy, Disease, desc(scRNA), desc(Xenium), desc(WES))

samples <- anno_sorted$Sample

# 构造一个全 0 的矩阵，10 个基因 × n 个样本
expr_mat <- matrix(0, nrow = 10, ncol = length(samples))
rownames(expr_mat) <- paste0("Gene", 1:10)
colnames(expr_mat) <- samples

# 提取第 2-6 列作为注释
anno_col <- anno_sorted[, 4:8]
rownames(anno_col) <- anno_sorted$Sample  

pheatmap(expr_mat,
         annotation_col = anno_col,
         annotation_colors = list(Disease = cols),
         breaks = c(-1, 0, 1))  # 避免全零报错

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

tumor_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_tumor_5_6.xlsx","Sheet 1")
normal_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_normal_5_6.xlsx","Sheet 1")
barcodes <- rbind(tumor_barcodes,normal_barcodes)
STS_pega@meta.data$barcode <- rownames(STS_pega@meta.data)
STS_pega_1 <- subset(STS_pega,barcode %in% barcodes$barcode)
rownames(barcodes) <- barcodes$barcode
barcodes <- barcodes[colnames(STS_pega_1),]
STS_pega_1$group <- barcodes$group
table <- as.data.frame.array(table(STS_pega_1$louvain_labels,STS_pega_1$group))

data_normal_Group_celltype <- read.csv("/cluster/home/yflu/STS/pegasus/data_normal_Group_celltype.csv")
barcodes_celltype <- tumor_barcodes
barcodes_celltype$barcodekey <- barcodes_celltype$barcode
barcodes_celltype$Group <- barcodes_celltype$group
barcodes_celltype$Celltype <- barcodes_celltype$group
barcodes_celltype <- rbind(barcodes_celltype[,c(4:6)],data_normal_Group_celltype[,c(1:3)])
rownames(barcodes_celltype) <- barcodes_celltype$barcodekey
barcodes_celltype <- barcodes_celltype[colnames(STS_pega_1),]
STS_pega_1$celltype <- barcodes_celltype$Celltype

DimPlot(STS_pega_1,group.by = "celltype",raster = T,reduction = "tsne")
STS_pega_1$celltype <- factor(STS_pega_1$celltype,levels = unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)])

tumor_sample <- anno_merged
meta <- STS_pega_1@meta.data

# 建立 Sample -> Disease 的映射表
sample2Disease <- setNames(as.character(tumor_sample$Disease), as.character(tumor_sample$Sample))

# 按照 Channel 填充 Disease
meta$Disease <- sample2Disease[as.character(meta$Channel)]
sample2Malignancy <- setNames(as.character(tumor_sample$Malignancy), as.character(tumor_sample$Sample))

meta$Malignancy <- sample2Malignancy[as.character(meta$Channel)]

# 更新回 Seurat 对象
STS_pega_1@meta.data <- meta

meta <- STS_pega_1@meta.data

#1D
# 统计每个 Disease - celltype 的数量和比例
cell_counts <- meta %>%
  group_by(Disease, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Disease) %>%
  mutate(prop = n / sum(n))

Disease_samples <- meta %>%
  select(Disease, Channel) %>%
  distinct() %>%
  group_by(Disease) %>%
  summarise(sample_n = n(), .groups = "drop")

# 2. 合并到 cell_counts
cell_counts <- cell_counts %>%
  left_join(Disease_samples, by = "Disease")
cell_counts$Disease <- factor(cell_counts$Disease,levels = levels(tumor_sample$Disease))

ggplot(cell_counts, aes(x = Disease, y = celltype)) +
  geom_point(aes(size = prop, fill = n), 
             shape = 21, color = "gray50", stroke = 0.5) +  # 黑色外框，0.3线宽
  scale_size_continuous(name = "Cell proprotion", range = c(1, 10)) +
  scale_fill_viridis(option = "inferno", direction = -1, name = "Cell counts") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(limits = rev)

metadata <- read_h5ad("/cluster3/yflu/STS/TMA/adata_TMA_merged_new.h5ad")
TMA_merged <- LoadH5Seurat("/cluster3/yflu/STS/TMA/adata_TMA_merged_new.h5seurat",meta.data = FALSE, misc = FALSE)
TMA_merged@meta.data <- metadata$obs

map_vec <- c(
  "Monocytes/macrophages" = "Mono/macro",
  "Vascular endothelial cells" = "Endothelial",
  "NK cells" = "T/NK cells",
  "Plasmacytoid dendritic cells" = "DC",
  "Lymphatic endothelial cells" = "Endothelial",
  "T cells" = "T/NK cells",
  "Myocytes" = "Other cells",
  "Schwann cell" = "Other cells",
  "Myoblasts" = "Other cells",
  "Schwann cells" = "Other cells",
  "Proliferating cells" =  "Other cells",
  "Tumor" =  "Tumor cells"
)

# 一次性替换
STS_pega_1@meta.data$Celltype_united <- 
  plyr::mapvalues(STS_pega_1@meta.data$celltype,
                  from = names(map_vec),
                  to   = unname(map_vec),
                  warn_missing = FALSE)
STS_pega_1@meta.data$Celltype_united <- as.character(STS_pega_1@meta.data$Celltype_united)

STS_pega_1.sce <- as.SingleCellExperiment(STS_pega_1)
TMA_merged.sce <- as.SingleCellExperiment(TMA_merged)

samplelist <- list(STS_pega_1.sce,TMA_merged.sce)
names(samplelist) <- c("scRNA-seq","Xenium")

fused_data = mergeSCE(samplelist)

global_hvgs = variableGenes(dat = fused_data, exp_labels = fused_data$study_id)

aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = fused_data,
                        study_id = fused_data$study_id,
                        cell_type = fused_data$celltype_new,
                        fast_version = TRUE)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(ggnewscale)
library(patchwork)

plot_alluvial_with_vertical_bar <- function(meta, cols, low_color="lightyellow", high_color="red") {
  
  if(!all(c("Disease", "Channel") %in% colnames(meta))) stop("meta must contain columns: 'Disease' and 'Channel'")
  
  # ---- 计算比例 ----
  sample_df <- meta %>%
    group_by(Disease) %>%
    summarise(n_sample = n_distinct(Channel), .groups="drop") %>%
    mutate(sample_prop = n_sample / sum(n_sample))
  
  cell_df <- meta %>%
    group_by(Disease) %>%
    summarise(n_cell = n(), .groups="drop") %>%
    mutate(cell_prop = n_cell / sum(n_cell))
  
  df <- left_join(sample_df, cell_df, by="Disease")
  
  df_long <- df %>%
    pivot_longer(cols=c(sample_prop, cell_prop),
                 names_to="Type", values_to="Prop") %>%
    mutate(Type = factor(Type,
                         levels=c("sample_prop","cell_prop"),
                         labels=c("Sample proportion","Cell proportion")))
  
  df_long$Disease <- factor(df_long$Disease, levels = names(cols))
  
  # ---- 平均每样本细胞数 ----
  avg_cells <- meta %>%
    group_by(Disease, Channel) %>%
    summarise(n_cells=n(), .groups="drop") %>%
    group_by(Disease) %>%
    summarise(avg_cells=mean(n_cells), .groups="drop") %>%
    mutate(relative = avg_cells / mean(avg_cells))
  avg_cells$Disease <- factor(avg_cells$Disease, levels = names(cols))
  
  # ---- 绘制 alluvial ----
  p_alluvial <- ggplot(df_long,
                       aes(x=Type, stratum=Disease, alluvium=Disease,
                           y=Prop, fill=Disease)) +
    geom_flow(stat="alluvium", lode.guidance="frontback",
              color="darkgray", alpha=0.7) +
    geom_stratum(alpha=0.9, color="black") +
    scale_fill_manual(values=cols) +
    theme_minimal() +
    ylab("Proportion") +
    xlab("") +
    ggtitle("Sample vs Cell proportions")
  
  # ---- 提取右侧 stratum 高度 ----
  stratum_data <- ggplot_build(p_alluvial)$data[[2]]  # geom_stratum
  right_stratum <- stratum_data %>%
    filter(x==2) %>%
    select(stratum, ymin, ymax) %>%
    distinct() %>%
    left_join(avg_cells, by=c("stratum"="Disease")) %>%
    mutate(Disease=stratum)
  
  # ---- 绘制右侧竖直条带，颜色范围固定0-2 ----
  p_bar <- ggplot(right_stratum) +
    geom_rect(aes(xmin=0, xmax=1, ymin=ymin, ymax=ymax, fill=relative), color=NA) +
    scale_fill_gradient(low=low_color, high=high_color, name="Relative\nAvg Cells", limits=c(0,3)) +
    theme_void() +
    theme(legend.position="right")
  
  # ---- 拼图，左右排列 ----
  p_final <- p_alluvial | p_bar + plot_layout(widths=c(4,1))
  
  return(p_final)
}

meta_xe <- TMA_merged@meta.data
meta_xe$Channel <- meta_xe$Sample

plot_alluvial_with_vertical_bar(meta,cols)
P2 <- plot_alluvial_with_vertical_bar(meta_xe,cols)

plot_celltype_compare_stack <- function(meta, meta_xe, cols) {
  library(dplyr)
  library(ggplot2)
  
  # 计算 meta 的 Celltype_united 比例
  df1 <- meta %>%
    group_by(Celltype_united) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prop = n / sum(n), dataset = "meta")
  
  # 计算 meta_xe 的 Celltype_united 比例
  df2 <- meta_xe %>%
    group_by(Celltype_united) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prop = n / sum(n), dataset = "meta_xe")
  
  df <- bind_rows(df1, df2)
  
  # 保持 Celltype_united 的一致顺序
  df$Celltype_united <- factor(df$Celltype_united, levels = unique(df$Celltype_united))
  
  # 绘制并列 stack bar
  p <- ggplot(df, aes(x = dataset, y = prop, fill = Celltype_united)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  # 分隔线
    scale_fill_manual(values = cols) +
    theme_minimal() +
    ylab("Proportion") +
    xlab("") +
    ggtitle("Celltype proportion in meta vs meta_xe") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_blank()
    )
  
  return(p)
}
meta <- STS_pega_1@meta.data
meta <- subset(meta,Disease != "Peritumor")
meta$Celltype_united <- factor(meta$Celltype_united,levels = unique(meta$Celltype_united)[c(8,3,12,1,6,9,2,7,4,5,11,10)])
meta_xe$Celltype_united <- factor(meta_xe$Celltype_united,levels = unique(meta$Celltype_united)[c(8,3,12,1,6,9,2,7,4,5,11,10)])

plot_celltype_compare_stack(meta, meta_xe, cols=cols_normal)

aurocs_SC_XE_250924 <- readRDS("/cluster3/yflu/STS/development/aurocs_SC_XE_250924.rds")
aurocs_SC_XE <- aurocs_SC_XE_250924[c(1:12),c(13:24)]
rownames(aurocs_SC_XE)[12] <- "scRNA-seq|Tumor cells"
pheatmap(aurocs_SC_XE[c(12,4,10,8,3,1,9,11,6,2,5,7),c(12,4,10,8,3,1,9,11,6,2,5,7)],cluster_rows = F,cluster_cols = F)

cellnumber <- as.data.frame.array(table(STS_pega_1@meta.data$Channel))
cellnumber$X <- rownames(cellnumber)
colnames(cellnumber)[1] <- "counts"
#cellnumber <- subset(cellnumber,X %in% samplenames_TNK)

cellnumber <- subset(cellnumber,counts > 100)

cellnumber_compare <- as.data.frame.array(table(STS_pega_1@meta.data$Channel,STS_pega_1@meta.data$celltype))
cellnumber_compare <- cellnumber_compare[cellnumber$X,]
rownames(cellnumber) <- cellnumber$X
cellnumber <- cellnumber[rownames(cellnumber_compare),]
cellnumber_compare_1 <- cellnumber_compare
for (i in 1:length(rownames(cellnumber_compare_1))) {
  cellnumber_compare_1[i,] <- cellnumber_compare_1[i,]/cellnumber$counts[i]
}
metadisease <- tumor_sample
rownames(metadisease) <- metadisease$Sample
metadisease_1 <- metadisease[rownames(cellnumber_compare),]
cellnumber_compare_1 <- cbind(cellnumber_compare_1,metadisease_1$Disease,metadisease_1$Malignancy)
colnames(cellnumber_compare_1)[c(18,19)] <- c("Disease","Malignancy")

diseasecounts <- as.data.frame(table(cellnumber_compare_1$Disease))

cellnumber_compare_1_1 <- subset(cellnumber_compare_1,Disease %in% as.character(subset(diseasecounts,Freq > 2)$Var1))
cellnumber_compare_1_1 <- rbind(subset(cellnumber_compare_1_1,Malignancy != "Peritumor"&Tumor > 0),
                                subset(cellnumber_compare_1_1,Malignancy == "Peritumor"))
plots <- plot_cellnumber_compare(cellnumber_compare_1_1,
                                 group_col = "Disease",
                                 mal_col = "Malignancy",
                                 start_col = 1,
                                 end_col = NULL,
                                 hide_ns = FALSE,
                                 step = 0.05,
                                 cols = cols[as.character(unique(cellnumber_compare_1_1$Disease)[order(unique(cellnumber_compare_1_1$Disease))])])

# 看第一个图
plots[[1]]+plots[[2]]
wrap_plots(plots)

#UMAP
metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
# Convert("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad", dest = "h5seurat", overwrite = F)
# 
# f <- H5File$new("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat", "r+")
# groups <- f$ls(recursive = TRUE)
# 
# for (name in groups$name[grepl("categories", groups$name)]) {
#   names <- strsplit(name, "/")[[1]]
#   names <- c(names[1:length(names) - 1], "levels")
#   new_name <- paste(names, collapse = "/")
#   f[[new_name]] <- f[[name]]
# }
# 
# for (name in groups$name[grepl("codes", groups$name)]) {
#   names <- strsplit(name, "/")[[1]]
#   names <- c(names[1:length(names) - 1], "values")
#   new_name <- paste(names, collapse = "/")
#   f[[new_name]] <- f[[name]]
#   grp <- f[[new_name]]
#   grp$write(args = list(1:grp$dims), value = grp$read() + 1)
# }
# 
# f$close_all()
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

tumor_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_tumor_5_6.xlsx","Sheet 1")
normal_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_normal_5_6.xlsx","Sheet 1")
barcodes <- rbind(tumor_barcodes,normal_barcodes)
STS_pega@meta.data$barcode <- rownames(STS_pega@meta.data)
STS_pega_1 <- subset(STS_pega,barcode %in% barcodes$barcode)
rownames(barcodes) <- barcodes$barcode
barcodes <- barcodes[colnames(STS_pega_1),]
STS_pega_1$group <- barcodes$group
DimPlot(STS_pega_1,group.by = "group",raster = F,reduction = "tsne")
table <- as.data.frame.array(table(STS_pega_1$louvain_labels,STS_pega_1$group))

data_normal_Group_celltype <- read.csv("/cluster/home/yflu/STS/pegasus/data_normal_Group_celltype.csv")
barcodes_celltype <- tumor_barcodes
barcodes_celltype$barcodekey <- barcodes_celltype$barcode
barcodes_celltype$Group <- barcodes_celltype$group
barcodes_celltype$Celltype <- barcodes_celltype$group
barcodes_celltype <- rbind(barcodes_celltype[,c(4:6)],data_normal_Group_celltype[,c(1:3)])
rownames(barcodes_celltype) <- barcodes_celltype$barcodekey
barcodes_celltype <- barcodes_celltype[colnames(STS_pega_1),]
STS_pega_1$celltype <- barcodes_celltype$Celltype

DimPlot(STS_pega_1,group.by = "celltype",raster = F,reduction = "tsne")
STS_pega_1$celltype <- factor(STS_pega_1$celltype,levels = unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)])
DimPlot(STS_pega_1,group.by = "celltype",raster = F,reduction = "tsne",
        cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(2)[c(1)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(6)[c(1:2)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(6)[c(1:2)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(2)[c(1)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1:4)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(9,10)])(2)[c(1)]))
DimPlot(STS_pega_1,group.by = "celltype",raster = F,reduction = "umap",
        cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(2)[c(1)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(6)[c(1:2)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(6)[c(1:2)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(2)[c(1)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1:4)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(9,10)])(2)[c(1)]))

anno_merged <- read.xlsx("/cluster/home/yflu/STS/pegasus/anno_merged.xlsx","Sheet2")

aurocs_Disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")

Disease_order <- as.data.frame(table(anno_merged$Disease))
Disease_order <- Disease_order[order(Disease_order$Freq,decreasing = T),]
Disease_order <- as.character(Disease_order$Var1)
aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
Diseasenames <- rownames(aurocs_Disease)[p$tree_row$order]
Diseasenames <- substr(Diseasenames,5,nchar(Diseasenames))

cols <- c(
  "Lipoblastoma"        = "#8DD3C7",
  "Synovial sarcoma"    = "#97D7C5",
  "Liposarcoma"         = "#A1DBC3",
  "Spindle cell tumor"  = "#ACDFC1",
  "Infantile fibrosarcoma" = "#B6E3BF",
  "Aggressive fibromatosis" = "#C0E7BD",
  "IMT"                 = "#CBEBBC",
  "ASPS"                = "#D5EFBA",
  "NF"                  = "#FFFFB3",
  "Schwannoma"          = "#F5F5B8",
  "Lymphangioma"        = "#FCCDE5",
  "Hemangioma"          = "#F0D1E1",
  "KHE"                 = "#E4D5DD",
  "MPNST"               = "#B3DE69",
  "Undifferentiated sarcoma" = "#80B1D3",
  "Angiosarcoma"        = "#FDB462",
  "RMS"                 = "#ADD8E6",
  "Pecoma"              = "#BC80BD",
  "EWS"            = "#BE8FBE",
  "MRT"                 = "#C09EBF",
  "Peritumor" = "#F18EAC"
)

anno_merged$Disease <- factor(anno_merged$Disease,names(cols))
anno_merged$Malignancy <- factor(anno_merged$Malignancy,c("Benign","Malignant","Peritumor"))

anno_sorted <- anno_merged %>%
  arrange(Malignancy, Disease, desc(scRNA), desc(Xenium), desc(WES))

samples <- anno_sorted$Sample

# 构造一个全 0 的矩阵，10 个基因 × n 个样本
expr_mat <- matrix(0, nrow = 10, ncol = length(samples))
rownames(expr_mat) <- paste0("Gene", 1:10)
colnames(expr_mat) <- samples

# 提取第 2-6 列作为注释
anno_col <- anno_sorted[, 4:8]
rownames(anno_col) <- anno_sorted$Sample  

pheatmap(expr_mat,
         annotation_col = anno_col,
         annotation_colors = list(Disease = cols),
         breaks = c(-1, 0, 1))  # 避免全零报错
metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

tumor_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_tumor_5_6.xlsx","Sheet 1")
normal_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_normal_5_6.xlsx","Sheet 1")
barcodes <- rbind(tumor_barcodes,normal_barcodes)
STS_pega@meta.data$barcode <- rownames(STS_pega@meta.data)
STS_pega_1 <- subset(STS_pega,barcode %in% barcodes$barcode)
rownames(barcodes) <- barcodes$barcode
barcodes <- barcodes[colnames(STS_pega_1),]
STS_pega_1$group <- barcodes$group
table <- as.data.frame.array(table(STS_pega_1$louvain_labels,STS_pega_1$group))

data_normal_Group_celltype <- read.csv("/cluster/home/yflu/STS/pegasus/data_normal_Group_celltype.csv")
barcodes_celltype <- tumor_barcodes
barcodes_celltype$barcodekey <- barcodes_celltype$barcode
barcodes_celltype$Group <- barcodes_celltype$group
barcodes_celltype$Celltype <- barcodes_celltype$group
barcodes_celltype <- rbind(barcodes_celltype[,c(4:6)],data_normal_Group_celltype[,c(1:3)])
rownames(barcodes_celltype) <- barcodes_celltype$barcodekey
barcodes_celltype <- barcodes_celltype[colnames(STS_pega_1),]
STS_pega_1$celltype <- barcodes_celltype$Celltype

STS_pega_1$celltype <- factor(STS_pega_1$celltype,levels = unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)])

tumor_sample <- anno_merged
meta <- STS_pega_1@meta.data

# 建立 Sample -> Disease 的映射表
sample2disease <- setNames(as.character(tumor_sample$Disease), as.character(tumor_sample$Sample))

# 按照 Channel 填充 Disease
meta$Disease <- sample2disease[as.character(meta$Channel)]

# 更新回 Seurat 对象
STS_pega_1@meta.data <- meta

meta <- STS_pega_1@meta.data

metadata <- read_h5ad("/cluster3/yflu/STS/TMA/adata_TMA_merged_new.h5ad")
TMA_merged <- LoadH5Seurat("/cluster3/yflu/STS/TMA/adata_TMA_merged_new.h5seurat",meta.data = FALSE, misc = FALSE)
TMA_merged@meta.data <- metadata$obs

map_vec <- c(
  "Monocytes/macrophages" = "Mono/macro",
  "Vascular endothelial cells" = "Endothelial",
  "NK cells" = "T/NK cells",
  "Plasmacytoid dendritic cells" = "DC",
  "Lymphatic endothelial cells" = "Endothelial",
  "T cells" = "T/NK cells",
  "Myocytes" = "Other cells",
  "Schwann cell" = "Other cells",
  "Myoblasts" = "Other cells",
  "Schwann cells" = "Other cells",
  "Proliferating cells" =  "Other cells"
)

# 一次性替换
STS_pega_1@meta.data$Celltype_united <- 
  plyr::mapvalues(STS_pega_1@meta.data$celltype,
                  from = names(map_vec),
                  to   = unname(map_vec),
                  warn_missing = FALSE)
STS_pega_1@meta.data$Celltype_united <- as.character(STS_pega_1@meta.data$Celltype_united)

STS_pega_1.sce <- as.SingleCellExperiment(STS_pega_1)
TMA_merged.sce <- as.SingleCellExperiment(TMA_merged)

samplelist <- list(STS_pega_1.sce,TMA_merged.sce)
names(samplelist) <- c("scRNA-seq","Xenium")

fused_data = mergeSCE(samplelist)

global_hvgs = variableGenes(dat = fused_data, exp_labels = fused_data$study_id)

#1C
aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = fused_data,
                        study_id = fused_data$study_id,
                        cell_type = fused_data$Celltype_united,
                        fast_version = TRUE)
saveRDS(aurocs,"/cluster3/yflu/STS/development/aurocs_SC_XE_250924.rds")

#1F
plot_cellnumber_compare <- function(df, group_col = "Disease", mal_col = "Malignancy",
                                    start_col = 1, end_col = NULL,
                                    hide_ns = TRUE, step = 0.05,
                                    cols = NULL, dodge_width = 0.8, box_width = 0.6) {
  if (is.null(end_col)) {
    end_col <- ncol(df) - 2  # 默认最后两列是 Disease / Malignancy
  }
  plot_list <- list()
  
  # 如果提供了颜色，保证 group_col 是 factor 且顺序一致
  if (!is.null(cols)) {
    df[[group_col]] <- factor(df[[group_col]], levels = names(cols))
  }
  
  for (i in start_col:end_col) {
    tmp <- df[, c(i, which(colnames(df) %in% c(group_col, mal_col)))]
    colnames(tmp)[1] <- "proportion"
    tmp[[mal_col]] <- as.character(tmp[[mal_col]])
    
    diseasenames <- sort(unique(tmp[[mal_col]]))
    my_comparisons <- combn(diseasenames, 2, simplify = FALSE)
    
    my_comparisons_sig <- list()
    for (j in 1:length(my_comparisons)) {
      por1 <- tmp[tmp[[mal_col]] == my_comparisons[[j]][1], "proportion"]
      por2 <- tmp[tmp[[mal_col]] == my_comparisons[[j]][2], "proportion"]
      if (length(por1) > 1 & length(por2) > 1) {
        test <- t.test(por1, por2)
        if (is.na(test$p.value)) test$p.value <- 1
        if (test$p.value < 0.05) {
          my_comparisons_sig <- append(my_comparisons_sig, list(my_comparisons[[j]]))
        }
      }
    }
    
    # 关键修改：interaction，确保 group 宽度一致
    P <- ggplot(tmp, aes(x = !!sym(mal_col), y = proportion, fill = !!sym(group_col))) +
      stat_boxplot(geom = "errorbar", width = 0.15, position = position_dodge(width = dodge_width)) +
      geom_boxplot(position = position_dodge(width = dodge_width), width = box_width, outlier.color = "white") +
      geom_jitter(shape = 21, size = 1.5, position = position_dodge(width = dodge_width)) +
      ggtitle(colnames(df)[i]) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )
    
    # 应用自定义配色
    if (!is.null(cols)) {
      P <- P + scale_fill_manual(values = cols, drop = FALSE)
    }
    
    # 显著性标注
    if (length(my_comparisons_sig) > 0) {
      ymax <- max(tmp$proportion, na.rm = TRUE)
      label_y <- ymax + (seq_along(my_comparisons_sig)) * (ymax * step)
      
      P <- P + stat_compare_means(
        comparisons = my_comparisons_sig,
        method = "t.test",
        label = "p.signif",
        hide.ns = hide_ns,
        label.y = label_y
      )
    }
    
    plot_list[[colnames(df)[i]]] <- P
  }
  return(plot_list)
}

anno_merged <- read.xlsx("/cluster/home/yflu/STS/pegasus/anno_merged.xlsx","Sheet2")

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')

Disease_order <- as.data.frame(table(anno_merged$Disease))
Disease_order <- Disease_order[order(Disease_order$Freq,decreasing = T),]
Disease_order <- as.character(Disease_order$Var1)

Diseasenames <- rownames(aurocs_disease)[p$tree_row$order]
Diseasenames <- substr(Diseasenames,5,nchar(Diseasenames))

cols <- c(
  "Lipoblastoma"        = "#8DD3C7",
  "Synovial sarcoma"    = "#97D7C5",
  "Liposarcoma"         = "#A1DBC3",
  "Spindle cell tumor"  = "#ACDFC1",
  "Infantile fibrosarcoma" = "#B6E3BF",
  "Aggressive fibromatosis" = "#C0E7BD",
  "IMT"                 = "#CBEBBC",
  "ASPS"                = "#D5EFBA",
  "NF"                  = "#FFFFB3",
  "Schwannoma"          = "#F5F5B8",
  "Lymphangioma"        = "#FCCDE5",
  "Hemangioma"          = "#F0D1E1",
  "KHE"                 = "#E4D5DD",
  "MPNST"               = "#B3DE69",
  "Undifferentiated sarcoma" = "#80B1D3",
  "Angiosarcoma"        = "#FDB462",
  "RMS"                 = "#ADD8E6",
  "Pecoma"              = "#BC80BD",
  "EWS/PNET"            = "#BE8FBE",
  "MRT"                 = "#C09EBF",
  "Peritumor" = "#F18EAC"
)

anno_merged$Disease <- factor(anno_merged$Disease,names(cols))
anno_merged$Malignancy <- factor(anno_merged$Malignancy,c("Benign","Malignant","Peritumor"))

anno_sorted <- anno_merged %>%
  arrange(Malignancy, Disease, desc(scRNA), desc(Xenium), desc(WES))

samples <- anno_sorted$Sample

# 构造一个全 0 的矩阵，10 个基因 × n 个样本
expr_mat <- matrix(0, nrow = 10, ncol = length(samples))
rownames(expr_mat) <- paste0("Gene", 1:10)
colnames(expr_mat) <- samples

# 提取第 2-6 列作为注释
anno_col <- anno_sorted[, 4:8]
rownames(anno_col) <- anno_sorted$Sample  

pheatmap(expr_mat,
         annotation_col = anno_col,
         annotation_colors = list(Disease = cols),
         breaks = c(-1, 0, 1))  # 避免全零报错

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

tumor_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_tumor_5_6.xlsx","Sheet 1")
normal_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_normal_5_6.xlsx","Sheet 1")
barcodes <- rbind(tumor_barcodes,normal_barcodes)
STS_pega@meta.data$barcode <- rownames(STS_pega@meta.data)
STS_pega_1 <- subset(STS_pega,barcode %in% barcodes$barcode)
rownames(barcodes) <- barcodes$barcode
barcodes <- barcodes[colnames(STS_pega_1),]
STS_pega_1$group <- barcodes$group
table <- as.data.frame.array(table(STS_pega_1$louvain_labels,STS_pega_1$group))

data_normal_Group_celltype <- read.csv("/cluster/home/yflu/STS/pegasus/data_normal_Group_celltype.csv")
barcodes_celltype <- tumor_barcodes
barcodes_celltype$barcodekey <- barcodes_celltype$barcode
barcodes_celltype$Group <- barcodes_celltype$group
barcodes_celltype$Celltype <- barcodes_celltype$group
barcodes_celltype <- rbind(barcodes_celltype[,c(4:6)],data_normal_Group_celltype[,c(1:3)])
rownames(barcodes_celltype) <- barcodes_celltype$barcodekey
barcodes_celltype <- barcodes_celltype[colnames(STS_pega_1),]
STS_pega_1$celltype <- barcodes_celltype$Celltype

DimPlot(STS_pega_1,group.by = "celltype",raster = T,reduction = "tsne")
STS_pega_1$celltype <- factor(STS_pega_1$celltype,levels = unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)])

tumor_sample <- anno_merged
meta <- STS_pega_1@meta.data

# 建立 Sample -> Disease 的映射表
sample2Disease <- setNames(as.character(tumor_sample$Disease), as.character(tumor_sample$Sample))

# 按照 Channel 填充 Disease
meta$Disease <- sample2Disease[as.character(meta$Channel)]
sample2Malignancy <- setNames(as.character(tumor_sample$Malignancy), as.character(tumor_sample$Sample))

meta$Malignancy <- sample2Malignancy[as.character(meta$Channel)]

# 更新回 Seurat 对象
STS_pega_1@meta.data <- meta

meta <- STS_pega_1@meta.data

# 统计每个 Disease - celltype 的数量和比例
cell_counts <- meta %>%
  group_by(Channel, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Channel) %>%
  mutate(prop = n / sum(n))

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

cell_counts <- subset(cell_counts,Channel %in% rownames(anno_sample_cluster_extended))

library(dplyr)
stromal_types <- unique(as.character(cell_counts$celltype))[c(2:6,14,16,17)]

stromal_summary <- cell_counts %>%
  filter(as.character(celltype) %in% stromal_types) %>%
  group_by(Channel) %>%
  summarise(
    celltype = "stromal",
    n = sum(n),
    prop = sum(prop),
    .groups = "drop"
  )

stromal_summary

cell_counts_new <- bind_rows(
  cell_counts,
  stromal_summary
)

cell_counts_tumor_fibro <- subset(cell_counts_new,celltype %in% c("Tumor","stromal"))
cell_counts_tumor_fibro <- as.data.frame(cell_counts_tumor_fibro)
cell_counts_tumor_fibro$celltype <- as.character(cell_counts_tumor_fibro$celltype)
cell_counts_tumor_fibro$Channel <- as.character(cell_counts_tumor_fibro$Channel)

cell_counts_tumor_fibro_complete <- cell_counts_tumor_fibro %>%
  complete(Channel, celltype, fill = list(n = 0, prop = 0))

cell_counts_tumor_fibro_wide <- cell_counts_tumor_fibro_complete %>%
  pivot_wider(
    id_cols = Channel,
    names_from = celltype,
    values_from = prop,
    values_fill = 0
  )

cell_counts_tumor_fibro_wide <- as.data.frame(cell_counts_tumor_fibro_wide)
rownames(cell_counts_tumor_fibro_wide) <- cell_counts_tumor_fibro_wide$Channel
cell_counts_tumor_fibro_wide <- cbind(cell_counts_tumor_fibro_wide,anno_sample_cluster_extended[rownames(cell_counts_tumor_fibro_wide),])

colnames(cell_counts_tumor_fibro_wide)[c(1,4)] <- c("Samples","Diseases")

ggscatter(
  cell_counts_tumor_fibro_wide, 
  x = "stromal", 
  y = "Tumor", 
  color = "Diseases",
  size = 5,
  shape = "Malignancy",
  xlab = "stromal_prop", 
  ylab = "Tumor_prop"
  
) +
  # Add overall regression line and confidence interval
  geom_smooth(
    method = "lm", 
    aes(x = stromal, y = Tumor), 
    color = "black", 
    se = TRUE
  ) +
  # Add overall correlation coefficient
  stat_cor(
    aes(x = stromal, y = Tumor), 
    method = "spearman", 
  ) +
  scale_color_manual(values = cols) +
  theme_classic()

ggscatter(
  subset(cell_counts_tumor_fibro_wide,Malignancy == "Malignant"), 
  x = "stromal", 
  y = "Tumor", 
  color = "Diseases",
  size = 5,
  shape = "Malignancy",
  xlab = "stromal_prop", 
  ylab = "Tumor_prop"
  
) +
  # Add overall regression line and confidence interval
  geom_smooth(
    method = "lm", 
    aes(x = stromal, y = Tumor), 
    color = "black", 
    se = TRUE
  ) +
  # Add overall correlation coefficient
  stat_cor(
    aes(x = stromal, y = Tumor), 
    method = "spearman", 
  ) +
  scale_color_manual(values = cols) +
  theme_classic()

plot_gene_average_boxplot <- function(expr_df, sample_anno, target_diseases, genes,
                                      disease_order, disease_colors,
                                      test = "t") {
  library(dplyr)
  library(ggplot2)
  
  # -----------------------
  # 1. 样本分组
  # -----------------------
  sample_anno <- sample_anno %>%
    mutate(Group = ifelse(Diseases %in% target_diseases, "Target", "Others"))
  
  group_samples <- as.character(sample_anno$Samples)
  
  # -----------------------
  # 2. 计算每个样本的平均表达
  # -----------------------
  expr_sub <- expr_df[genes, group_samples, drop=FALSE]
  expr_avg <- colMeans(expr_sub, na.rm = TRUE)
  
  plot_df <- data.frame(
    Sample = names(expr_avg),
    Expression = expr_avg,
    Group = sample_anno$Group[match(names(expr_avg), sample_anno$Samples)],
    Diseases = sample_anno$Diseases[match(names(expr_avg), sample_anno$Samples)]
  )
  
  plot_df$Group <- factor(plot_df$Group, levels = c("Target", "Others"))
  
  # -----------------------
  # 3. 统计检验
  # -----------------------
  group1 <- plot_df$Expression[plot_df$Group=="Target"]
  group2 <- plot_df$Expression[plot_df$Group=="Others"]
  
  if (test=="t") {
    p_val <- t.test(group1, group2)$p.value
  } else if (test=="wilcox") {
    p_val <- wilcox.test(group1, group2)$p.value
  } else {
    stop("test must be 't' or 'wilcox'")
  }
  
  get_sig <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  sig_label <- get_sig(p_val)
  y_max <- max(plot_df$Expression) * 1.15
  
  # -----------------------
  # 4. 绘图
  # -----------------------
  box_fill <- c("Target" = "#E64B35", "Others" = "#4DBBD5")
  names(disease_colors) <- disease_order
  
  p <- ggplot(plot_df, aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(fill = Diseases), shape = 21, size = 2, width = 0.2, alpha = 0.8, color = "gray") +
    geom_segment(aes(x = 1, xend = 2, y = y_max, yend = y_max), inherit.aes = FALSE) +
    geom_text(aes(x = 1.5, y = y_max*1.02, label = sig_label), inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c(box_fill, disease_colors)) +
    theme_classic() +
    labs(title = paste("Average expression of", paste(genes, collapse=", "),
                       "in target vs others"))
  
  print(p)
  
  return(p)
}
colnames(cell_counts_tumor_fibro_wide)[c(1,4)] <- c("Samples","Diseases")
disease_order <- names(cols)
p1 = plot_gene_average_boxplot(t(cell_counts_tumor_fibro_wide[,c(2,3)]),
                               cell_counts_tumor_fibro_wide, 
                               target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                               genes = "stromal",
                               disease_order = disease_order,
                               disease_colors = cols,
                               test="t")

p2 = plot_gene_average_boxplot(t(cell_counts_tumor_fibro_wide[,c(2,3)]),
                               cell_counts_tumor_fibro_wide, 
                               target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                               genes = "Tumor",
                               disease_order = disease_order,
                               disease_colors = cols,
                               test="t")
p2+p1

plot_gene_pair_boxplot <- function(expr_df, sample_anno, genes,
                                   disease_order, disease_colors,
                                   test = "t") {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  if (length(genes) != 2) stop("genes must contain exactly two gene names for paired comparison.")
  
  # -----------------------
  # 1. 提取两个基因的表达
  # -----------------------
  group_samples <- intersect(colnames(expr_df), sample_anno$Samples)
  if (length(group_samples) == 0)
    stop("No overlapping sample names between expr_df and sample_anno$Samples.")
  
  expr_sub <- expr_df[genes, group_samples, drop = FALSE]
  plot_df <- as.data.frame(t(expr_sub))
  colnames(plot_df) <- genes
  plot_df$Sample <- rownames(plot_df)
  
  # 加入注释信息
  plot_df <- plot_df %>%
    left_join(sample_anno[, c("Samples", "Diseases")], by = c("Sample" = "Samples"))
  
  # -----------------------
  # 2. 转成长格式
  # -----------------------
  plot_df_long <- pivot_longer(plot_df, cols = all_of(genes),
                               names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = factor(Gene, levels = genes))
  
  # -----------------------
  # 3. 配对统计检验
  # -----------------------
  g1 <- plot_df[[genes[1]]]
  g2 <- plot_df[[genes[2]]]
  
  if (test == "t") p_val <- t.test(g1, g2, paired = TRUE)$p.value
  else if (test == "wilcox") p_val <- wilcox.test(g1, g2, paired = TRUE)$p.value
  else stop("test must be 't' or 'wilcox'")
  
  get_sig <- function(p) {
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  }
  sig_label <- get_sig(p_val)
  y_max <- max(plot_df_long$Expression, na.rm = TRUE) * 1.15
  
  # -----------------------
  # 4. 准备配对连线数据
  # -----------------------
  connect_df <- plot_df_long %>%
    select(Sample, Gene, Expression) %>%
    pivot_wider(names_from = Gene, values_from = Expression)
  
  xstart_col <- 1
  xend_col <- 2
  ystart_col <- genes[1]
  yend_col <- genes[2]
  
  # -----------------------
  # 5. 绘图
  # -----------------------
  box_fill <- c(setNames("#E64B35", genes[1]), setNames("#4DBBD5", genes[2]))
  names(disease_colors) <- disease_order
  
  p <- ggplot(plot_df_long, aes(x = Gene, y = Expression)) +
    # boxplot
    geom_boxplot(aes(fill = Gene), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    # 配对连线
    geom_segment(data = connect_df,
                 aes(x = xstart_col, xend = xend_col,
                     y = !!as.name(ystart_col), yend = !!as.name(yend_col)),
                 color = "gray50", alpha = 0.6, linewidth = 0.4) +
    # 点
    geom_point(aes(fill = Diseases), shape = 21, size = 2, color = "gray") +
    # 显著性标注
    geom_segment(aes(x = 1, xend = 2, y = y_max, yend = y_max), inherit.aes = FALSE) +
    geom_text(aes(x = 1.5, y = y_max * 1.02, label = sig_label),
              inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c(box_fill, disease_colors)) +
    scale_x_discrete(limits = genes) +
    theme_classic(base_size = 14) +
    labs(title = paste("Paired expression comparison:", genes[1], "vs", genes[2]),
         subtitle = paste0("Paired ", test, " test, P = ", signif(p_val, 3)))
  
  print(p)
  invisible(p)
}

p1 = plot_gene_pair_boxplot(t(subset(cell_counts_tumor_fibro_wide,Malignancy == "Malignant")[,c(2,3)]),
                            cell_counts_tumor_fibro_wide, 
                            genes = c("Tumor","stromal"),
                            disease_order = disease_order,
                            disease_colors = cols,
                            test="t")

p2 = plot_gene_pair_boxplot(t(subset(cell_counts_tumor_fibro_wide,Malignancy == "Benign")[,c(2,3)]),
                            cell_counts_tumor_fibro_wide, 
                            genes = c("Tumor","stromal"),
                            disease_order = disease_order,
                            disease_colors = cols,
                            test="t")
p1+p2

#1E
metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

tumor_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_tumor_5_6.xlsx","Sheet 1")
normal_barcodes <- openxlsx::read.xlsx("/cluster/home/yflu/STS/separated/CNVSCORE/group_all_normal_5_6.xlsx","Sheet 1")
barcodes <- rbind(tumor_barcodes,normal_barcodes)
STS_pega@meta.data$barcode <- rownames(STS_pega@meta.data)
STS_pega_1 <- subset(STS_pega,barcode %in% barcodes$barcode)
rownames(barcodes) <- barcodes$barcode
barcodes <- barcodes[colnames(STS_pega_1),]
STS_pega_1$group <- barcodes$group
table <- as.data.frame.array(table(STS_pega_1$louvain_labels,STS_pega_1$group))

data_normal_Group_celltype <- read.csv("/cluster/home/yflu/STS/pegasus/data_normal_Group_celltype.csv")
barcodes_celltype <- tumor_barcodes
barcodes_celltype$barcodekey <- barcodes_celltype$barcode
barcodes_celltype$Group <- barcodes_celltype$group
barcodes_celltype$Celltype <- barcodes_celltype$group
barcodes_celltype <- rbind(barcodes_celltype[,c(4:6)],data_normal_Group_celltype[,c(1:3)])
rownames(barcodes_celltype) <- barcodes_celltype$barcodekey
barcodes_celltype <- barcodes_celltype[colnames(STS_pega_1),]
STS_pega_1$celltype <- barcodes_celltype$Celltype

STS_pega_1$celltype <- factor(STS_pega_1$celltype,levels = unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)])

sample_counts <- as.data.frame.array(table(STS_pega_1$Channel,STS_pega_1$celltype))
sample_counts <- sample_counts[rownames(anno_sample_cluster_extended),]

sample_counts$sum <- rowSums(sample_counts)

prop <- sample_counts
prop[, colnames(prop) != "sum"] <- prop[, colnames(prop) != "sum"] / prop$sum

prop <- prop[rownames(anno_sample_cluster_extended), ]

library(tidyverse)

df_long <- prop %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "CellType",
    values_to = "Proportion"
  ) %>%
  left_join(
    anno_sample_cluster_extended %>%
      rownames_to_column("Sample"),
    by = "Sample"
  )
diff_df <- df_long %>%
  group_by(CellType, Malignancy) %>%
  summarise(mean_prop = mean(Proportion), .groups = "drop") %>%
  pivot_wider(
    names_from = Malignancy,
    values_from = mean_prop
  ) %>%
  mutate(diff = Malignant - Benign)

stat_df <- df_long %>%
  group_by(CellType) %>%
  summarise(
    p_value = t.test(
      Proportion ~ Malignancy,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

res <- diff_df %>%
  left_join(stat_df, by = "CellType") %>%
  arrange(FDR)
res <- res %>%
  mutate(sig = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    TRUE ~ ""
  ))
res <- res[-12,]
diff_mat <- res %>%
  select(CellType, diff) %>%
  column_to_rownames("CellType") %>%
  as.matrix()
res <- as.data.frame(res)
rownames(res) <- res$CellType
res <- res[as.character(unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)]),]

library(pheatmap)

max_abs <- max(abs(diff_mat))

sig_label <- ifelse(
  res$FDR < 0.001, "***",
  ifelse(res$FDR < 0.01, "**",
         ifelse(res$FDR < 0.05, "*", ""))
)

sig_mat <- matrix(
  sig_label,
  nrow = nrow(diff_mat),
  ncol = ncol(diff_mat)
)
rownames(sig_mat) <- rownames(diff_mat)
colnames(sig_mat) <- colnames(diff_mat)

pheatmap(
  diff_mat[as.character(unique(STS_pega_1$celltype)[c(8,3,14,1,17,12,6,16,11,9,2,7,10,4,5,13,15)]),c(1,1)],
  cluster_rows = F,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  breaks = seq(-max_abs, max_abs, length.out = 101),
  border_color = NA,
  display_numbers = cbind(sig_mat,sig_mat),
  number_color = "black",
  fontsize_row = 10,
  main = "Difference in cell proportions (Malignant − Benign)"
)