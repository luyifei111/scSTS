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
library(ggpubr)
library(ggcorrplot)
library(tidyr)

metadata <- read_h5ad("/cluster3/yflu/STS/development/data_development_umap.h5ad")
STS.pega.development <- LoadH5Seurat("/cluster3/yflu/STS/development/data_development_umap.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.development@meta.data <- metadata$obs
Idents(STS.pega.development) <- STS.pega.development$louvain_labels
celltype <- read.xlsx("/cluster3/yflu/STS/development/integrated/celltype.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS.pega.development)
STS.pega.development <- RenameIdents(STS.pega.development, new.cluster.ids)
STS.pega.development$celltype_new <- Idents(STS.pega.development)

aurocs_disease_250210 <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
order <- names(table(STS.pega.development$celltype_new))[c(1,7,5,18,8,12,4,10,6,3,13,2,16,9,17,11,14,15)]
order <- order[c(2:4,14,15,12,13,6,11,7:10,1,5,16:18)]
aurocs_disease <- readRDS("aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
#plotHeatmap(aurocs_1)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

samplenames <- read.csv("/cluster/home/yflu/STS/pegasus/STS_sample_select_24.5.5.csv")
sample_tumors <- unique(anno_sample_cluster_extended$Disease)
samplenames_tumor <- sample_tumors
samplenames_tumor_disease <- sample_tumors
samplenames_tumor_disease <- as.data.frame(samplenames_tumor_disease)
rownames(samplenames_tumor_disease) <- sample_tumors
colnames(samplenames_tumor_disease) <- "disease"
rownames(samplenames_tumor_disease) <- paste("STS|",rownames(samplenames_tumor_disease),sep = "")

p = pheatmap(aurocs_disease[,paste("Development|",rev(order),sep = "")],clustering_distance_rows = 'euclidean',cluster_cols = F)

rb.genes <- rownames(STS.pega.development)[grep("^RP[SL]",rownames(STS.pega.development))]
mt.genes <- rownames(STS.pega.development)[grep("^MT-",rownames(STS.pega.development))]

intersect_genes <- readRDS("/cluster3/yflu/STS/development/intersect_genes.rds")

i=1
markers <- openxlsx::read.xlsx("/cluster3/yflu/STS/development/development_markers_250213.xlsx",paste(order[i],"|up",sep = ""))
markers <- markers[order(markers$log2FC,decreasing = T),]
markers <- subset(markers,!(featurekey %in% rb.genes))
markers <- subset(markers,!(featurekey %in% mt.genes))
markers <- subset(markers,featurekey %in% intersect_genes)
markers_select <- markers$featurekey[1:2]

order_1 <- order
order_1[2] <- "Muscle_adipose MSC"

for (i in 2:length(order_1)) {
  markers <- openxlsx::read.xlsx("/cluster3/yflu/STS/development/development_markers_250213.xlsx",paste(order_1[i],"|up",sep = ""))
  markers <- markers[order(markers$log2FC,decreasing = T),]
  markers <- subset(markers,!(featurekey %in% markers_select))
  markers <- subset(markers,!(featurekey %in% rb.genes))
  markers <- subset(markers,!(featurekey %in% mt.genes))
  markers <- subset(markers,featurekey %in% intersect_genes)
  markers_select_1 <- markers$featurekey[1:2]
  markers_select <- c(markers_select,markers_select_1)
  print(i)
}

STS.pega.development$celltype_new <- factor(STS.pega.development$celltype_new,levels = order)
P = DotPlot(STS.pega.development,features = markers_select,group.by = "celltype_new",cols = c("#4575B4", "#D73027"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
P

data <- P$data
order_id <- rev(levels(data$id))
data$id <- as.character(data$id)
data$id <- factor(data$id,levels = order_id)
levels(data$features.plot) <- levels(data$features.plot)

diseasenames <- unique(as.character(STS.pega.development$celltype_new))
diseasenames <- as.data.frame(cbind(diseasenames,diseasenames))
colnames(diseasenames) <- c("Var1","Var2")
rownames(diseasenames) <- paste("Development|",diseasenames$Var1,sep="")
diseasenames <- diseasenames[paste("Development|",order,sep = ""),]
names(cols) <- diseasenames$Var1
data$cols <- as.character(data$id)

for (i in 1:648) {
  data$cols[i] <- cols[data$id[i]]
}
ggplot(data,aes(x=features.plot,y=id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_tile(data=data, aes(x=-0.5, y=id, fill=id), width=1, height=1)+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))+
  labs(x=NULL,y=NULL)

library(circlize)
library(RColorBrewer)

plot_similarity_circos_final <- function(df, group_colors, celltype_colors, group_order=NULL,
                                         plot_title=NULL, show_colorbar=TRUE){
  
  library(circlize)
  library(RColorBrewer)
  
  names(celltype_colors) <- sub("^.*\\|", "", names(celltype_colors))
  
  if(!is.null(group_order)){
    df$group <- factor(df$group, levels = group_order)
    df <- df[order(df$group, match(df$cell_type, unique(df$cell_type))), ]
  }
  
  n <- nrow(df)
  
  sim_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
  sim_breaks <- seq(0, 1, length.out = length(sim_colors))
  sim_col_fun <- colorRamp2(sim_breaks, sim_colors)
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after = rep(0.3, n),
    cell.padding = c(0,0,0,0),
    track.margin = c(0.002, 0.002)
  )
  circos.initialize(factors=df$cell_type, xlim=cbind(rep(0,n),rep(1,n)))
  
  ## ===== Outer: celltype ring =====
  circos.trackPlotRegion(df$cell_type, rep(0,n), ylim=c(0,1), track.height=0.06,
                         bg.border=NA,    # ✅ 恢复：不绘制外框线
                         panel.fun=function(x,y){
                           nm <- get.cell.meta.data("sector.index")
                           col <- celltype_colors[nm]; if(is.na(col)) col <- "grey80"
                           xl <- get.cell.meta.data("xlim")
                           circos.rect(xl[1],0,xl[2],1,col=col,border=NA)
                         })
  
  ## ===== Group ring =====
  circos.trackPlotRegion(df$cell_type, rep(0,n), ylim=c(0,1), track.height=0.1,
                         bg.border=NA,    # ✅ 恢复默认
                         panel.fun=function(x,y){
                           nm <- get.cell.meta.data("sector.index")
                           i <- which(df$cell_type==nm)
                           xl <- get.cell.meta.data("xlim")
                           circos.rect(xl[1],0,xl[2],1,col=group_colors[df$group[i]],border=NA)
                         })
  
  ## ===== Main similarity ring =====
  circos.trackPlotRegion(df$cell_type, rep(0.5,n), ylim=c(0,1), track.height=0.7,
                         bg.border="grey70",   # ✅ 只给本圈加灰框线
                         panel.fun=function(x,y){
                           nm <- get.cell.meta.data("sector.index")
                           i <- which(df$cell_type==nm)
                           
                           circos.rect(0,0,1,1, col=adjustcolor(group_colors[df$group[i]], alpha=0.15), border=NA)
                           
                           left <- 0.22; right <- 0.78; bottom <- 0
                           top <- df$similarity[i]
                           col <- sim_col_fun(top)
                           
                           circos.polygon(c(left,right,0.5),c(bottom,bottom,top),col=col,border=NA)
                           
                           circos.points(0.5, top, 
                                         pch=16, 
                                         cex=2, 
                                         col=col)
                         })
  
  ## ===== labels =====
  for(i in seq_len(n)){
    circos.text(0.5,1.95, df$cell_type[i],
                sector.index=df$cell_type[i],
                facing="clockwise",
                adj=c(0.5,0),
                cex=0.7)
  }
  
  ## ===== group legend =====
  legend("topright", legend=names(group_colors),
         col=group_colors, pch=16, cex=0.8, bty="n")
  
  ## ===== title =====
  if(!is.null(plot_title)) title(plot_title,cex.main=1.5,font.main=2)
  
  ## ===== color bar =====
  if(show_colorbar){
    par(xpd=TRUE)
    usr <- par("usr")
    y_offset <- usr[3] - 0.22
    
    rect_x <- seq(0,1,length.out=length(sim_colors)+1)
    for(i in seq_along(sim_colors)){
      rect(rect_x[i],y_offset,rect_x[i+1],y_offset+0.04,
           col=sim_colors[i],border=NA,xpd=TRUE)
    }
    axis(1,at=seq(0,1,0.2),labels=seq(0,1,0.2),pos=y_offset,cex.axis=0.8)
    mtext("Similarity",side=1,line=3)
  }
}

aurocs_disease_250210 <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")

diseases <- p$tree_row$labels[p$tree_row$order]

group_disease <- as.data.frame(p$tree_row$labels[p$tree_row$order])
colnames(group_disease) <- "disease"
group_disease$group <- c(rep("Fibro",7),'ASPS',rep("Peripheral neural",2),rep("Endothelial",3),rep("Undifferentiated",2),"Endothelial","RMS",rep("Neural-like",3))

i=1
group_seleceted <- subset(group_disease,group == unique(group_disease$group)[i])

df <- data.frame(
  cell_type = order,
  group = c(rep("Mesenchymal",3),rep("Supporting",2),rep("Endothelial",2),"Mesenchymal",rep("Myogenic",5),rep("Mesenchymal",2),rep("Neural",3)),
  similarity = as.numeric(colMeans(aurocs_disease[group_seleceted$disease,paste("Development|",order,sep = "")]))
)

# 组颜色
group_colors <- c(
  Mesenchymal=colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1)],
  Supporting=colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1)],
  Endothelial=colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1)],
  Myogenic=colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(20)[c(1)],
  Neural=colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(12)[c(1)]
)

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(20)[c(1:5)],
         colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(9,10)],
         colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(12)[c(1:3)])
names(cols) <- paste("Development|",order,sep = "")

group_order <- c("Mesenchymal","Supporting","Endothelial","Myogenic","Neural")


plot_similarity_circos_final(df, group_colors, cols, group_order,
                             plot_title = paste(unique(group_disease$group)[i],"aurocs circos"))
paste(unique(group_disease$group)[i],"aurocs circos.pdf")

plot_compare_feature_two_groups <- function(
    seurat_obj,
    features,
    group_col,
    groupA_levels,
    data_type = c("meta", "assay"),
    assay_name = "RNA",
    feature_name = NULL,
    line_height = 1.0,
    color_groupA = "#E64B35",  # 新增参数
    color_groupB = "#4DBBD5"   # 新增参数
){
  library(ggplot2)
  library(dplyr)
  
  data_type <- match.arg(data_type)
  
  # ---- 1. 提取 feature 数据 ----
  if(data_type == "meta"){
    mat <- seurat_obj@meta.data[, features, drop = FALSE]
  } else {
    mat <- FetchData(seurat_obj, vars = features, assay = assay_name)
  }
  
  # 确保 numeric
  mat <- as.data.frame(lapply(mat, function(x){
    if(is.factor(x)) x <- as.character(x)
    if(is.logical(x)) x <- as.numeric(x)
    as.numeric(x)
  }))
  
  # 多 feature 取平均
  value <- if(ncol(mat) > 1) rowMeans(mat, na.rm = TRUE) else mat[[1]]
  
  # ---- 2. 分组 ----
  group <- seurat_obj@meta.data[[group_col]]
  group_binary <- ifelse(group %in% groupA_levels, "GroupA", "GroupB")
  group_binary <- factor(group_binary, levels = c("GroupA","GroupB"))
  
  df <- data.frame(value, group_binary) |> na.omit()
  
  # ---- 3. Wilcoxon test ----
  p <- wilcox.test(value ~ group_binary, data = df)$p.value
  p_tag <- ifelse(p < 0.001,"***", ifelse(p < 0.01,"**", ifelse(p < 0.05,"*","ns")))
  
  # ---- 4. 密度（原始范围） ----
  densA <- density(df$value[df$group_binary=="GroupA"])
  densB <- density(df$value[df$group_binary=="GroupB"])
  
  dA <- data.frame(x = densA$x, y = densA$y, group = "GroupA")
  dB <- data.frame(x = densB$x, y = densB$y, group = "GroupB")
  dens_df <- rbind(dA, dB)
  
  min_val <- min(dens_df$x)
  
  # ---- 5. 中位数 ----
  mA <- median(df$value[df$group_binary=="GroupA"])
  mB <- median(df$value[df$group_binary=="GroupB"])
  
  # ---- 6. bar 和点固定位置 ----
  med_pts <- data.frame(
    x = c(mA, mB),
    y = c(-max(dens_df$y)*0.05, -max(dens_df$y)*0.1),
    group = c("GroupA","GroupB")
  )
  
  bar_df <- data.frame(
    x_start = min_val,
    x_end = c(mA, mB),
    y = c(-max(dens_df$y)*0.05, -max(dens_df$y)*0.1),
    group = c("GroupA","GroupB")
  )
  
  # ---- 7. 图标题 ----
  if(is.null(feature_name)){
    feature_name <- if(length(features)==1) features else paste0("Mean(", paste(features, collapse=","), ")")
  }
  
  # ---- 8. 绘图 ----
  p_plot <- ggplot() +
    geom_line(data = dens_df, aes(x=x, y=y, color=group), size = 1.2) +
    geom_segment(data = bar_df,
                 aes(x=x_start, xend=x_end, y=y, yend=y, color=group),
                 linewidth = 2, lineend = "round") +
    geom_point(data = med_pts, aes(x=x, y=y, color=group), size = 3) +
    annotate("text", x = median(c(mA,mB)), y = max(dens_df$y)*0.95,
             label = paste0("p = ", signif(p,2)," ",p_tag), size = 5) +
    theme_classic() +
    labs(title = feature_name,
         x = "Expression value",
         y = "Density (shared)") +
    scale_color_manual(values = c("GroupA" = color_groupA, "GroupB" = color_groupB)) +
    theme(legend.position="none")
  
  return(p_plot)
}

TMA_merged <- readRDS("/cluster3/yflu/STS/TMA/TMA_merged_origin_score.rds")
fits = readRDS("/cluster3/yflu/STS/development/regression_sim/lrFoetalClustersV4_development_TMA.RDS")
log_regression_250916 <- readRDS("/cluster3/yflu/STS/TMA/log_regression_250916.rds")

pp = do.call(cbind,log_regression_250916)
colnames(pp) = names(fits)
pp = (1+exp(-pp))**-1
colnames(pp) <- gsub("[ /-]", ".", colnames(pp))
#df = cbind(df,pp[rownames(df),])

#rownames(STS.pega.tumor@meta.data)[1:10]
#rownames(pp)[1:10]

TMA_merged@meta.data <- cbind(TMA_merged@meta.data,pp)
TMA_merged_tumor <- subset(TMA_merged,Celltype_united == "Tumor cells")

p1 = plot_compare_feature_two_groups(
  TMA_merged_tumor,
  features = c("Ganglion.cells","Neural.progenitors","Fetal.fibroblasts"),
  group_col = "Disease",
  groupA_levels = levels(TMA_merged_tumor$Disease)[c(18:20)], # group A
  data_type = "meta",
  color_groupA = "#E64B35",
  color_groupB = "#4DBBD5"
)

p2 = plot_compare_feature_two_groups(
  TMA_merged_tumor,
  features = c("Ganglion cells","Neural progenitors","Fetal fibroblasts"),
  group_col = "Disease",
  groupA_levels = levels(TMA_merged_tumor$Disease)[c(18:20)], # group A
  data_type = "assay",
  assay_name = "AUCell",
  color_groupA = "#E64B35",
  color_groupB = "#4DBBD5"
)
p1 / p2

TMA_merged_tumor_score_avag <- AverageExpression(TMA_merged_tumor,features = rownames(TMA_merged_tumor@assays$AUCell),assays = "AUCell",group.by = "Sample")
TMA_merged_tumor_score_avag <- TMA_merged_tumor_score_avag$AUCell
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

plot_density_with_peak_general_disease_scaled <- function(
    seurat_obj,
    features,
    color_vec,
    disease,
    data_type = c("meta", "assay"),
    assay_name = "RNA",
    sample_frac = 0.1,
    seed = 123,
    line_height = 0.3
) {
  library(ggplot2)
  library(reshape2)
  library(ggbeeswarm)
  library(dplyr)
  
  data_type <- match.arg(data_type)
  
  # 1. 获取数据 ----------------------------------------------------------
  if (data_type == "meta") {
    data_all <- seurat_obj@meta.data[, features, drop = FALSE]
  } else if (data_type == "assay") {
    data_all <- GetAssayData(seurat_obj, assay = assay_name, slot = "data")[features, , drop = FALSE]
    data_all <- as.data.frame(t(as.matrix(data_all)))
  }
  
  data_all$cell_id <- rownames(data_all)
  data_long <- melt(data_all, id.vars = "cell_id", variable.name = "feature", value.name = "value")
  
  # 2. 对每个 feature 单独 scale（使用全部细胞）-------------------------
  data_long <- data_long %>%
    group_by(feature) %>%
    mutate(value_scaled = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    ungroup()
  
  # 3. 只保留指定 disease 的细胞用于绘图 ------------------------------
  cell_use <- rownames(seurat_obj)[seurat_obj$Disease %in% disease]
  data_plot_all <- data_long %>% filter(cell_id %in% cell_use)
  
  # 4. 计算峰值和中心（基于所有绘图细胞）-------------------------------
  peak_df <- data_plot_all %>%
    group_by(feature) %>%
    summarise(
      npoints = n(),
      peak = if (n() >= 2) {
        d <- density(value_scaled)
        d$x[which.max(d$y)]
      } else {
        median(value_scaled, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  # 5. factor 顺序和颜色匹配 -------------------------------------------
  feature_levels <- features[features %in% unique(data_plot_all$feature)]
  data_plot_all$feature <- factor(data_plot_all$feature, levels = feature_levels)
  peak_df$feature <- factor(peak_df$feature, levels = feature_levels)
  color_vec_plot <- color_vec[names(color_vec) %in% feature_levels]
  
  # 6. 可视化前再采样 ---------------------------------------------------
  set.seed(seed)
  data_plot_sample <- data_plot_all %>%
    group_by(feature) %>%
    sample_frac(sample_frac, replace = FALSE) %>%
    ungroup()
  
  # 7. 获取每个 feature 的 ycenter（采样后估计）
  temp_plot <- ggplot(data_plot_sample, aes(x = value_scaled, y = feature)) +
    geom_quasirandom(groupOnX = FALSE, varwidth = TRUE, orientation = "y")
  gb <- ggplot_build(temp_plot)
  points_df <- gb$data[[1]]
  ycenter_df <- points_df %>%
    group_by(group) %>%
    summarise(ycenter = median(y), .groups = "drop")
  peak_df <- left_join(peak_df, ycenter_df, by = c("feature" = "group"))
  
  # 8. 绘图 --------------------------------------------------------------
  p <- ggplot(data_plot_sample, aes(x = value_scaled, y = feature, color = feature)) +
    geom_quasirandom(groupOnX = FALSE, varwidth = TRUE, alpha = 0.5, size = 0.5, orientation = "y") +
    geom_segment(data = peak_df,
                 aes(x = peak, xend = peak, y = ycenter - line_height, yend = ycenter + line_height),
                 color = "black", size = 1) +
    scale_color_manual(values = color_vec_plot) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          legend.position = "none") +
    labs(x = "Scaled expression")
  
  return(p)
}

old_names <- colnames(TMA_merged_tumor@meta.data)[31:48]

# 提取新名字
new_names <- cellnames$Var1

# 统一格式以便匹配（去掉点号、空格、连字符等）
normalize <- function(x) {
  gsub("[[:punct:] ]", "", x)  # 去掉所有标点和空格，包括 / . - _
}

# 创建匹配关系
match_idx <- match(normalize(old_names), normalize(new_names))

# 检查匹配情况
data.frame(old_names, matched_to = new_names[match_idx])

colnames(TMA_merged_tumor@meta.data)[31:48] <- new_names[match_idx]
colnames(TMA_merged_tumor@meta.data)[31:48] <- paste("Development|",colnames(TMA_merged_tumor@meta.data)[31:48],sep="")

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(20)[c(1:5)],
         colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(9,10)],
         colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(12)[c(1:3)])
names(cols) <- paste("Development|",order,sep = "")

meta_cols <- rev(names(cols)[c(11,16,10,12)])  # 对应你提供的 celltype 列

p1 = plot_density_with_peak_by_disease(
  seurat_obj = TMA_merged_tumor,
  features = meta_cols,
  color_vec = cols,
  disease = "RMS",
  data_type = "meta",
  sample_frac = 0.01
)

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(20)[c(1:5)],
         colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(9,10)],
         colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(12)[c(1:3)])
names(cols) <- paste("Development|",order,sep = "")

meta_cols <- substr(meta_cols,13,nchar(meta_cols))
names(cols) <- substr(names(cols),13,nchar(names(cols)))

p2 = plot_density_with_peak_by_disease(
  seurat_obj = TMA_merged_tumor,
  features = meta_cols,
  color_vec = cols,
  disease = "RMS",
  data_type = "assay",
  sample_frac = 0.01,
  assay_name = "AUCell"
)
p1+p2

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

sample_anno <- as.data.frame(table(TMA_merged_tumor$Sample,TMA_merged_tumor$Disease))
sample_anno <- subset(sample_anno,Freq > 0)
sample_anno <- sample_anno[,-3]
colnames(sample_anno) <- c("Samples","Diseases")
rownames(sample_anno) <- sample_anno$Samples

disease_order = c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS","NF","SWN",
                  "LYM","HE","KHE","MPNST","US","AS","RMS","PECOMA","EWS","MRT")
disease_colors = c(
  "#8DD3C7","#97D7C5","#A1DBC3","#ACDFC1","#B6E3BF","#C0E7BD","#CBEBBC","#D5EFBA",
  "#FFFFB3","#F5F5B8","#FCCDE5","#F0D1E1","#E4D5DD","#B3DE69","#80B1D3","#FDB462",
  "#ADD8E6","#BC80BD","#BE8FBE","#C09EBF"
)

TMA_merged_tumor_score_avag <- AverageExpression(TMA_merged_tumor,features = rownames(TMA_merged_tumor@assays$AUCell),assays = "AUCell",group.by = "Sample")
TMA_merged_tumor_score_avag <- TMA_merged_tumor_score_avag$AUCell

p3 = plot_gene_average_boxplot(TMA_merged_tumor_score_avag,
                          sample_anno, 
                          target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                          genes = order[c(3,5,7,4)],
                          disease_order = disease_order,
                          disease_colors = disease_colors,
                          test="t")
p4 = plot_gene_average_boxplot(TMA_merged_tumor_score_avag,
                          sample_anno, 
                          target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                          genes = order[c(16,13,16,10)],
                          disease_order = disease_order,
                          disease_colors = disease_colors,
                          test="t")
p3+p4
aurocs_channel_250210 <- readRDS("/cluster3/yflu/STS/development/aurocs_channel_250210.rds")
aurocs_channel_250210 <- as.data.frame(aurocs_channel_250210)
aurocs_channel_250210 <- aurocs_channel_250210[c(1:78),c(79:96)]
aurocs_channel_250210_1 <- aurocs_channel_250210
rownames(aurocs_channel_250210_1) <- substr(rownames(aurocs_channel_250210_1),5,nchar(rownames(aurocs_channel_250210_1)))
colnames(aurocs_channel_250210_1) <- substr(colnames(aurocs_channel_250210_1),13,nchar(colnames(aurocs_channel_250210_1)))

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

disease_order = c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS","NF","SWN",
                  "LYM","HE","KHE","MPNST","US","AS","RMS","PECOMA","EWS","MRT")
diseases <- substr(diseases,5,nchar(diseases))

anno_sample_cluster_extended$Samples <- rownames(anno_sample_cluster_extended)
colnames(anno_sample_cluster_extended)[1] <- "Diseases"

disease_map <- setNames(disease_order, diseases)

# 替换列中的疾病名
anno_sample_cluster_extended$Diseases <- disease_map[anno_sample_cluster_extended$Diseases]

p1 = plot_gene_average_boxplot(t(aurocs_channel_250210_1),
                               anno_sample_cluster_extended, 
                               target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                               genes = order[c(3,5,7,4)],
                               disease_order = disease_order,
                               disease_colors = disease_colors,
                               test="t")
p2 = plot_gene_average_boxplot(t(aurocs_channel_250210_1),
                               anno_sample_cluster_extended, 
                               target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                               genes = order[c(16,13,16,10)],
                               disease_order = disease_order,
                               disease_colors = disease_colors,
                               test="t")
p1+p2+p3+p4
