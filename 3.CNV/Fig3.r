library(data.table)
library(ggplot2)
library(RColorBrewer)
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
library(maftools)
library(pbapply)
library(S4Vectors)

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5ad")
STS.pega.tumor <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.tumor@meta.data <- metadata$obs

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
#samplenames <- subset(anno_sample_cluster_extended,Malignancy == "Malignant")
#samplenames <- rownames(samplenames)
samplenames <- rownames(anno_sample_cluster_extended)

hmm_regions_path <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat",sep = "")
celltypepath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/celltype/",samplenames,".celltype.integrated.xlsx",sep = "")
outpath <- paste("/cluster3/yflu/STS/WES_CNV/infer_segments_old/",samplenames,"_infer_segments.seg",sep = "")
cnvobjpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/run.final.infercnv_obj",sep = "")

gene_order_file= read.delim("/cluster3/yflu/RT/inferCNV/hg38_gencode_v27.txt",header = F)

genes <- rownames(STS.pega.tumor)
i = 1
cnvobj <- readRDS(cnvobjpath[i])
colnames(cnvobj@expr.data) <- paste(samplenames[i],"-",substr(colnames(cnvobj@expr.data),1,16),sep = "")
cells <- intersect(colnames(STS.pega.tumor),colnames(cnvobj@expr.data))
cnv_mean <- as.data.frame(rowMeans(cnvobj@expr.data[,cells]))
colnames(cnv_mean) <- "mean_cnv"
cnv_mean$gene <- rownames(cnv_mean)

# 对所有目标 genes 建立完整向量，缺失补0
cnv_mean_full <- data.frame(
  gene = genes,
  mean_cnv = cnv_mean$mean_cnv[match(genes, cnv_mean$gene)]
)

# NA 补 0
cnv_mean_full$mean_cnv[is.na(cnv_mean_full$mean_cnv)] <- 1

# 保证行名和 STS.pega.tumor 一致
rownames(cnv_mean_full) <- cnv_mean_full$gene
cnv_mean_full <- cnv_mean_full[ , "mean_cnv", drop = FALSE]
colnames(cnv_mean_full) <- samplenames[i]

for (i in 2:length(samplenames)) {
  cnvobj <- readRDS(cnvobjpath[i])
  colnames(cnvobj@expr.data) <- paste(samplenames[i],"-",substr(colnames(cnvobj@expr.data),1,16),sep = "")
  cells <- intersect(colnames(STS.pega.tumor),colnames(cnvobj@expr.data))
  cnv_mean <- as.data.frame(rowMeans(cnvobj@expr.data[,cells]))
  colnames(cnv_mean) <- "mean_cnv"
  cnv_mean$gene <- rownames(cnv_mean)
  
  # 对所有目标 genes 建立完整向量，缺失补0
  cnv_mean_full_1 <- data.frame(
    gene = genes,
    mean_cnv = cnv_mean$mean_cnv[match(genes, cnv_mean$gene)]
  )
  
  # NA 补 0
  cnv_mean_full_1$mean_cnv[is.na(cnv_mean_full_1$mean_cnv)] <- 1
  
  # 保证行名和 STS.pega.tumor 一致
  rownames(cnv_mean_full_1) <- cnv_mean_full_1$gene
  cnv_mean_full_1 <- cnv_mean_full_1[ , "mean_cnv", drop = FALSE]
  colnames(cnv_mean_full_1) <- samplenames[i]
  cnv_mean_full <- cbind(cnv_mean_full,cnv_mean_full_1)
  print(i)
}

rna_average <- AverageExpression(STS.pega.tumor,group.by = "Channel")
rna_average <- as.data.frame(rna_average$RNA)
rna_average[1:5,1:5]

library(future.apply)

library(pbapply)

calc_cross_gene_cor_parallel <- function(mat1, mat2, method = c("pearson", "spearman"), n_cores = parallel::detectCores() - 1) {
  # 选择相关性方法
  method <- match.arg(method)
  suppressPackageStartupMessages({
    library(parallel)
  })
  
  # 对齐基因和样本
  common_genes   <- intersect(rownames(mat1), rownames(mat2))
  common_samples <- intersect(colnames(mat1), colnames(mat2))
  
  mat1 <- mat1[common_genes, common_samples, drop = FALSE]
  mat2 <- mat2[common_genes, common_samples, drop = FALSE]
  
  genes <- common_genes
  n_genes <- length(genes)
  
  message("Start parallel cross-matrix correlation: ", n_genes, " x ", n_genes, " using ", n_cores, " cores...")
  
  # 初始化结果矩阵
  cor_mat <- matrix(NA, n_genes, n_genes, dimnames = list(genes, genes))
  pval_mat <- matrix(NA, n_genes, n_genes, dimnames = list(genes, genes))
  
  # 分块处理：每个核心计算部分行
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  # 导出变量
  clusterExport(cl, c("mat1", "mat2", "method", "genes"), envir = environment())
  clusterEvalQ(cl, library(stats))
  
  # 定义计算函数
  compute_row <- function(i) {
    x <- as.numeric(mat1[i, ])
    n_genes <- length(genes)
    cor_row <- numeric(n_genes)
    pval_row <- numeric(n_genes)
    for (j in seq_len(n_genes)) {
      y <- as.numeric(mat2[j, ])
      if (sd(x) == 0 || sd(y) == 0) {
        cor_row[j] <- NA
        pval_row[j] <- NA
      } else {
        test <- suppressWarnings(cor.test(x, y, method = method))
        cor_row[j] <- unname(test$estimate)
        pval_row[j] <- test$p.value
      }
    }
    return(list(cor = cor_row, pval = pval_row))
  }
  
  # 并行计算
  results <- parLapply(cl, seq_len(n_genes), compute_row)
  
  # 合并结果
  for (i in seq_len(n_genes)) {
    cor_mat[i, ] <- results[[i]]$cor
    pval_mat[i, ] <- results[[i]]$pval
  }
  
  message("Finished computing cross-gene correlation.")
  return(list(cor = cor_mat, pval = pval_mat))
}

# 用法示例：
res <- calc_cross_gene_cor_parallel(rna_average, cnv_mean_full, method = "pearson", n_cores = 20)
cor_matrix  <- res$cor
pval_matrix <- res$pval

saveRDS(res,"/cluster3/yflu/STS/WES_CNV/cnv_correlations.rds")
###B
cnv_correlations_intersect <- readRDS("/cluster3/yflu/STS/WES_CNV/cnv_correlations_intersect.rds")
ordered_genes <- intersect(gene_order_file$V1,rownames(cnv_correlations_intersect$cor))

cor <- cnv_correlations_intersect$cor[rev(ordered_genes),ordered_genes]
p <- cnv_correlations_intersect$p[rev(ordered_genes),ordered_genes]

cor_1 <- cor
cor_1[p > 0.05] <- 0

row_sums <- rowSums(cor_1)
col_sums <- colSums(cor_1)

# 保留行列和都不为0的基因
keep_genes <- which(row_sums != 0 | col_sums != 0)
cor_1 <- cor_1[keep_genes, keep_genes]

rownames(gene_order_file) <- gene_order_file$V1
gene_order_file_intersect <- gene_order_file[ordered_genes,]

gene_order_file_intersect$chr_num <- as.numeric(sub("chr", "", gene_order_file_intersect$V2))

# 按染色体顺序和原始顺序排序
gene_order_file_intersect <- gene_order_file_intersect[order(gene_order_file_intersect$chr_num, seq_len(nrow(gene_order_file_intersect))), ]

gene_order_file_intersect$group <- unlist(
  lapply(split(seq_len(nrow(gene_order_file_intersect)), gene_order_file_intersect$chr_num), function(idx) {
    ceiling(seq_along(idx) / 50)
  })
)

common_genes <- intersect(rownames(cor_1), rownames(gene_order_file_intersect))
cor_1 <- cor_1[common_genes, common_genes]
gene_info <- gene_order_file_intersect[common_genes, ]

# 生成 group 标签（如 chr1_group1）
gene_info$group_label <- paste0(gene_info$V2, "_group", gene_info$group)

# 唯一的分组标签
groups <- unique(gene_info$group_label)

# 初始化结果矩阵
group_cor <- matrix(NA, nrow = length(groups), ncol = length(groups),
                    dimnames = list(groups, groups))

# 循环计算每组之间的平均值
for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    g1 <- gene_info$group_label == groups[i]
    g2 <- gene_info$group_label == groups[j]
    submat <- cor_1[g1, g2, drop = FALSE]
    group_cor[i, j] <- mean(submat, na.rm = TRUE)
  } 
  print(paste(i,j))
}

# 查看结果
group_cor[1:5, 1:5]

colors = colorRampPalette(brewer.pal(8,'RdBu'))(90)
colors1 = colorRampPalette(brewer.pal(8,'RdBu'))(100)

start_color <- colors1[20]
end_color <- "#FFFFFF"

# 生成渐变函数
grad_fun <- colorRampPalette(c(start_color, end_color))

# 生成10个颜色
colors_white <- grad_fun(3)

start_color_1 <- "#FFFFFF"
end_color_1 <- colors[84]

# 生成渐变函数
grad_fun_1 <- colorRampPalette(c(start_color_1, end_color_1))

# 生成10个颜色
colors_white_1 <- grad_fun_1(3)

color = c(colors1[1:19],colors_white,colors_white_1,colors[84:90])

anno_group <- as.data.frame(rownames(group_cor))
colnames(anno_group) <- "group"
anno_group$chr <- sub("_.*", "", anno_group$group)
rownames(anno_group) <- anno_group$group

anno_group_1 <- as.data.frame(anno_group[,-1])
colnames(anno_group_1) <- "chr"
rownames(anno_group_1) <- anno_group$group

pheatmap::pheatmap(group_cor[,rev(rownames(group_cor))],cluster_rows = F,cluster_cols = F,color = rev(color),
                   show_rownames = F,show_colnames = F,annotation_row = anno_group_1,annotation_col = anno_group_1)

gene_info <- gene_order_file_intersect
common_genes <- intersect(rownames(cor_1), rownames(gene_info))
cor_1 <- cor_1[common_genes, common_genes]
gene_info <- gene_info[common_genes, ]
gene_info$group_label <- paste("chr",gene_info$chr_num,"_group",gene_info$group,sep = "")

# 按 group_label 拆分
groups <- unique(gene_info$group_label)

# 初始化结果表
group_sign_counts <- data.frame(
  group = groups,
  positive_count = NA,
  negative_count = NA
)

# 计算每个 group 内部的正负相关数
for (g in groups) {
  genes_in_group <- gene_info$group_label == g
  submat <- cor_1[genes_in_group, genes_in_group, drop = FALSE]
  
  # 去掉对角线（自相关=1）
  submat[lower.tri(submat, diag = TRUE)] <- NA
  
  group_sign_counts[group_sign_counts$group == g, "positive_count"] <- sum(submat > 0, na.rm = TRUE)
  group_sign_counts[group_sign_counts$group == g, "negative_count"] <- sum(submat < 0, na.rm = TRUE)
}

# 查看结果
group_sign_counts
rownames(group_sign_counts) <- group_sign_counts$group
group_sign_counts <- group_sign_counts[rev(rownames(group_cor)),]

library(ggplot2)
library(dplyr)
library(tidyr)

plot_df <- group_sign_counts %>%
  mutate(negative_count = -negative_count) %>%
  pivot_longer(cols = c("positive_count", "negative_count"),
               names_to = "type", values_to = "count")

ggplot(plot_df, aes(x = group, y = count, fill = type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("positive_count" = "#B2182B", 
                               "negative_count" = "#2166AC"),
                    labels = c("Positive", "Negative")) +
  scale_y_continuous(
    limits = c(-100, NA),                     # y轴范围
    breaks = seq(-100, max(plot_df$count, na.rm = TRUE), by = 200),  # 🔹显示刻度
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Group", y = "Count", fill = "Correlation sign",
       title = "Positive and negative correlation counts per group") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
anno_group <- anno_group[rev(rownames(group_cor)),]
anno_group$index <- seq_len(nrow(anno_group))

# 找出每个染色体最后一个 group 的 index
chr_breaks <- tapply(anno_group$index, anno_group$chr, max)

# 绘图（假设前面已有 ggplot）
plot_df$group <- factor(plot_df$group,levels = rev(rownames(group_cor)))
ggplot(plot_df[,], aes(x = group, y = count, fill = type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("positive_count" = "#E64B35", 
                               "negative_count" = "#4DBBD5"),
                    labels = c("Negative", "Positive")) +
  scale_y_continuous(
    limits = c(-100, NA),
    breaks = seq(-100, max(plot_df$count, na.rm = TRUE), by = 200)
  ) +
  geom_vline(xintercept = chr_breaks + 0.5, color = "gray20", linetype = "dashed", linewidth = 0.3) +  # 🔹加分割线
  labs(x = "Group", y = "Count", fill = "Correlation sign") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

gene_order_file
chr_len <- as.data.frame(cbind(c(1:22),c(1:22)))
i=1
len <- subset(gene_order_file,V2 == paste("chr",i,sep = ""))
len <- max(len$V4) - min(len$V3)
chr_len$V2[i] <- len
for (i in 2:22) {
  len <- subset(gene_order_file,V2 == paste("chr",i,sep = ""))
  len <- max(len$V4) - min(len$V3)
  chr_len$V2[i] <- len
}

seg_merged <- read.delim("/cluster3/yflu/STS/WES_CNV/infer_segments/segments_merged.seg")

i=1
seg <- subset(seg_merged,Sample == samplenames[i])
chr_statistic <- as.data.frame(cbind(c(1:22),c(1:22)))
colnames(chr_statistic) <- c(paste(samplenames[i],"amp",sep = ""),paste(samplenames[i],"del",sep = ""))
for (j in 1:22) {
  chr <- paste("chr",j,sep = "")
  seg_sub <- subset(seg,Chromosome == chr)
  if (max(seg_sub$Segment_Mean > 0)) {
    seg_amp <- subset(seg_sub,Segment_Mean > 0)
    seg_amp$len <- seg_amp$End-seg_amp$Start
    chr_statistic[j,1] <- sum(seg_amp$len)/chr_len$V2[j]
  }else {
    chr_statistic[j,1] <- 0
  }
  if (max(seg_sub$Segment_Mean < 0)) {
    seg_amp <- subset(seg_sub,Segment_Mean < 0)
    seg_amp$len <- seg_amp$End-seg_amp$Start
    chr_statistic[j,2] <- sum(seg_amp$len)/chr_len$V2[j]
  }else {
    chr_statistic[j,2] <- 0
  }
}

for (i in 2:length(samplenames)) {
  seg <- subset(seg_merged,Sample == samplenames[i])
  chr_statistic_1 <- as.data.frame(cbind(c(1:22),c(1:22)))
  colnames(chr_statistic_1) <- c(paste(samplenames[i],"amp",sep = ""),paste(samplenames[i],"del",sep = ""))
  for (j in 1:22) {
    chr <- paste("chr",j,sep = "")
    seg_sub <- subset(seg,Chromosome == chr)
    if (max(seg_sub$Segment_Mean > 0)) {
      seg_amp <- subset(seg_sub,Segment_Mean > 0)
      seg_amp$len <- seg_amp$End-seg_amp$Start
      chr_statistic_1[j,1] <- sum(seg_amp$len)/chr_len$V2[j]
    }else {
      chr_statistic_1[j,1] <- 0
    }
    if (max(seg_sub$Segment_Mean < 0)) {
      seg_amp <- subset(seg_sub,Segment_Mean < 0)
      seg_amp$len <- seg_amp$End-seg_amp$Start
      chr_statistic_1[j,2] <- sum(seg_amp$len)/chr_len$V2[j]
    }else {
      chr_statistic_1[j,2] <- 0
    }
  }
  chr_statistic <- cbind(chr_statistic,chr_statistic_1)
}
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
disease_counts <- as.data.frame(table(anno_sample_cluster_extended$Disease))
anno_sub <- subset(anno_sample_cluster_extended,Disease %in% as.character(subset(disease_counts,Freq > 2)$Var1))

chr_statistic_amp <- chr_statistic[,2*(1:length(samplenames))-1]
chr_statistic_del <- chr_statistic[,2*(1:length(samplenames))]
colnames(chr_statistic_amp) <- samplenames
colnames(chr_statistic_del) <- samplenames

chr_statistic_amp <- chr_statistic_amp[,rownames(anno_sub)]
chr_statistic_del <- chr_statistic_del[,rownames(anno_sub)]

calc_amp_proportion_by_disease <- function(chr_statistic_amp, anno_sub) {
  # 确保样本名匹配
  common_samples <- intersect(colnames(chr_statistic_amp), rownames(anno_sub))
  if(length(common_samples) == 0) stop("chr_statistic_amp 和 anno_sub 没有共同样本名！")
  
  # 提取匹配部分
  chr_statistic_amp_sub <- chr_statistic_amp[, common_samples, drop = FALSE]
  anno_sub_match <- anno_sub[common_samples, , drop = FALSE]
  
  # 按 Disease 分组
  disease_groups <- split(colnames(chr_statistic_amp_sub), anno_sub_match$Disease)
  
  # 计算每行在每个疾病中的阳性比例
  prop_matrix <- sapply(disease_groups, function(samples) {
    rowMeans(chr_statistic_amp_sub[, samples, drop = FALSE])
  })
  
  # 保持原行名
  rownames(prop_matrix) <- rownames(chr_statistic_amp)
  
  return(prop_matrix)
}

# 使用示例
amp_prop_by_disease <- calc_amp_proportion_by_disease(chr_statistic_amp, anno_sub)
del_prop_by_disease <- calc_amp_proportion_by_disease(chr_statistic_del, anno_sub)

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))
common_disease <- intersect(labels,colnames(del_prop_by_disease))
amp_prop_by_disease <- amp_prop_by_disease[,common_disease[c(2:6,1:10)]]
del_prop_by_disease <- del_prop_by_disease[,common_disease[c(2:6,1:10)]]

###A
library(ggplot2)
library(reshape2)

plot_amp_del_by_disease_shaded <- function(amp_prop, del_prop, top_n = 20, disease_colors = NULL) {
  # 取前 top_n 行
  n_plot <- min(top_n, nrow(amp_prop))
  amp_sub <- amp_prop[1:n_plot, , drop = FALSE]
  del_sub <- del_prop[1:n_plot, , drop = FALSE]
  
  # 转长格式
  df_amp <- melt(amp_sub)
  colnames(df_amp) <- c("Chr", "Disease", "Proportion")
  df_amp$Type <- "Amp"
  
  df_del <- melt(del_sub)
  colnames(df_del) <- c("Chr", "Disease", "Proportion")
  df_del$Proportion <- -df_del$Proportion
  df_del$Type <- "Del"
  
  df_all <- rbind(df_amp, df_del)
  df_all$Chr <- factor(df_all$Chr, levels = rownames(amp_sub))
  
  # 背景色数据（白灰交替）
  chr_levels <- levels(df_all$Chr)
  bg_df <- data.frame(
    xmin = as.numeric(chr_levels) - 0.5,
    xmax = as.numeric(chr_levels) + 0.5,
    ymin = min(df_all$Proportion),
    ymax = max(df_all$Proportion),
    fill_bg = rep(c("#FFFFFF", "#F0F0F0"), length.out = length(chr_levels))
  )
  
  # 绘图
  p <- ggplot() +
    # 背景色
    geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = bg_df$fill_bg, inherit.aes = FALSE, show.legend = FALSE) +
    # 柱状图
    geom_bar(data = df_all, aes(x = Chr, y = Proportion, fill = Disease),
             stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Chr", y = "Proportion", title = "Amp/Del proportion by Disease") +
    geom_hline(yintercept = 0, color = "black")
  
  # 如果提供了 Disease 颜色，使用自定义颜色
  if(!is.null(disease_colors)) {
    # 确保 Disease 因子顺序与颜色长度匹配
    disease_levels <- levels(df_all$Disease)
    if(length(disease_colors) < length(disease_levels)) {
      warning("disease_colors长度小于Disease水平数，将自动重复填充")
      disease_colors <- rep(disease_colors, length.out = length(disease_levels))
    }
    names(disease_colors) <- disease_levels
    p <- p + scale_fill_manual(values = disease_colors)
  }
  
  return(p)
}

# 使用示例
library(RColorBrewer)

# 自定义 Disease 颜色
disease_colors <- c(
  colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(6,7)],
  colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1)],
  colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(2,3)],
  colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(5)],
  colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
  "#ADD8E6",
  colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(2:3)]
)

p <- plot_amp_del_by_disease_shaded(amp_prop_by_disease, del_prop_by_disease, top_n = 22, disease_colors = disease_colors)
print(p)

###D
library(S4Vectors)
infer.gistic <- readRDS("infer.gistic.rds")

recurrent_cnv <- infer.gistic@cytoband.summary
recurrent_cnv <- separate(
  data = recurrent_cnv,
  col = Wide_Peak_Limits,        # 待分割的列
  into = c("Chromosome", "Start_End"),  # 分割后的列名（先拆分成两列）
  sep = ":",                     # 按冒号首次分割
  remove = FALSE                 # 保留原始列
)

# 第二次分割（拆分起止位置）
recurrent_cnv <- separate(
  data = recurrent_cnv,
  col = Start_End,               # 对中间列继续分割
  into = c("Start", "End"),      # 最终需要的两列
  sep = "-",                     # 按短横线分割
  convert = TRUE                 # 自动转换为数值型
)
recurrent_cnv$Cytoband <- paste(recurrent_cnv$Cytoband,recurrent_cnv$Variant_Classification,sep = "_")
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

gene_list <- List()
gene_data <- infer.gistic@data
for (i in 1:nrow(recurrent_cnv)) {
  gene_list <- append(gene_list,list(unique(subset(gene_data,Cytoband == recurrent_cnv$Unique_Name[i])$Hugo_Symbol)))
  names(gene_list)[i+1] <- recurrent_cnv$Cytoband[i]
}
gene_list[[2]] <- gene_list[[1]]
gene_list <- gene_list[-1]

STS.integrated.pega <- readRDS("/cluster3/yflu/STS/pegasus/STS.integrated.pega_20240507.rds")
STS.integrated.pega[["pca_regressed"]] <- NULL
STS.integrated.pega[["pca_regressed_harmony"]] <- NULL
STS.integrated.pega[["pca_regressed_harmony_knn_distances"]] <- NULL
STS.integrated.pega[["pca_regressed_harmony_knn_indices"]] <- NULL

STS.integrated.pega_score <- irGSEA.score(object = STS.integrated.pega, assay = "RNA", slot = "data", 
                                          seeds = 123, ncores = 1,msigdb=F, 
                                          custom = T, geneset = gene_list, method = c("AUCell", "UCell", "singscore",
                                                                                      "ssgsea")[c(1)], 
                                          kcdf = 'Gaussian')
STS_score <- AverageExpression(STS.integrated.pega_score,assays = "AUCell",group.by = "Channel")
STS_score <- as.data.frame(STS_score$AUCell)
STS_score <- as.data.frame(t(STS_score))

colnames(STS_score) <- gsub("-", "_", colnames(STS_score), fixed = TRUE)
saveRDS(STS_score,"STS_score.rds")

library(ggplot2)
library(RColorBrewer)
library(ggsignif)

library(ggplot2)
library(RColorBrewer)
library(ggsignif)

plot_cnv_box_by_malignancy <- function(STS_score_sub, plot_title = "CNV score distribution by Malignancy") {
  df <- STS_score_sub
  df$Malignancy <- factor(df$Malignancy, levels = c("Benign", "Malignant"))
  
  if (!is.factor(df$Disease)) {
    df$Disease <- factor(df$Disease)
  }
  
  # ---------------------
  # 颜色设置
  # ---------------------
  box_colors <- c(
    "Benign" = "#ADD8E6",   # 浅蓝
    "Malignant" = "#FFC0CB" # 浅粉
  )
  
  disease_levels <- levels(df$Disease)
  disease_palette <- c(
    colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
    colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
    colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
    colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
    colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
    colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
    "#ADD8E6",
    colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)]
  )
  disease_colors <- setNames(disease_palette[seq_along(disease_levels)], disease_levels)
  
  # ---------------------
  # Wilcoxon 检验
  # ---------------------
  wilcox_res <- wilcox.test(CNV_score ~ Malignancy, data = df)
  p_value <- wilcox_res$p.value
  # 星号表示法
  p_label <- if (p_value < 0.001) {
    "***"
  } else if (p_value < 0.01) {
    "**"
  } else if (p_value < 0.05) {
    "*"
  } else {
    "ns"
  }
  
  # ---------------------
  # 绘图
  # ---------------------
  p <- ggplot(df, aes(x = Malignancy, y = CNV_score, fill = Malignancy)) +
    geom_boxplot(outlier.shape = NA, color = "gray40", alpha = 0.7) +
    geom_jitter(aes(color = Disease),
                width = 0.15, size = 2.5, alpha = 0.9) +
    geom_signif(
      comparisons = list(c("Benign", "Malignant")),
      annotations = p_label,
      y_position = max(df$CNV_score) * 1.05,
      tip_length = 0.03,
      textsize = 5
    ) +
    scale_fill_manual(values = box_colors) +
    scale_color_manual(values = disease_colors) +
    theme_bw() +
    labs(
      x = "Malignancy",
      y = "CNV score",
      title = plot_title
    ) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}
i=1
STS_score_sub <- as.data.frame(STS_score[,i])
rownames(STS_score_sub) <- rownames(STS_score)
STS_score_sub <- STS_score_sub[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),]
STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),])
STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)
STS_score_sub$Malignancy <- ifelse(STS_score_sub$Disease == "RMS", "Malignant", "Benign")
rownames(STS_score_sub) <- rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant"))
colnames(STS_score_sub)[1] <- "CNV_score"

p <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])
for (i in 2:ncol(STS_score)) {
  STS_score_sub <- as.data.frame(STS_score[,i])
  rownames(STS_score_sub) <- rownames(STS_score)
  STS_score_sub <- STS_score_sub[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),]
  STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),])
  STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)
  STS_score_sub$Malignancy <- ifelse(STS_score_sub$Disease == "RMS", "Malignant", "Benign")
  rownames(STS_score_sub) <- rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant"))
  colnames(STS_score_sub)[1] <- "CNV_score"
  
  p_1 <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])
  p = p + p_1
}

p

calc_p_by_malignancy_matrix <- function(STS_score, anno_sample_cluster_extended, group_col = "Malignancy") {
  # STS_score: 样本 × 特征的数据框或矩阵，行名为样本名
  # anno_sample_cluster_extended: 样本信息数据框，行名为样本名，必须包含 group_col
  # group_col: 分组列名，例如 "Malignancy"
  
  # 确保样本顺序一致
  common_samples <- intersect(rownames(STS_score), rownames(anno_sample_cluster_extended))
  STS_score <- STS_score[common_samples, , drop = FALSE]
  anno <- anno_sample_cluster_extended[common_samples, , drop = FALSE]
  
  res <- data.frame(
    Variable = colnames(STS_score),
    p_value = NA,
    significance = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(colnames(STS_score))) {
    col <- colnames(STS_score)[i]
    tmp_df <- data.frame(
      value = STS_score[, col],
      group = factor(anno[[group_col]])
    )
    
    # 只在有两组时进行 Wilcoxon 检验
    if (length(unique(tmp_df$group)) == 2) {
      test <- wilcox.test(value ~ group, data = tmp_df)
      p <- test$p.value
    } else {
      p <- NA
    }
    
    sig <- if (is.na(p)) {
      NA
    } else if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      "ns"
    }
    
    res$p_value[i] <- p
    res$significance[i] <- sig
  }
  
  return(res)
}

# ---------------------
# 使用示例
# ---------------------
p_table <- calc_p_by_malignancy_matrix(STS_score, anno_sample_cluster_extended)
head(p_table)

write.xlsx(p_table,"p_table.xlsx")

i=1
STS_score_sub <- as.data.frame(STS_score[,i])
rownames(STS_score_sub) <- rownames(STS_score)
colnames(STS_score_sub)[1] <- "CNV_score"
STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(STS_score),])
STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)
STS_score_sub <- subset(STS_score_sub,Malignancy == "Malignant")

p <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])
for (i in 2:ncol(STS_score)) {
  STS_score_sub <- as.data.frame(STS_score[,i])
  rownames(STS_score_sub) <- rownames(STS_score)
  colnames(STS_score_sub)[1] <- "CNV_score"
  STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(STS_score),])
  STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)
  
  p_1 <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])
  p = p + p_1
}

STS_score_RMS <- STS_score[rownames(subset(anno_sample_cluster_extended,Disease == "RMS")),]
STS_score_EWS <- STS_score[rownames(subset(anno_sample_cluster_extended,Disease == "EWS/PNET")),]

STS_score_RMS_sd <- sapply(STS_score_RMS, sd, na.rm = TRUE)
STS_score_EWS_sd <- sapply(STS_score_EWS, sd, na.rm = TRUE)
sd <- as.data.frame(cbind(STS_score_RMS_sd,STS_score_EWS_sd))

p_table_anno <- read.xlsx("p_table_anno.xlsx","Sheet1")
rownames(p_table_anno) <- p_table_anno$Variable
p_table_anno[order(p_table_anno$Variable),]

vars <- p_table_anno$Variable

# 提取染色体号、臂、位置数字
library(dplyr)
library(stringr)

sorted_vars <- data.frame(
  var = vars,
  chr = as.numeric(str_extract(vars, "^[0-9]+")),
  arm = str_extract(vars, "(?<=\\d)(p|q)"),
  pos = as.numeric(str_extract(vars, "(?<=p|q)[0-9]+\\.?[0-9]*"))
) %>%
  mutate(
    # 为排序创建辅助列
    arm_order = ifelse(arm == "p", 1, 2),
    # p臂数字大者靠前（降序），q臂数字小者靠前（升序）
    pos_order = ifelse(arm == "p", -pos, pos)
  ) %>%
  arrange(chr, arm_order, pos_order) %>%
  pull(var)

# 更新排序结果
p_table_anno$Variable <- factor(p_table_anno$Variable, levels = sorted_vars)

# 查看排序结果
sorted_vars
p_table_anno <- p_table_anno[sorted_vars,]

library(ggplot2)
library(reshape2)
library(dplyr)

plot_p_table_heatmap <- function(p_table_anno) {
  # 选取后五列
  df <- p_table_anno[, c("Variable", "malignancy", "RMS", "EWS", "ITCC_RMS", "COG_EWS")]
  
  # 转换为长格式
  df_long <- melt(df, id.vars = "Variable", variable.name = "Category", value.name = "Significance")
  
  # 判断Amp或Del
  df_long$Type <- ifelse(grepl("Amp", df_long$Variable), "Amp",
                         ifelse(grepl("Del", df_long$Variable), "Del", NA))
  
  # 定义颜色映射
  color_map <- c(
    # Amp（红色系）
    "***_Amp" = "#8B0000",   # 深红
    "**_Amp"  = "#CD5C5C",   # 中红
    "*_Amp"   = "#F4A6A6",   # 浅红
    
    # Del（蓝色系）
    "***_Del" = "#08306B",   # 深蓝
    "**_Del"  = "#4292C6",   # 中蓝
    "*_Del"   = "#9ECAE1",   # 浅蓝
    
    # 其他
    "ns"       = "grey95",
    "reversed" = "grey60"
  )
  
  # 生成颜色键
  df_long$Significance_Type <- ifelse(df_long$Significance %in% c("***", "**", "*"),
                                      paste0(df_long$Significance, "_", df_long$Type),
                                      df_long$Significance)
  
  # 只在非 ns 和 reversed 时显示文字
  df_long$Label <- ifelse(df_long$Significance %in% c("ns", "reversed"), "", df_long$Significance)
  
  # 根据背景深浅确定文字颜色
  df_long$Label_Color <- sapply(df_long$Significance_Type, function(x) {
    if (x %in% c("***_Amp", "**_Amp", "***_Del", "**_Del", "reversed")) {
      "white"
    } else {
      "black"
    }
  })
  
  # 调整 y 轴顺序（反转）
  df_long$Category <- factor(df_long$Category,
                             levels = rev(c("malignancy", "RMS", "EWS", "ITCC_RMS", "COG_EWS")))
  
  # 绘图
  p <- ggplot(df_long, aes(x = Variable, y = Category, fill = Significance_Type)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Label, color = Label_Color), size = 3) +
    scale_color_identity() +
    scale_fill_manual(values = color_map, na.value = "white") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

# 示例：
p <- plot_p_table_heatmap(p_table_anno)
print(p)

i=match("2q33.1_Amp",colnames(STS_score))
STS_score_sub <- as.data.frame(STS_score[,i])
rownames(STS_score_sub) <- rownames(STS_score)
colnames(STS_score_sub)[1] <- "CNV_score"
STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(STS_score),])
STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)

p1 <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])

STS_score_sub <- as.data.frame(STS_score[,i])
rownames(STS_score_sub) <- rownames(STS_score)
STS_score_sub <- STS_score_sub[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),]
STS_score_sub <- cbind(STS_score_sub,anno_sample_cluster_extended[rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant")),])
STS_score_sub$Disease <- factor(STS_score_sub$Disease,levels = labels)
STS_score_sub$Malignancy <- ifelse(STS_score_sub$Disease == "RMS", "Benign", "Malignant")
rownames(STS_score_sub) <- rownames(subset(anno_sample_cluster_extended,Malignancy == "Malignant"))
colnames(STS_score_sub)[1] <- "CNV_score"

p2 <- plot_cnv_box_by_malignancy(STS_score_sub, plot_title = colnames(STS_score)[i])
p1
p2
gene_list_cnv <- gene_list
gene_list_cnv[["2p25.1_Amp"]]
gene_list_cnv[["2q33.1_Amp"]]
gene_list_cnv[["21q22.13_Amp"]]
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_list <- gene_list_cnv[["21q22.13_Amp"]]

annot <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = gene_list,
  mart = ensembl
)
table(annot$gene_biotype)

library(gprofiler2)
library(PANTHER.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_list <- subset(annot,gene_biotype == "protein_coding")$external_gene_name

# SYMBOL -> ENTREZID
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_list,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# ENTREZID -> PANTHER CLASS_TERM
protein_class_raw <- mapIds(PANTHER.db,
                            keys = entrez_ids,
                            column = "CLASS_TERM",
                            keytype = "ENTREZ",
                            multiVals = function(x) paste(x, collapse=";"))  # 多值合并

# 构建表格，对齐 gene_list
protein_class <- sapply(entrez_ids, function(x) {
  if (!is.na(x) && x %in% names(protein_class_raw)) {
    protein_class_raw[x]
  } else {
    NA
  }
})

df <- data.frame(
  Gene = gene_list,
  EntrezID = entrez_ids,
  ProteinClass = protein_class,
  stringsAsFactors = FALSE
)

df

cnv_genes_list <- list(annot,df)
saveRDS(cnv_genes_list,"21q22.13_Amp_cnv_genes_list.rds")
write.xlsx(df,"21q22.13_Amp_cnv_genes.xlsx")

library(ggplot2)
library(dplyr)
library(patchwork)

plot_biotype_group_pies <- function(annot, annot_protein,
                                    biotype_colors, group_colors,
                                    marker_genes = NULL) {
  # --- 第1个饼图：gene_biotype 分布 ---
  p1_data <- annot %>%
    count(gene_biotype)
  
  p1 <- ggplot(p1_data, aes(x = "", y = n, fill = gene_biotype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = "Gene biotype distribution") +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_fill_manual(values = biotype_colors)
  
  # --- 第2个饼图：protein_coding 基因 Group 分布 ---
  p2_data <- annot_protein %>%
    count(Group)
  
  if (!is.null(marker_genes)) {
    # 按 marker_genes 全局顺序拼接，每个基因换行
    group_labels <- annot_protein %>%
      filter(Gene %in% marker_genes) %>%
      mutate(Gene = factor(Gene, levels = marker_genes)) %>%
      arrange(Gene) %>%
      group_by(Group) %>%
      summarise(markers_in_group = paste(Gene, collapse = "\n"), .groups = "drop")
    
    p2_data <- p2_data %>%
      left_join(group_labels, by = "Group") %>%
      mutate(markers_in_group = ifelse(is.na(markers_in_group), "", markers_in_group))
  } else {
    p2_data$markers_in_group <- ""
  }
  
  p2 <- ggplot(p2_data, aes(x = "", y = n, fill = Group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = "Group distribution (protein-coding genes)") +
    geom_text(aes(label = markers_in_group),
              position = position_stack(vjust = 0.5), size = 3,
              lineheight = 0.9, hjust = 0.5) +  # 控制行高和居中
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_fill_manual(values = group_colors)
  
  # --- 拼图返回 ---
  return(p1 | p2)
}

colors_4 <- c(
  "#BEBADA",  # 淡紫色
  "#8DD3C7",  # 青绿色
  "#FDB462",  # 橙黄色
  "#80B1D3"   # 天蓝色
)
colors_10 <- c(
  "#FB8072", # 红
  "#FDB462", # 橙
  "#FFED6F", # 黄
  "#B3DE69", # 亮绿
  "#8DD3C7", # 青绿
  "#80B1D3", # 蓝
  "#BEBADA", # 淡紫
  "#BC80BD", # 紫
  "#FCCDE5", # 粉
  "#D9D9D9"  # 灰白中性
)

cnv_marker_list_all <- readRDS("/cluster3/yflu/STS/WES_CNV/cnv_marker_list_all.rds")

cnv_2p25_1_Amp_cnv_genes_list <- readRDS("/cluster3/yflu/STS/WES_CNV/2p25.1_Amp_cnv_genes_list.rds")
cnv_2p25_1_Amp_cnv_genes_list[[2]] <- read.xlsx("/cluster3/yflu/STS/WES_CNV/2p25.1_Amp_cnv_genes.xlsx","Sheet1")

cnv_marker_list_all <- readRDS("/cluster3/yflu/STS/WES_CNV/cnv_marker_list_all.rds")
cnv_2p25_1_markers <- cnv_marker_list_all[["2p25_1_Amp"]]
cnv_2p25_1_markers <- subset(cnv_2p25_1_markers,gene %in% gene_list_cnv[["2p25.1_Amp"]]&cluster == 1&avg_log2FC>0)

cnv_2p25_1_markers_top <- cnv_2p25_1_markers$gene[1:10]

cnv_2p25_1_Amp_cnv_genes_list[[1]]$gene_biotype <- factor(cnv_2p25_1_Amp_cnv_genes_list[[1]]$gene_biotype,
                                                          levels = rev(c("protein_coding","lncRNA","miRNA","snoRNA")))
cnv_2p25_1_Amp_cnv_genes_list[[2]]$Group <- factor(cnv_2p25_1_Amp_cnv_genes_list[[2]]$Group,
                                                   levels = rev(unique(cnv_2p25_1_Amp_cnv_genes_list[[2]]$Group)[order(unique(cnv_2p25_1_Amp_cnv_genes_list[[2]]$Group))]))

plot_biotype_group_pies(annot = cnv_2p25_1_Amp_cnv_genes_list[[1]],
                        annot_protein = cnv_2p25_1_Amp_cnv_genes_list[[2]],
                        biotype_colors = colors_4,
                        group_colors = rev(colors_10),
                        marker_genes = cnv_2p25_1_markers_top)

cnv_2q33_1_Amp_cnv_genes_list <- readRDS("/cluster3/yflu/STS/WES_CNV/2q33.1_Amp_cnv_genes_list.rds")
cnv_2q33_1_Amp_cnv_genes_list[[2]] <- read.xlsx("/cluster3/yflu/STS/WES_CNV/2q33.1_Amp_cnv_genes.xlsx","Sheet1")

cnv_2q33_1_markers <- cnv_marker_list_all[["2q33_1_Amp"]]
cnv_2q33_1_markers <- subset(cnv_2q33_1_markers,gene %in% gene_list_cnv[["2q33.1_Amp"]]&cluster == 1&avg_log2FC>0)

cnv_2q33_1_markers_top <- cnv_2q33_1_markers$gene[1:10]

cnv_2q33_1_Amp_cnv_genes_list[[1]]$gene_biotype <- factor(cnv_2q33_1_Amp_cnv_genes_list[[1]]$gene_biotype,
                                                          levels = rev(c("protein_coding","lncRNA","miRNA","snoRNA")))
cnv_2q33_1_Amp_cnv_genes_list[[2]]$Group <- factor(cnv_2q33_1_Amp_cnv_genes_list[[2]]$Group,
                                                   levels = rev(unique(cnv_2q33_1_Amp_cnv_genes_list[[2]]$Group)[order(unique(cnv_2q33_1_Amp_cnv_genes_list[[2]]$Group))]))

plot_biotype_group_pies(annot = cnv_2q33_1_Amp_cnv_genes_list[[1]],
                        annot_protein = cnv_2q33_1_Amp_cnv_genes_list[[2]],
                        biotype_colors = colors_4[-c(1,2)],
                        group_colors = rev(colors_10[-c(3,5)]),
                        marker_genes = cnv_2q33_1_markers_top)

cnv_21q22.13_Amp_cnv_genes_list <- readRDS("/cluster3/yflu/STS/WES_CNV/21q22.13_Amp_cnv_genes_list.rds")
cnv_21q22.13_Amp_cnv_genes_list[[2]] <- read.xlsx("/cluster3/yflu/STS/WES_CNV/21q22.13_Amp_cnv_genes.xlsx","Sheet1")

cnv_21q22.13_markers <- cnv_marker_list_all[["21q22_13_Amp"]]
cnv_21q22.13_markers <- subset(cnv_21q22.13_markers,gene %in% cnv_21q22.13_Amp_cnv_genes_list[[2]]$Gene&cluster == 1&avg_log2FC>0)

cnv_21q22.13_markers_top <- cnv_21q22.13_markers$gene[1:10]

cnv_21q22.13_Amp_cnv_genes_list[[1]]$gene_biotype <- factor(cnv_21q22.13_Amp_cnv_genes_list[[1]]$gene_biotype,
                                                            levels = rev(c("protein_coding","lncRNA","miRNA","snoRNA","transcribed_processed_pseudogene")))
cnv_21q22.13_Amp_cnv_genes_list[[2]]$Group <- factor(cnv_21q22.13_Amp_cnv_genes_list[[2]]$Group,
                                                     levels = rev(unique(cnv_21q22.13_Amp_cnv_genes_list[[2]]$Group)[order(unique(cnv_21q22.13_Amp_cnv_genes_list[[2]]$Group))]))

plot_biotype_group_pies(annot = cnv_21q22.13_Amp_cnv_genes_list[[1]],
                        annot_protein = cnv_21q22.13_Amp_cnv_genes_list[[2]],
                        biotype_colors = c("darkgray",colors_4),
                        group_colors = rev(c(colors_10[1],"orange",colors_10[-c(1,3,5)])),
                        marker_genes = cnv_21q22.13_markers_top)

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)], colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)], colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)], colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)], colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)], colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)], "#ADD8E6", colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 假设 cnv_2p25_1_markers_top 是你要画的基因列表
# 获取表达矩阵（基因 x 细胞）
expr_mat <- as.data.frame(as.matrix(GetAssayData(STS.integrated.pega, slot = "data")))
expr_mat$Gene <- rownames(expr_mat)

# 只保留目标基因
expr_mat_1 <- expr_mat %>% filter(Gene %in% cnv_2p25_1_markers_top[1:5])

# 转为长格式：Gene x Cell x Expression
expr_long <- expr_mat_1 %>%
  pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression")

# 加上细胞对应 Disease 信息
meta <- STS.integrated.pega@meta.data

STS_meta <- readRDS("/cluster3/yflu/STS/WES_CNV/STS_meta.rds")

expr_long <- expr_long %>%
  left_join(STS_meta %>% rownames_to_column("Cell") %>% dplyr::select(Cell, "2p25_1_Amp"), by = "Cell")

expr_long$Gene <- factor(expr_long$Gene, levels = cnv_2p25_1_markers_top[1:5])
expr_long$Disease <- factor(expr_long$"2p25_1_Amp")
#expr_long$Disease <- factor(expr_long$Disease, levels = labels)


# 绘制 violin plot
ggplot(expr_long, aes(x = Disease, y = Expression, fill = Disease)) +
  geom_violin(trim = TRUE, scale = "width") +
  facet_wrap(~Gene, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = c("#FDB462","#B3DE69")) +
  labs(x = "Disease", y = "Expression", fill = "Disease") +
  theme(
    panel.background = element_blank(),       # 背景去掉
    panel.grid = element_blank(),             # 网格线去掉
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(0.2, "cm")
  )

cnvscorepath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.xlsx",sep = "")
i=1
cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
cnv_mean <- cnv_score %>%
  group_by(group2) %>%
  summarise(mean_cnvscore = mean(cnvscore, na.rm = TRUE))
cnv_score_sample <- cnv_mean$mean_cnvscore[2] - cnv_mean$mean_cnvscore[1]
for (i in 2:length(samplenames)) {
  cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
  cnv_mean <- cnv_score %>%
    group_by(group2) %>%
    summarise(mean_cnvscore = mean(cnvscore, na.rm = TRUE))
  cnv_score_sample_1 <- cnv_mean$mean_cnvscore[2] - cnv_mean$mean_cnvscore[1]
  cnv_score_sample <- rbind(cnv_score_sample,cnv_score_sample_1)
  print(i)
}
cnv_score_sample <- as.data.frame(cnv_score_sample)
rownames(cnv_score_sample) <- samplenames

cnv_score_sample <- cbind(cnv_score_sample,anno_sample_cluster_extended)

ggboxplot(
  cnv_score_sample,
  x = "Malignancy",
  y = "V1",
  fill = "Malignancy",
  palette = "Set2",
  add = "jitter",              # 添加散点，展示所有样本点
  shape = 21,                  # 圆点
  size = 2,                    # 点大小
  alpha = 0.7                  # 透明度
) +
  stat_compare_means(
    method = "wilcox.test",    # 非参数检验
    label = "p.signif"         # 显著性星号（*, **, ***）
  ) +
  theme_bw() +
  labs(
    x = "Malignancy",
    y = "CNV score (V1)",
    title = "CNV score distribution by malignancy"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )

cnv_score_sample$Disease <- factor(cnv_score_sample$Disease,levels = labels)

colnames(cnv_score_sample)[1] <- "CNV_score"
plot_cnv_box_by_malignancy(cnv_score_sample)

top_genes <- rbind(cnv_2p25_1_markers_top,cnv_2q33_1_markers_top,cnv_21q22.13_markers_top)

top_genes_avg_exp <- AverageExpression(STS.pega.tumor,features = c(top_genes[1,],top_genes[2,],top_genes[3,]),group.by = "Channel")
top_genes_avg_exp <- as.data.frame(top_genes_avg_exp$RNA)

cnv_anno <- as.data.frame(c(top_genes[1,],top_genes[2,],top_genes[3,]))
cnv_anno$cnv <- c(rep("2p25_1_amp",10),rep("2q33_1_amp",10),rep("21q22_13_amp",10))
cnv_anno <- as.data.frame(cnv_anno[,-1])
rownames(cnv_anno) <- c(top_genes[1,],top_genes[2,],top_genes[3,])
colnames(cnv_anno) <- "cnv"

sample_cnv <- as.data.frame(table(STS_meta$Channel,STS_meta$"2p25_1_Amp"))
sample_cnv <- subset(sample_cnv,Freq >0)
rownames(sample_cnv) <- sample_cnv$Var1
sample_cnv <- sample_cnv[,-3]
colnames(sample_cnv) <- c("sample","2p25_1_Amp")
sample_cnv_1 <- as.data.frame(table(STS_meta$Channel,STS_meta$"2q33_1_Amp"))
sample_cnv_1 <- subset(sample_cnv_1,Freq >0)
rownames(sample_cnv_1) <- sample_cnv_1$Var1
sample_cnv_1 <- sample_cnv_1[rownames(sample_cnv),]
sample_cnv <- cbind(sample_cnv,sample_cnv_1[,-c(1,3)])
colnames(sample_cnv)[3] <- c("2q33_1_Amp")

sample_cnv_1 <- as.data.frame(table(STS_meta$Channel,STS_meta$"21q22_13_Amp"))
sample_cnv_1 <- subset(sample_cnv_1,Freq >0)
rownames(sample_cnv_1) <- sample_cnv_1$Var1
sample_cnv_1 <- sample_cnv_1[rownames(sample_cnv),]
sample_cnv <- cbind(sample_cnv,sample_cnv_1[,-c(1,3)])
colnames(sample_cnv)[4] <- c("21q22_13_Amp")

sample_anno <- cbind(anno_sample_cluster_extended[,],sample_cnv[rownames(anno_sample_cluster_extended),-1])
saveRDS(sample_anno,"/cluster3/yflu/STS/WES_CNV/sample_anno_CNV.rds")

disease_colors <- c(
  colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[1:8],
  colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[1:2],
  colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[1:3],
  colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[1],
  colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[1],
  colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[1],
  "#ADD8E6",
  colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[1:3]
)
# Disease 名称需与你的 sample_anno$Disease levels 顺序一致
names(disease_colors) <- labels
# 2️⃣ Malignancy 颜色
malignancy_colors <- c(
  "Benign" = "#ADD8E6",     # 浅蓝
  "Malignant" = "#FFC0CB"   # 浅粉
)
# 3️⃣ CNV 状态列颜色：0 白色，1 三种颜色（每个CNV不同）
amp_colors <- list(
  "2p25_1_Amp"   = c("0" = "white", "1" = "#66C2A5"),
  "2q33_1_Amp"   = c("0" = "white", "1" = "#FC8D62"),
  "21q22_13_Amp" = c("0" = "white", "1" = "#8DA0CB")
)
# --- 合并所有注释颜色 ---
annotation_colors <- c(
  list(Disease = disease_colors),
  list(Malignancy = malignancy_colors),
  amp_colors
)

row_colors <- c(
  "2p25_1_amp" = "#66C2A5",
  "2q33_1_amp" = "#FC8D62",
  "21q22_13_amp" = "#8DA0CB"
)

row_annotation_colors <- list(cnv = row_colors)

colors_combined = colorRampPalette(brewer.pal(8,'RdBu'))(100)

start_color_1 <- "#FFFFFF"
end_color_1 <- colors_combined[60]

# 生成渐变函数
grad_fun_1 <- colorRampPalette(c(start_color_1, end_color_1))
blues <- grad_fun_1(10)

start_color_2 <- colors_combined[40]
end_color_2 <- "#FFFFFF"

# 生成渐变函数
grad_fun_2 <- colorRampPalette(c(start_color_2, end_color_2))
reds <- grad_fun_2(10)

# --- 绘图 ---
pheatmap::pheatmap(
  top_genes_avg_exp,
  cluster_rows = FALSE,
  scale = "row",
  color = rev(c(colors_combined[1:40],reds,blues,colors_combined[60:100])),
  annotation_row = cnv_anno,
  annotation_col = sample_anno[,rev(colnames(sample_anno))],
  annotation_colors = c(annotation_colors, row_annotation_colors),
  border_color = NA,
  show_colnames = TRUE,
  show_rownames = TRUE,
  clustering_distance_cols = "canberra"
)

sample_anno$test1 <- 0
sample_anno$test2 <- seq(1:78)

samplenames <- rownames(anno_sample_cluster_extended)
cnvscorepath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.xlsx",sep = "")
i=1
cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
cnv_mean <- cnv_score %>%
  group_by(group2) %>%
  summarise(mean_cnvscore = mean(cnvscore, na.rm = TRUE))
cnv_score_sample <- cnv_mean$mean_cnvscore[2] - cnv_mean$mean_cnvscore[1]
for (i in 2:length(samplenames)) {
  cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
  cnv_mean <- cnv_score %>%
    group_by(group2) %>%
    summarise(mean_cnvscore = mean(cnvscore, na.rm = TRUE))
  cnv_score_sample_1 <- cnv_mean$mean_cnvscore[2] - cnv_mean$mean_cnvscore[1]
  cnv_score_sample <- rbind(cnv_score_sample,cnv_score_sample_1)
  print(i)
}
cnv_score_sample <- as.data.frame(cnv_score_sample)
rownames(cnv_score_sample) <- samplenames

pheatmap::pheatmap(
  t(cbind(cnv_score_sample, cnv_score_sample)),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  color = viridisLite::magma(100),
  annotation_col = sample_anno[, rev(colnames(sample_anno))[-c(1, 2)]],
  annotation_colors = annotation_colors,
  border_color = NA
)

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

STS_group <- as.data.frame(STS_meta$"2p25_1_Amp")
STS_group$barcodekey <- rownames(STS_meta)
rownames(STS_group) <- STS_group$barcodekey
AUC_mtx <- fread("/cluster3/yflu/STS/scenic/auc/all_auc_40_transposed.csv",sep = ",",header = T)
rownames(AUC_mtx) <- AUC_mtx$V1

AUC_mtx_1 <- AUC_mtx

AUC_mtx <- AUC_mtx[,-1]
AUC_mtx <- as.matrix(as.data.frame(AUC_mtx))
rownames(AUC_mtx) <- rownames(AUC_mtx_1)
rownames(STS_group) <- STS_group$barcodekey
group_df <- STS_group[, "STS_meta$\"2p25_1_Amp\"", drop = FALSE]
colnames(group_df) <- "Group"
group_vec <- as.numeric(group_df$Group)
names(group_vec) <- rownames(group_df)  # 如果有行名

# 再运行 calcRSS
rss <- calcRSS(AUC_mtx, group_vec)
colnames(rss) <- paste("2p25_1_Amp",colnames(rss),sep = "_")

STS_group <- as.data.frame(STS_meta$"21q22_13_Amp")
STS_group$barcodekey <- rownames(STS_meta)
rownames(STS_group) <- STS_group$barcodekey
group_df <- STS_group[, "STS_meta$\"21q22_13_Amp\"", drop = FALSE]
colnames(group_df) <- "Group"
group_vec <- as.numeric(group_df$Group)
rss_1 <- calcRSS(AUC_mtx, group_vec)
colnames(rss_1) <- paste("21q22_13_Amp",colnames(rss_1),sep = "_")
rss <- cbind(rss, rss_1)

STS_group <- as.data.frame(STS_meta$"2q33_1_Amp")
STS_group$barcodekey <- rownames(STS_meta)
rownames(STS_group) <- STS_group$barcodekey
group_df <- STS_group[, "STS_meta$\"2q33_1_Amp\"", drop = FALSE]
colnames(group_df) <- "Group"
group_vec <- as.numeric(group_df$Group)
rss_1 <- calcRSS(AUC_mtx, group_vec)
colnames(rss_1) <- paste("2q33_1_Amp",colnames(rss_1),sep = "_")
rss <- cbind(rss, rss_1)
rss <- as.data.frame(rss)

rss$gene <- rownames(rss)
rss <- rss[,c(1,2,5,6,3,4,7)]

plot_rss <- function(
    rss,
    highlight_regulons = NULL,
    regulon_colors = NULL,   # ← 为每个高亮 regulon 指定颜色
    group_colors = NULL,     # ← 为每个 group 指定颜色
    title = "Regulon RSS ranking"
) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(ggrepel)
  
  # --- 数据整理 ---
  rss_df <- as.data.frame(rss)
  rss_df$Regulon <- rownames(rss_df)
  rss_long <- rss_df %>%
    pivot_longer(
      cols = -Regulon,
      names_to = "Group",
      values_to = "RSS"
    )  
  rss_long <- rss_long %>%
    dplyr::group_by(Group) %>%
    dplyr::arrange(desc(RSS)) %>%
    dplyr::mutate(Rank = row_number())
  
  # --- 提取需要标注的行 ---
  label_df <- NULL
  if (!is.null(highlight_regulons)) {
    label_df <- rss_long %>%
      filter(Regulon %in% highlight_regulons)
    if (nrow(label_df) == 0) {
      warning("⚠️ No regulons matched 'highlight_regulons'. Check naming with rownames(rss).")
    }
  }
  
  # --- 基础曲线 ---
  p <- ggplot(rss_long, aes(x = Rank, y = RSS, color = Group)) +
    geom_line(size = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = "Rank of regulons",
      y = "RSS score",
      color = "Group"
    ) +
    theme(
      legend.position = "top",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  # --- 设置 Group 颜色 ---
  if (!is.null(group_colors)) {
    p <- p + scale_color_manual(values = group_colors)
  }
  
  # --- 标注部分 ---
  if (!is.null(label_df) && nrow(label_df) > 0) {
    # 如果指定了 regulon 颜色，则按颜色绘制点
    if (!is.null(regulon_colors)) {
      label_df$LabelColor <- regulon_colors[label_df$Regulon]
      p <- p +
        geom_point(
          data = label_df,
          aes(x = Rank, y = RSS),
          size = 3,
          color = label_df$LabelColor
        ) +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = Regulon),
          color = label_df$LabelColor,
          size = 4,
          nudge_y = 0.005,
          max.overlaps = Inf
        )
    } else {
      # 否则用默认黑色
      p <- p +
        geom_point(data = label_df, aes(x = Rank, y = RSS), size = 3, color = "black") +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = Regulon),
          size = 4,
          color = "black",
          nudge_y = 0.005,
          max.overlaps = Inf
        )
    }
  }
  
  return(p)
}

my_colors <- c("2p25_1_Amp_0" = "#CCEADF", "2p25_1_Amp_1" = "#66C2A5",
               "2q33_1_Amp_0" = "#FED9CA", "2q33_1_Amp_1" = "#FC8D62",
               "21q22_13_Amp_0" = "#D9DFED", "21q22_13_Amp_1" = "#8DA0CB")
my_regulon_colors <- c(
  "SOX11(+)" = "#66C2A5",
  "OLIG1(+)" = "#8DA0CB",
  "OLIG2(+)" = "#8DA0CB",
  "RUNX1(+)" = "#8DA0CB",
  "SIM2(+)" = "#8DA0CB"
)
# 绘制 RSS 曲线并标注指定 regulons
plot_rss(
  rss[,-7],
  highlight_regulons = c("SOX11(+)", "OLIG1(+)", "OLIG2(+)","RUNX1(+)","SIM2(+)"),
  group_colors = my_colors,
  regulon_colors = my_regulon_colors,
  title = "CNV TF RSS ranking curve"
)

regulon_targets <- read.csv("/cluster3/yflu/STS/scenic/all_reg_cell_40.csv")
colnames(regulon_targets) <- c(regulon_targets[2,c(1,2)],regulon_targets[1,-c(1,2)])
regulon_targets <- regulon_targets[-c(1,2),]

sox11 <- subset(regulon_targets,TF == "SOX11")
olig1 <- subset(regulon_targets,TF == "OLIG1")
olig2 <- subset(regulon_targets,TF == "OLIG2")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

parse_targetgenes <- function(tg_string) {
  if(is.na(tg_string) || tg_string == "") return(NULL)
  
  # 去掉最外面的方括号
  tg_clean <- str_remove_all(tg_string, "^\\[|\\]$")
  
  # 拆分成单个 ('Gene', Score) 条目
  tg_items <- str_split(tg_clean, "\\), \\(")[[1]]
  
  # 去掉每个条目的多余符号
  tg_items <- str_replace_all(tg_items, "[\\(\\)' ]", "")
  
  # 按逗号拆分成 gene 和 score
  df <- data.frame(
    Gene = sapply(str_split(tg_items, ","), `[`, 1),
    Score = as.numeric(sapply(str_split(tg_items, ","), `[`, 2))
  )
  
  return(df)
}

# 对整列展开
sox11_long <- sox11 %>%
  dplyr::mutate(row_id = row_number()) %>%  # 保留行号
  dplyr::group_by(row_id, TF = TF) %>%
  group_modify(~ {
    df <- parse_targetgenes(.x$TargetGenes)
    if(is.null(df)) return(NULL)
    df
  }) %>%
  ungroup()
olig1_long <- olig1 %>%
  dplyr::mutate(row_id = row_number()) %>%  # 保留行号
  dplyr::group_by(row_id, TF = TF) %>%
  group_modify(~ {
    df <- parse_targetgenes(.x$TargetGenes)
    if(is.null(df)) return(NULL)
    df
  }) %>%
  ungroup()
olig2_long <- olig2 %>%
  dplyr::mutate(row_id = row_number()) %>%  # 保留行号
  dplyr::group_by(row_id, TF = TF) %>%
  group_modify(~ {
    df <- parse_targetgenes(.x$TargetGenes)
    if(is.null(df)) return(NULL)
    df
  }) %>%
  ungroup()

sox11_targets <- unique(sox11_long$Gene)
olig1_targets <- unique(olig1_long$Gene)
olig2_targets <- unique(olig2_long$Gene)

degs <- cnv_marker_list_all[["2p25_1_Amp"]]
degs <- subset(degs,cluster == "1")
rownames(degs) <- degs$gene
degs_sub <- degs %>% 
  filter(gene %in% sox11_targets & (avg_log2FC > 1|avg_log2FC < -1))

# 2️⃣ 提取基因名
sox11_targets <- degs_sub$gene

degs <- cnv_marker_list_all[["21q22_13_Amp"]]
degs <- subset(degs,cluster == "1")
rownames(degs) <- degs$gene
degs_sub <- degs %>% 
  filter(gene %in% olig1_targets & (avg_log2FC > 1|avg_log2FC < -1))
olig1_targets <- degs_sub$gene

degs_sub <- degs %>% 
  filter(gene %in% olig2_targets & (avg_log2FC > 1|avg_log2FC < -1))
olig2_targets <- degs_sub$gene

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

genelist_sox11 <- bitr(sox11_targets, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_olig1 <- bitr(olig1_targets, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_olig2 <- bitr(olig2_targets, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go_sox11 <- enrichGO(genelist_sox11$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,keyType = 'ENTREZID',readable = T)
go_olig1 <- enrichGO(genelist_olig1$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,keyType = 'ENTREZID',readable = T)
go_olig2 <- enrichGO(genelist_olig2$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,keyType = 'ENTREZID',readable = T)
go_sox11_res <- go_sox11@result 
go_olig1_res <- go_olig1@result 
go_olig2_res <- go_olig2@result 

go_sox11_res <- subset(go_sox11_res,p.adjust < 0.05&qvalue < 0.05)
go_olig1_res <- subset(go_olig1_res,p.adjust < 0.05&qvalue < 0.05)
go_olig2_res <- subset(go_olig2_res,p.adjust < 0.05&qvalue < 0.05)

gene_ids <- as.character(genelist_sox11$ENTREZID)

# KEGG 富集分析
kegg_sox11 <- enrichKEGG(
  gene         = gene_ids,
  organism     = "hsa",        # hsa 表示人
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kegg_sox11_res <- kegg_sox11@result
kegg_sox11_res <- subset(kegg_sox11_res,p.adjust < 0.05&qvalue < 0.05)

gene_ids <- as.character(genelist_olig1$ENTREZID)

# KEGG 富集分析
kegg_olig1 <- enrichKEGG(
  gene         = gene_ids,
  organism     = "hsa",        # hsa 表示人
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_olig1_res <- kegg_olig1@result
kegg_olig1_res <- subset(kegg_olig1_res,p.adjust < 0.05&qvalue < 0.05)

gene_ids <- as.character(genelist_olig2$ENTREZID)

# KEGG 富集分析
kegg_olig2 <- enrichKEGG(
  gene         = gene_ids,
  organism     = "hsa",        # hsa 表示人
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_olig2_res <- kegg_olig2@result
kegg_olig2_res <- subset(kegg_olig2_res,p.adjust < 0.05&qvalue < 0.05)

library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(purrr)

get_genes_list <- function(go_res = NULL, kegg_res = NULL) {
  
  # =========================
  # 处理 GO
  # =========================
  go_list <- list()
  if (!is.null(go_res)) {
    
    # GO ID -> 基因 SYMBOL
    get_go_genesymbols <- function(go_id) {
      # 获取该 term 下所有 EntrezID（包括间接注释）
      entrez_ids <- tryCatch({
        suppressMessages(
          AnnotationDbi::select(org.Hs.eg.db,
                                keys = go_id,
                                keytype = "GOALL",  # 包含所有层级的 GO 注释
                                columns = "ENTREZID")$ENTREZID
        )
      }, error = function(e) NULL)
      
      if (is.null(entrez_ids) || length(entrez_ids) == 0)
        return(character(0))
      
      syms <- suppressMessages(
        bitr(entrez_ids,
             fromType = "ENTREZID",
             toType = "SYMBOL",
             OrgDb = org.Hs.eg.db)
      )
      
      unique(syms$SYMBOL)
    }
    
    go_genes <- map(go_res$ID, get_go_genesymbols)
    go_list <- setNames(go_genes, go_res$Description)
  }
  
  # =========================
  # 处理 KEGG
  # =========================
  kegg_list <- list()
  if (!is.null(kegg_res)) {
    
    get_kegg_genesymbols <- function(kegg_id){
      pathway <- tryCatch(KEGGREST::keggGet(kegg_id)[[1]], error = function(e) NULL)
      if (is.null(pathway) || is.null(pathway$GENE))
        return(character(0))
      
      entrez_ids <- pathway$GENE[seq(1, length(pathway$GENE), 2)]
      syms <- suppressMessages(
        bitr(entrez_ids,
             fromType = "ENTREZID",
             toType = "SYMBOL",
             OrgDb = org.Hs.eg.db)
      )
      unique(syms$SYMBOL)
    }
    
    kegg_genes <- map(kegg_res$ID, get_kegg_genesymbols)
    kegg_list <- setNames(kegg_genes, kegg_res$Description)
  }
  
  # =========================
  # 合并输出
  # =========================
  combined_list <- c(go_list, kegg_list)
  return(combined_list)
}


run_gsea_for_list <- function(degs, genes_list, pvalueCutoff = 1) {
  # 准备 geneList
  gene_list <- degs$avg_log2FC
  names(gene_list) <- degs$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # 对每个 term 做 GSEA
  gsea_results <- map(names(genes_list), function(term_name){
    geneset <- genes_list[[term_name]]
    
    # 跳过空基因集
    if(length(geneset) == 0) return(NULL)
    
    term_df <- data.frame(term = term_name, gene = geneset, stringsAsFactors = FALSE)
    
    res <- tryCatch({
      enr <- GSEA(geneList = gene_list,
                  TERM2GENE = term_df,
                  verbose = FALSE,
                  pvalueCutoff = pvalueCutoff)
      
      if(is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
      
      df <- as.data.frame(enr)
      # 有些列可能为空，选取常规列
      df[, c("ID","Description","NES","pvalue","p.adjust")]
    }, error = function(e) NULL)
    
    return(res)
  })
  
  # 去掉 NULL
  gsea_results <- gsea_results[!sapply(gsea_results, is.null)]
  
  # 合并
  gsea_df <- if(length(gsea_results) > 0) bind_rows(gsea_results) else data.frame()
  
  return(gsea_df)
}

genes_list_sox11 <- get_genes_list(go_res = go_sox11_res, kegg_res = kegg_sox11_res)
degs <- cnv_marker_list_all[["2p25_1_Amp"]]
degs <- subset(degs,cluster == "1")
rownames(degs) <- degs$gene

gsea_df_sox11 <- run_gsea_for_list(degs, genes_list_sox11, pvalueCutoff = 1)
gsea_df_sox11 <- subset(gsea_df_sox11,p.adjust < 0.05)

genes_sox11 <- subset(go_sox11_res,Description %in% rownames(gsea_df_sox11))
gene_list <- strsplit(genes_sox11$geneID, "/")
gene_vec_sox11 <- unlist(gene_list)
gene_vec_sox11 <- unique(gene_vec_sox11)

go_sox11_res_sub <- subset(go_sox11_res,Description %in% gsea_df_sox11$Description)

genes_list_olig1 <- get_genes_list(go_res = go_olig1_res, kegg_res = kegg_olig1_res)
degs <- cnv_marker_list_all[["21q22_13_Amp"]]
degs <- subset(degs,cluster == "1")
rownames(degs) <- degs$gene

gsea_df_olig1 <- run_gsea_for_list(degs, genes_list_olig1, pvalueCutoff = 1)
gsea_df_olig1 <- subset(gsea_df_olig1,p.adjust < 0.05)

genes_olig1 <- subset(go_olig1_res,Description %in% rownames(gsea_df_olig1))
gene_list <- strsplit(genes_olig1$geneID, "/")
gene_vec_olig1 <- unlist(gene_list)
gene_vec_olig1 <- unique(gene_vec_olig1)

go_olig1_res_sub <- subset(go_olig1_res,Description %in% gsea_df_olig1$Description)

genes_list_olig2 <- get_genes_list(go_res = go_olig2_res, kegg_res = kegg_olig2_res)
degs <- cnv_marker_list_all[["21q22_13_Amp"]]
degs <- subset(degs,cluster == "1")
rownames(degs) <- degs$gene

gsea_df_olig2 <- run_gsea_for_list(degs, genes_list_olig2, pvalueCutoff = 1)
gsea_df_olig2 <- subset(gsea_df_olig2,p.adjust < 0.05)

genes_olig2 <- subset(go_olig2_res,Description %in% rownames(gsea_df_olig2))
gene_list <- strsplit(genes_olig2$geneID, "/")
gene_vec_olig2 <- unlist(gene_list)
gene_vec_olig2 <- unique(gene_vec_olig2)

go_olig2_res_sub <- subset(go_olig2_res,Description %in% gsea_df_olig2$Description)

plot_go_gene_network <- function(gene_vec, go_res, degs, gsea_df=NULL, highlight_gene="SOX11") {
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(ggnewscale)
  
  # -------------------------
  # 构建边
  # -------------------------
  # GO–gene edges
  go_long <- go_res %>%
    separate_rows(geneID, sep="/") %>%
    rename(Gene=geneID)
  edges_go_gene <- go_long %>% dplyr::select(from=Description, to=Gene)
  
  # highlight_gene–gene edges
  edges_center <- data.frame(from=highlight_gene, to=setdiff(gene_vec, highlight_gene))
  
  # 合并
  edges <- bind_rows(edges_go_gene, edges_center)
  
  # -------------------------
  # 构建节点属性
  # -------------------------
  nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
  nodes <- nodes %>%
    mutate(
      type = case_when(
        name == highlight_gene ~ "highlight",
        name %in% go_res$Description ~ "GO_term",
        TRUE ~ "Gene"
      ),
      node_size = case_when(
        name == highlight_gene ~ 12,
        TRUE ~ 5
      )
    ) %>%
    # 给 gene 节点添加 avg_log2FC
    left_join(degs %>% select(gene, avg_log2FC), by=c("name"="gene")) %>%
    # 给 GO term 节点添加 NES
    left_join(gsea_df %>% select(ID, NES), by=c("name"="ID"))
  
  # 如果 highlight_gene 的 avg_log2FC 是 NA，设置为0
  nodes$avg_log2FC[is.na(nodes$avg_log2FC) & nodes$type=="highlight"] <- 0
  
  # -------------------------
  # 构建 igraph
  # -------------------------
  g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  layout <- create_layout(as_tbl_graph(g), layout="fr")
  
  # -------------------------
  # 绘图
  # -------------------------
  ggraph(layout) +
    geom_edge_link(alpha = 0.4, color="grey70") +
    
    # gene 节点（包括高亮） 蓝白红
    geom_node_point(
      data = layout %>% filter(type %in% c("Gene","highlight")),
      aes(x=x, y=y, color=avg_log2FC, size=node_size),
      shape=16
    ) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0, name="log2FC") +
    
    # 新 color scale 用于 GO term
    ggnewscale::new_scale_color() +
    
    # GO term 节点 黄-白-紫
    geom_node_point(
      data = layout %>% filter(type=="GO_term"),
      aes(x=x, y=y, color=NES, size=node_size),
      shape=15
    ) +
    scale_color_gradient2(low="yellow", mid="white", high="purple", midpoint=0, name="NES") +
    
    geom_node_text(aes(x=x, y=y, label=name), repel=TRUE, size=3) +
    scale_size_identity() +
    theme_void()
}

degs_sox11 <- cnv_marker_list_all[["2p25_1_Amp"]]
degs_sox11 <- subset(degs_sox11,gene %in% gene_vec_sox11&cluster == "1")
rownames(degs_sox11) <- degs_sox11$gene
gsea_df_sox11

set.seed(123)
plot_go_gene_network(
  gene_vec = gene_vec_sox11,
  go_res = go_sox11_res_sub,
  degs = degs_sox11,
  gsea_df = gsea_df_sox11,
  highlight_gene = "SOX11"
)

gene_vec_olig1 <- c(gene_vec_olig1,"OLIG1")
degs_olig1 <- cnv_marker_list_all[["21q22_13_Amp"]]
degs_olig1 <- subset(degs_olig1,gene %in% gene_vec_olig1&cluster == "1")
rownames(degs_olig1) <- degs_olig1$gene

set.seed(123)
plot_go_gene_network(
  gene_vec = gene_vec_olig1,
  go_res = go_olig1_res_sub,
  degs = degs_olig1,
  gsea_df = gsea_df_olig1,
  highlight_gene = "OLIG1"
)

gene_vec_olig2 <- c(gene_vec_olig2,"OLIG2")
degs_olig2 <- cnv_marker_list_all[["21q22_13_Amp"]]
degs_olig2 <- subset(degs_olig2,gene %in% gene_vec_olig2&cluster == "1")
rownames(degs_olig2) <- degs_olig2$gene

set.seed(120)
plot_go_gene_network(
  gene_vec = gene_vec_olig2,
  go_res = go_olig2_res_sub,
  degs = degs_olig2,
  gsea_df = gsea_df_olig2,
  highlight_gene = "OLIG2"
)
