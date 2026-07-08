neighbor <- read.csv("celltype_neighbor_proportion_radius40um.csv")

library(umap)
library(cluster)
library(stats)
library(tibble)
library(dplyr)
# ===============================
# 参数设置
# ===============================
max_k <- 15            # silhouette 检测最大聚类数
sample_size <- 5000    # 每次采样细胞数
num_repeats <- 50       # 多次采样次数

feature_matrix <- neighbor %>%
  select(starts_with("prop_")) %>%   # 只取比例列
  as.matrix()
rownames(feature_matrix) <- neighbor$barcode
# 去掉 NA 或 NaN
feature_matrix[is.na(feature_matrix)] <- 0
feature_matrix[is.nan(feature_matrix)] <- 0
feature_matrix_clean <- feature_matrix

# ===============================
# 多次采样 + silhouette 估计最佳 k
# ===============================
sil_scores_all <- matrix(NA, nrow = max_k, ncol = num_repeats)

set.seed(42)
for (r in 1:num_repeats) {
  cat("采样轮次:", r, "\n")
  
  sample_idx <- sample(1:nrow(feature_matrix_clean), sample_size)
  feature_sample <- feature_matrix_clean[sample_idx, ]
  
  for (k in 2:max_k) {
    km <- kmeans(feature_sample, centers = k, nstart = 10)
    ss <- silhouette(km$cluster, dist(feature_sample))
    sil_scores_all[k, r] <- mean(ss[,3])
  }
}

# 计算平均 silhouette
sil_scores_mean <- rowMeans(sil_scores_all, na.rm = TRUE)
plot(2:max_k, sil_scores_mean[2:max_k], type="b", xlab="k", ylab="Average silhouette", main="Silhouette")
best_k <- which.max(sil_scores_mean)
cat("建议的最佳 k =", best_k, "\n")

# ===============================
# 全量数据 k-means
# ===============================

best_k = 5

cat("进行全量数据 k-means 聚类...\n")
km_full <- kmeans(feature_matrix_clean, centers = best_k, nstart = 20)
cluster_labels <- km_full$cluster

# 保存结果
feature_matrix_clean <- as.data.frame(feature_matrix_clean)
feature_matrix_clean$cluster <- cluster_labels
feature_matrix_clean$barcode <- neighbor$barcode

cluster_prop_mean <- feature_matrix_clean %>%
  group_by(cluster) %>%
  summarise(
    across(
      starts_with("prop_"),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

head(cluster_prop_mean)

library(dplyr)
library(tibble)
library(pheatmap)

# raw proportion
sig_raw <- cluster_prop_mean %>%
  column_to_rownames("cluster") %>%
  as.matrix()

# row-scaled version
sig_scaled <- t(scale(t(sig_raw)))   # row Z-score

sig_text_pct <- matrix(sprintf("%.1f%%", sig_raw * 100),
                       nrow = nrow(sig_raw),
                       ncol = ncol(sig_raw))
rownames(sig_text_pct) <- rownames(sig_raw)
colnames(sig_text_pct) <- colnames(sig_raw)

pheatmap(sig_scaled,
         display_numbers = sig_text_pct,
         fontsize_number = 6,
         clustering_method = "ward.D2",
         main = "Niche composition (row-scaled, numbers = raw proportion)")
feature_matrix_clean <- cbind(neighbor[,c(2:7)],feature_matrix_clean)
#write.csv(feature_matrix_clean,"feature_matrix_clean.csv")
feature_matrix_clean <- read.csv("feature_matrix_clean.csv")
library(dplyr)

# 计算每个 Sample 的每个 cluster 的细胞比例
cluster_prop <- feature_matrix_clean %>%
  group_by(Sample, cluster) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(cluster_prop = cell_count / sum(cell_count)) %>%
  ungroup()

# 查看结果
head(cluster_prop)

library(tidyr)

cluster_prop_wide <- cluster_prop %>%
  select(Sample, cluster, cluster_prop) %>%
  pivot_wider(names_from = cluster, values_from = cluster_prop, values_fill = 0)

head(cluster_prop_wide)

anno_TMA <- openxlsx::read.xlsx("/cluster3/yflu/STS/TMA/order.xlsx","STS_anno")
rownames(anno_TMA) <- anno_TMA$Sample
anno_TMA <- anno_TMA[unique(cluster_prop$Sample),]

cluster_prop_wide$Sample <- as.character(cluster_prop_wide$Sample)

anno_TMA$Sample <- rownames(anno_TMA)

dat <- cluster_prop_wide %>%
  left_join(anno_TMA, by = "Sample")

head(dat)

library(dplyr)
library(purrr)

cluster_cols <- setdiff(colnames(cluster_prop_wide), "Sample")

stat_res <- map_dfr(cluster_cols, function(cl){
  
  df <- dat %>% select(Sample, Malignancy, all_of(cl)) %>%
    rename(prop = all_of(cl))
  
  # Wilcoxon
  w <- wilcox.test(prop ~ Malignancy, data = df)
  
  # 计算组均值
  mean_benign <- mean(df$prop[df$Malignancy == "Benign"], na.rm = TRUE)
  mean_malignant <- mean(df$prop[df$Malignancy == "Malignant"], na.rm = TRUE)
  
  tibble(
    cluster = cl,
    mean_benign = mean_benign,
    mean_malignant = mean_malignant,
    diff_malignant_minus_benign = mean_malignant - mean_benign,
    p_value = w$p.value
  )
})

# FDR校正
stat_res <- stat_res %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

stat_res

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

dat <- cluster_prop_wide %>%
  left_join(anno_TMA, by = "Sample")

full_names <- c("Hemangioma", "KHE", "Schwannoma", "MPNST", "Undifferentiated sarcoma",
                "RMS", "MRT", "IMT", "Angiosarcoma", "EWS/PNET",
                "NF", "Aggressive fibromatosis", "Liposarcoma", "Spindle cell tumor", "ASPS",
                "Infantile fibrosarcoma", "Synovial sarcoma", "Lipoblastoma", "Pecoma", "Lymphangioma")

# 对应缩写
abbreviations <- c("HE","KHE","SWN","MPNST","US","RMS","MRT","IMT","AS","EWS",
                   "NF","AF","LPS","SCT","ASPS","IFS","SS","LPB","PECOMA","LYM")

name_map <- setNames(full_names,abbreviations)
dat$Disease <- as.character(dat$Disease)
dat$Disease <- name_map[dat$Disease]
dat$Disease <- as.character(dat$Disease)
aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]

library(pheatmap)
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))
dat$Disease <- factor(dat$Disease,levels = labels)

i=1
dat_1 <- dat
colnames(dat_1)[i+1] <- "CNV_score"
p <- plot_cnv_box_by_malignancy(dat_1, plot_title = paste("niche",i,sep = "_"))
p
for (i in 2:5) {
  dat_1 <- dat
  colnames(dat_1)[i+1] <- "CNV_score"
  p1 <- plot_cnv_box_by_malignancy(dat_1, plot_title = paste("niche",i,sep = "_"))
  p <- p+p1
}
p <- p + plot_layout(nrow = 1)
p

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

dat_sc <- as.data.frame(dat)
rownames(dat_sc) <- dat_sc$Sample
dat_sc <- dat_sc[rownames(anno_sample_cluster_extended),]
dat_sc <- na.omit(dat_sc)
rownames(anno_sample_cluster_extended)[10] <- "T947C"
dat_sc <- cbind(dat_sc,anno_sample_cluster_extended[rownames(dat_sc),c(3)])
colnames(dat_sc)[10] <- "Louvain"

colnames(dat_sc)[2:6] <- paste("niche",c(1:5),sep = "_")

dat_sc$Louvain <- as.character(dat_sc$Louvain)

i=1
dat_sc_1 <- dat_sc
colnames(dat_sc_1)[i+1] <- "niche"

P1 <- ggplot(dat_sc_1, mapping=aes(x=Louvain,y=niche,fill=Louvain))+ ##设置图形的纵坐标横坐标和分组
  stat_boxplot(mapping=aes(x=Louvain,y=niche),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
  geom_boxplot(aes(fill=Louvain),                             ##分组比较的变量
               position=position_dodge(0.8),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  ggtitle(colnames(dat_sc)[i+1])+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none",    
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
P1
louviannames <- unique(dat_sc$Louvain)
louviannames <- louviannames[order(louviannames)]
my_comparisons <- combn(louviannames,2,simplify = F)
my_comparisons_sig <- list()
for (i in 1:length(my_comparisons)) {
  por1 <- subset(dat_sc_1,Louvain == my_comparisons[[i]][1])$niche
  por2 <- subset(dat_sc_1,Louvain == my_comparisons[[i]][2])$niche
  
  if(length(por1) > 1 & length(por2) > 1){
    test <- t.test(por1, 
                   por2)
    if(is.na(test$p.value)) {
      test$p.value <- 1
    }
    if(test$p.value < 0.05){
      my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[i]]))
    }
  }
}
P1 <- P1 + stat_compare_means(comparisons=my_comparisons_sig,
                              label.y = seq(from=max(dat_sc_1$niche)+0.1, to=80, by=0.1),
                              method="t.test",
                              label="p.signif",hide.ns = T)
P1

for (i in 2:5) {
  dat_sc_1 <- dat_sc
  colnames(dat_sc_1)[i+1] <- "niche"
  
  P2 <- ggplot(dat_sc_1, mapping=aes(x=Louvain,y=niche,fill=Louvain))+ ##设置图形的纵坐标横坐标和分组
    stat_boxplot(mapping=aes(x=Louvain,y=niche),
                 geom ="errorbar",                             ##添加箱子的bar为最大、小值
                 width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
    geom_boxplot(aes(fill=Louvain),                             ##分组比较的变量
                 position=position_dodge(0.8),                 ##因为分组比较，需设组间距
                 width=0.6,                                    ##箱子的宽度
                 outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
    geom_jitter(aes(fill=Louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
    ggtitle(colnames(dat_sc)[i+1])+ #设置总的标题
    theme_bw()+ #背景变为白色
    theme(legend.position="none",    
          panel.grid.major = element_blank(), #不显示网格线
          panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  P2
  louviannames <- unique(dat_sc$Louvain)
  louviannames <- louviannames[order(louviannames)]
  my_comparisons <- combn(louviannames,2,simplify = F)
  my_comparisons_sig <- list()
  for (i in 1:length(my_comparisons)) {
    por1 <- subset(dat_sc_1,Louvain == my_comparisons[[i]][1])$niche
    por2 <- subset(dat_sc_1,Louvain == my_comparisons[[i]][2])$niche
    
    if(length(por1) > 1 & length(por2) > 1){
      test <- t.test(por1, 
                     por2)
      if(is.na(test$p.value)) {
        test$p.value <- 1
      }
      if(test$p.value < 0.05){
        my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[i]]))
      }
    }
  }
  P2 <- P2 + stat_compare_means(comparisons=my_comparisons_sig,
                                label.y = seq(from=max(dat_sc_1$niche)+0.1, to=80, by=0.1),
                                method="t.test",
                                label="p.signif",hide.ns = T)
  P1 <- P1+P2
}
P1

TMA_merged_NMF_score <- readRDS("/cluster3/yflu/STS/TMA/TMA_merged_NMF_score.rds")
TMA_merged_NMF_score_tumor <- subset(TMA_merged_NMF_score,celltype == "Tumor cells")

TMA_merged_NMF_score_tumor_avag <- AverageExpression(TMA_merged_NMF_score_tumor,assays = "AUCell",group.by = "Sample")
TMA_merged_NMF_score_tumor_avag <- as.data.frame(TMA_merged_NMF_score_tumor_avag$AUCell)
TMA_merged_NMF_score_tumor_avag <- t(TMA_merged_NMF_score_tumor_avag)

dat <- as.data.frame(dat)
rownames(dat) <- dat$Sample
TMA_merged_NMF_score_tumor_avag_dat <- cbind(TMA_merged_NMF_score_tumor_avag,dat[rownames(TMA_merged_NMF_score_tumor_avag),])

prog_cols <- paste0("program", LETTERS[1:8])   # programA–H
meta_cols <- c("Sample", "Disease", "Malignancy", "TMA")
niche_cols <- as.character(1:5)               # "1","2","3","4","5"

dat_cor <- TMA_merged_NMF_score_tumor_avag_dat %>%
  dplyr::select(all_of(prog_cols), all_of(niche_cols))

cor_mat <- cor(dat_cor[, prog_cols],
               dat_cor[, niche_cols],
               method = "pearson",
               use = "pairwise.complete.obs")

cor_mat

library(Hmisc)

res <- rcorr(as.matrix(dat_cor), type = "pearson")

# 提取 program vs niche p-values
p_mat <- res$P[prog_cols, niche_cols]
p_mat

p_vec <- as.vector(p_mat)
p_adj <- p.adjust(p_vec, method = "BH")
p_adj_mat <- matrix(p_adj, nrow = nrow(p_mat), dimnames = dimnames(p_mat))
p_adj_mat

library(pheatmap)

# 对称色阶（以最大绝对值为界）
max_abs <- max(abs(cor_mat), na.rm = TRUE)

bk <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(cor_mat,
         color = colorRampPalette(c("blue","white","red"))(100),
         breaks = bk,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_mat[p_adj_mat  < 0.05]  <- "*"
sig_mat[p_adj_mat  < 0.01]  <- "**"
sig_mat[p_adj_mat  < 0.001] <- "***"

rownames(sig_mat) <- rownames(p_mat)
colnames(sig_mat) <- colnames(p_mat)

pheatmap(cor_mat,
         color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
         breaks = bk,
         display_numbers = sig_mat,
         fontsize_number = 10,
         fontsize_row = 12,
         fontsize_col = 12,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F)

distance <- read.csv("tumor_cell_min_distance_wide.csv")
rownames(feature_matrix_clean) <- feature_matrix_clean$barcode
rownames(distance) <- distance$cell_id

library(dplyr)
library(tidyr)

cell_types <- colnames(distance)[3:17]   # 距离列

dist_long <- distance %>%
  pivot_longer(cols = all_of(cell_types),
               names_to = "target_celltype",
               values_to = "distance")

library(ggplot2)
library(ggridges)

library(dplyr)
library(purrr)

# 计算每个 cluster + target_celltype 的 density peak
peak_df <- dist_long %>%
  filter(!is.na(distance)) %>%
  group_by(target_celltype, cluster) %>%
  summarise(
    dens = list(density(distance, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    mode = map_dbl(dens, ~ .x$x[which.max(.x$y)]),
    peak_height = map_dbl(dens, ~ max(.x$y))
  ) %>%
  select(-dens)

library(ggplot2)
library(ggridges)

peak_df <- peak_df %>%
  mutate(label = sprintf("Peak=%.1f\nDensity=%.2f", mode, peak_height))

p <- ggplot(dist_long, aes(x = distance, y = factor(cluster), fill = factor(cluster))) +
  geom_density_ridges(alpha = 0.6, scale = 2, color = "white") +
  facet_wrap(~ target_celltype, scales = "free_x", ncol = 4) +
  theme_ridges() +
  theme_bw() +
  labs(x = "nearest distance", y = "cluster") +
  guides(fill = "none") +
  coord_cartesian(xlim = c(-50, 300)) + 
  geom_text(
    data = peak_df,
    aes(x = mode, y = as.numeric(cluster) + 0.25, 
        label = sprintf("Peak=%.1f", mode)),
    size = 2.5, color = "black"
  )  +
  geom_text(data = peak_df,
            aes(x = mode, y = as.numeric(cluster) + 0.3, label = label),
            size = 2.2)

TMA_merged_NMF_score_tumor_avag_dat

sample_mean_distance <- distance %>%
  group_by(Sample) %>%
  summarise(across(Adipocytes:Tumor.cells, 
                   ~ mean(.x[.x != 0], na.rm = TRUE)))

distance_nmf <- cbind(TMA_merged_NMF_score_tumor_avag_dat,sample_mean_distance)

prog_mat  <- distance_nmf[, 1:8]
dist_mat  <- distance_nmf[, 19:33]

# 保证都是数值
prog_mat <- as.data.frame(lapply(prog_mat, as.numeric))
dist_mat <- as.data.frame(lapply(dist_mat, as.numeric))

cor_mat <- cor(prog_mat, dist_mat, use = "pairwise.complete.obs", method = "spearman")

p_mat <- matrix(NA, nrow = ncol(prog_mat), ncol = ncol(dist_mat))
rownames(p_mat) <- colnames(prog_mat)
colnames(p_mat) <- colnames(dist_mat)

for(i in 1:ncol(prog_mat)){
  for(j in 1:ncol(dist_mat)){
    tmp <- cor.test(prog_mat[[i]], dist_mat[[j]], method = "spearman", use = "pairwise.complete.obs")
    p_mat[i, j] <- tmp$p.value
  }
}

density <- read.csv("/cluster3/yflu/STS/TMA/tumor_cell_weighted_density.csv")
rownames(density) <- density$cell_id

library(dplyr)

# 假设 density 是你的数据框
density_means <- density %>%
  group_by(Sample) %>%                  # 按 Sample 分组
  summarise(across(1:15, mean, na.rm = TRUE))  # 计算前15列的均值，忽略NA

# 查看结果
density_means <- as.data.frame(density_means)
rownames(density_means) <- density_means$Sample

density_means_nmf <- cbind(density_means,TMA_merged_NMF_score_tumor_avag_dat)

library(ggplot2)
library(reshape2)

# 提取两组数据
group1 <- density_means_nmf[, 2:16]
group2 <- density_means_nmf[, 17:24]

# 计算相关性矩阵
cor_mat <- cor(group1, group2, use = "pairwise.complete.obs", method = "spearman")

# 计算 p 值矩阵
p_mat <- matrix(NA, nrow = ncol(group1), ncol = ncol(group2))
rownames(p_mat) <- colnames(group1)
colnames(p_mat) <- colnames(group2)

for(i in 1:ncol(group1)) {
  for(j in 1:ncol(group2)) {
    test <- cor.test(group1[[i]], group2[[j]], method = "spearman")
    p_mat[i,j] <- test$p.value
  }
}

# 生成显著性标记
sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_mat[p_mat < 0.05] <- "*"
sig_mat[p_mat < 0.01] <- "**"
sig_mat[p_mat < 0.001] <- "***"
rownames(sig_mat) <- rownames(p_mat)
colnames(sig_mat) <- colnames(p_mat)

# 将矩阵转为长格式
cor_df <- melt(cor_mat)
colnames(cor_df) <- c("CellType", "Program", "Correlation")

sig_df <- melt(sig_mat)
colnames(sig_df) <- c("CellType", "Program", "Signif")

plot_df <- merge(cor_df, sig_df, by = c("CellType", "Program"))

# 画热图
ggplot(plot_df, aes(x = Program, y = CellType, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), color = "black", size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1,1), name = "Pearson r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  labs(x = "Program", y = "Cell Type", title = "Correlation between Cell Types and Programs")

feature_matrix_clean <- read.csv("feature_matrix_clean.csv")
rownames(feature_matrix_clean) <- feature_matrix_clean$barcode
density_cluster <- cbind(density,feature_matrix_clean[rownames(density),c(22,23)])

library(dplyr)

# 计算每个 Sample 和 cluster 的平均值
density_cluster_means <- density_cluster %>%
  group_by(Sample, cluster) %>%              # 按 Sample + cluster 分组
  summarise(across(1:15, mean, na.rm = TRUE))  # 计算前15列平均值

density_cluster_means <- as.data.frame(density_cluster_means)

TMA_merged_NMF_score_tumor_dat <- as.matrix(TMA_merged_NMF_score_tumor@assays$AUCell@data)
TMA_merged_NMF_score_tumor_dat <- as.data.frame(t(TMA_merged_NMF_score_tumor_dat))

TMA_merged_NMF_score_tumor_dat <- cbind(TMA_merged_NMF_score_tumor_dat,feature_matrix_clean[rownames(TMA_merged_NMF_score_tumor_dat),c(3,22,23)])

TMA_merged_NMF_score_tumor_dat_means <- TMA_merged_NMF_score_tumor_dat %>%
  group_by(Sample, cluster) %>%              # 按 Sample + cluster 分组
  summarise(across(1:8, mean, na.rm = TRUE))  # 计算前15列平均值

density_cluster_nmf <- cbind(density_cluster_means, TMA_merged_NMF_score_tumor_dat_means)
density_cluster_nmf <- density_cluster_nmf[,-c(2)]

library(dplyr)

# 去掉重复列名
df <- density_cluster_nmf[, !duplicated(colnames(density_cluster_nmf))]

clusters <- unique(df$cluster)
cor_list <- list()

for(cl in clusters){
  
  # 当前 cluster
  df_cl <- df %>% filter(cluster == cl)
  
  # 按 Sample 聚合
  df_mean <- df_cl %>%
    group_by(Sample) %>%
    summarise(across(1:15, mean, na.rm = TRUE),
              across(18:25, mean, na.rm = TRUE)) %>%
    ungroup()
  
  group1 <- df_mean[, 2:16]
  group2 <- df_mean[, 17:24]
  
  # 初始化矩阵
  cor_mat <- matrix(NA, ncol = ncol(group2), nrow = ncol(group1))
  p_mat   <- matrix(NA, ncol = ncol(group2), nrow = ncol(group1))
  
  rownames(cor_mat) <- colnames(group1)
  colnames(cor_mat) <- colnames(group2)
  rownames(p_mat)   <- colnames(group1)
  colnames(p_mat)   <- colnames(group2)
  
  # 两两相关性
  for(i in 1:ncol(group1)){
    for(j in 1:ncol(group2)){
      
      x <- group1[[i]]
      y <- group2[[j]]
      
      # 去掉 celltype=0 的 sample
      idx <- which(x != 0 & !is.na(x) & !is.na(y))
      
      # 至少 3 个样本才计算
      if(length(idx) >= 3){
        ct <- cor.test(x[idx], y[idx], method = "spearman")
        cor_mat[i, j] <- ct$estimate
        p_mat[i, j]   <- ct$p.value
      }
    }
  }
  
  cor_list[[as.character(cl)]] <- list(
    correlation = cor_mat,
    p_value = p_mat
  )
}

# 查看
names(cor_list)
cor_list[[1]]$correlation
cor_list[[1]]$p_value

library(dplyr)

TMA_prog_mean <- TMA_merged_NMF_score_tumor_dat_means %>%
  group_by(cluster, Sample) %>% 
  summarise(across(programA:programH, mean, na.rm = TRUE)) %>%
  ungroup()

library(dplyr)
library(tidyr)

prog_long <- TMA_prog_mean %>%
pivot_longer(cols = programA:programH,
names_to = "Program",
values_to = "Score")

prog_wide <- prog_long %>%
  pivot_wider(names_from = cluster, values_from = Score)

prog_plot_long <- prog_wide %>%
  pivot_longer(
    cols = -c(Sample, Program),
    names_to = "cluster",
    values_to = "Score"
  )

# 固定 cluster 顺序
cluster_levels <- c("1","2","3","4","5")
prog_plot_long$cluster <- factor(prog_plot_long$cluster, levels = cluster_levels)
prog_plot_long$Program <- factor(prog_plot_long$Program)

library(purrr)

pairwise_list <- prog_wide %>%
  group_by(Program) %>%
  group_split() %>%
  map(function(df){
    
    clusters <- colnames(df)[-(1:2)]  # Sample, Program
    
    res <- combn(clusters, 2, simplify = FALSE, FUN = function(cl){
      
      x <- df[[cl[1]]]
      y <- df[[cl[2]]]
      idx <- which(!is.na(x) & !is.na(y))
      
      if(length(idx) < 3) return(NULL)
      
      wt <- wilcox.test(x[idx], y[idx], paired = TRUE)
      
      data.frame(
        Program = unique(df$Program),
        group1 = cl[1],
        group2 = cl[2],
        n = length(idx),
        mean_group1 = mean(x[idx]),
        mean_group2 = mean(y[idx]),
        delta_mean = mean(x[idx] - y[idx]),
        p = wt$p.value
      )
    })
    
    bind_rows(res)
  }) %>% bind_rows() %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni"))


stat_df <- pairwise_list %>%
  filter(p.adj < 0.05) %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*"
    )
  )
stat_df <- stat_df %>%
  mutate(
    group1 = factor(group1, levels = cluster_levels),
    group2 = factor(group2, levels = cluster_levels),
    g1 = as.numeric(as.character(group1)),
    g2 = as.numeric(as.character(group2))
  ) %>%
  arrange(Program, g1, g2)
stat_df <- stat_df %>%
  left_join(
    prog_plot_long %>%
      group_by(Program) %>%
      summarise(ymax = max(Score, na.rm = TRUE)),
    by = "Program"
  ) %>%
  group_by(Program) %>%
  mutate(y.position = ymax * (1 + 0.1 * row_number())) %>%
  ungroup()
cluster_cols <- c(
  "1" = "#F4FBD2",
  "2" = "#7CAAD0",
  "3" = "#DD4030",
  "4" = "#FEE79B",
  "5" = "#4575B4"
)
library(ggplot2)
library(ggpubr)

ggplot(prog_plot_long, aes(cluster, Score, fill = cluster)) +
  
  geom_boxplot(outlier.shape = NA) +
  
  # facet_wrap 2行4列
  facet_wrap(~Program, scales = "free_y", nrow = 2, ncol = 4) +
  
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  
  stat_pvalue_manual(
    stat_df,
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    label = "label",
    inherit.aes = FALSE,
    hide.ns = TRUE
  ) +
  
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

density_cluster_means_long <- density_cluster_means %>%
  pivot_longer(
    cols = -c(Sample, cluster),
    names_to = "celltype",
    values_to = "density"
  )

density_cluster_wide <- density_cluster_means_long %>%
  pivot_wider(
    id_cols = c(Sample, celltype),   # 每个 Sample + celltype 是一行
    names_from = cluster,            # cluster 列的值变成列名
    values_from = density           # density 填入对应单元格
  )

# 查看效果
head(density_cluster_wide)

# 固定 cluster 顺序
cluster_levels <- c("1","2","3","4","5")
density_cluster_means_long$cluster <- factor(density_cluster_means_long$cluster, levels = cluster_levels)
density_cluster_means_long$celltype <- factor(density_cluster_means_long$celltype)

library(purrr)

density_pairwise_list <- density_cluster_wide %>%
  group_by(celltype) %>%
  group_split() %>%
  map(function(df){
    
    clusters <- colnames(df)[-(1:2)]
    
    res <- combn(clusters, 2, simplify = FALSE, FUN = function(cl){
      
      x <- df[[cl[1]]]
      y <- df[[cl[2]]]
      idx <- which(!is.na(x) & !is.na(y))
      
      if(length(idx) < 3) return(NULL)
      
      wt <- wilcox.test(x[idx], y[idx], paired = TRUE)
      
      data.frame(
        celltype = unique(df$celltype),
        group1 = cl[1],
        group2 = cl[2],
        n = length(idx),
        mean_group1 = mean(x[idx]),
        mean_group2 = mean(y[idx]),
        delta_mean = mean(x[idx] - y[idx]),
        p = wt$p.value
      )
    })
    
    bind_rows(res)
  }) %>% bind_rows() %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni"))


stat_df <- density_pairwise_list %>%
  filter(p.adj < 0.05) %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*"
    )
  )
stat_df <- stat_df %>%
  mutate(
    group1 = factor(group1, levels = cluster_levels),
    group2 = factor(group2, levels = cluster_levels),
    g1 = as.numeric(as.character(group1)),
    g2 = as.numeric(as.character(group2))
  ) %>%
  arrange(celltype, g1, g2)
stat_df <- stat_df %>%
  left_join(
    density_cluster_means_long %>%
      group_by(celltype) %>%
      summarise(ymax = max(density, na.rm = TRUE)),
    by = "celltype"
  ) %>%
  group_by(celltype) %>%
  mutate(y.position = ymax * (1 + 0.1 * row_number())) %>%
  ungroup()
cluster_cols <- c(
  "1" = "#F4FBD2",
  "2" = "#7CAAD0",
  "3" = "#DD4030",
  "4" = "#FEE79B",
  "5" = "#4575B4"
)
library(ggplot2)
library(ggpubr)

ggplot(density_cluster_means_long, aes(cluster, density, fill = cluster)) +
  
  geom_boxplot(outlier.shape = NA) +
  
  # facet_wrap 2行4列
  facet_wrap(~celltype, scales = "free_y", nrow = 3, ncol = 5) +
  
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  
  stat_pvalue_manual(
    stat_df,
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    label = "label",
    inherit.aes = FALSE,
    hide.ns = TRUE
  ) +
  
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
TMA_merged_NMF_score_tumor_dat <- cbind(TMA_merged_NMF_score_tumor_dat,density[rownames(TMA_merged_NMF_score_tumor_dat),c(16,17)])

program_cols <- paste0("program", LETTERS[1:8])

TMA_mean_base <- aggregate(
  TMA_merged_NMF_score_tumor_dat[, c(program_cols, "prop_Tumor.cells")],
  by = list(Sample = TMA_merged_NMF_score_tumor_dat$Sample),
  FUN = mean,
  na.rm = TRUE
)

library(dplyr)

# 找出所有细胞类型列
cell_cols <- setdiff(colnames(density), c("cell_id", "Sample"))

# 按 Sample 汇总均值
density_mean <- density %>%
  group_by(Sample) %>%
  summarise(across(all_of(cell_cols), ~ mean(.x, na.rm = TRUE)))
density_mean <- as.data.frame(density_mean)
rownames(density_mean) <- density_mean$Sample
rownames(TMA_mean_base) <- TMA_mean_base$Sample
density_nmf_mean <- cbind(TMA_mean_base[,c(-10)],density_mean[,c(-1)])

# 提取 program 列和细胞类型列
program_cols <- paste0("program", LETTERS[1:8])
celltype_cols <- setdiff(colnames(density_nmf_mean), c("Sample", program_cols))

# 初始化结果列表
cor_results <- list()

# 循环计算 Spearman 相关系数和 p 值
for (prog in program_cols) {
  for (cell in celltype_cols) {
    # 提取两列
    x <- density_nmf_mean[[prog]]
    y <- density_nmf_mean[[cell]]
    
    # 去掉 y == 0 的行，同时保持 x 对应行一致
    keep <- y != 0
    x_sub <- x[keep]
    y_sub <- y[keep]
    
    # 如果去掉之后没有数据，就跳过
    if (length(y_sub) < 2) next
    
    # 计算 Spearman 相关
    test <- cor.test(x_sub, y_sub, method = "spearman")
    
    # 保存结果
    cor_results[[length(cor_results) + 1]] <- data.frame(
      program = prog,
      celltype = cell,
      rho = test$estimate,
      p_value = test$p.value
    )
  }
}

# 合并所有结果
cor_results_df <- do.call(rbind, cor_results)

# 查看前几行
head(cor_results_df)

# 1️⃣ 添加显著性符号
cor_results_df <- cor_results_df %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))

# 2️⃣ 去掉 NA rho
heat_data <- cor_results_df %>% filter(!is.na(rho))
heat_data <- subset(heat_data,celltype != c("Liver.sinusoidal.endothelial.cells"))

# 3️⃣ 计算 rho 范围
rho_min <- min(heat_data$rho, na.rm = TRUE)
rho_max <- max(heat_data$rho, na.rm = TRUE)

# 4️⃣ 画热图
ggplot(heat_data, aes(x = celltype, y = program, fill = rho)) +
  geom_tile(color = "white") +                      # 热图格子
  geom_text(aes(label = sig), color = "black", size = 5) +  # 显著性符号
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(rho_min, rho_max),  # 负相关最蓝，正相关最红
    oob = scales::squish          # 超出范围自动压缩
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(fill = "Spearman rho")

plot_program_center_radial <- function(heat_data, program_name, celltypes){
  
  library(dplyr); library(ggplot2)
  
  # 过滤并按 celltypes 顺序排列
  df <- heat_data %>%
    filter(program == program_name, celltype %in% celltypes) %>%
    mutate(celltype = factor(celltype, levels = celltypes)) %>%
    arrange(celltype)
  
  n <- nrow(df)
  angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
  
  df <- df %>%
    mutate(
      angle = angles,
      rho_abs = abs(rho),
      line_len = (1 - rho)*2,
      x  = line_len * cos(angle),
      y  = line_len * sin(angle),
      cx = (line_len + 0.15) * cos(angle),
      cy = (line_len + 0.15) * sin(angle),
      sig_flag = sig != ""
    )
  
  ggplot(df) +
    
    # center
    geom_point(aes(x=0, y=0), size=40, shape=21, fill="grey85", stroke=0) +
    geom_text(aes(x=0, y=0, label=program_name), fontface="bold") +
    
    # lines
    geom_segment(aes(x=0, y=0, xend=x, yend=y), color="grey60") +
    
    # circles with legend
    geom_point(aes(x=cx, y=cy, fill=rho, size=rho_abs),
               shape=21, color=NULL, stroke=0) +
    
    # significance ring
    geom_point(data=df %>% filter(sig_flag),
               aes(x=cx, y=cy, size=rho_abs),
               shape=21, color="gold", fill=NA, stroke=2) +
    
    # labels
    geom_text(aes(x=cx*1.15, y=cy*1.15, label=celltype), size=3) +
    
    # color legend
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0,
                         limits=c(-1,1), name="Spearman rho") +
    
    # size legend
    scale_size_continuous(
      name = "|rho|",
      range = c(6, 30),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ) +
    
    coord_equal() +
    theme_void()
}

plot_program_center_radial(
  heat_data,
  program_name = "programA",
  celltypes = c("Tumor.cells","Fibroblasts","Smooth.muscle.cells","Pericytes","Endothelial","T.cells","B.cells","Plasma.cells","Mono.macro","DC"
                ,"Mast.cells")
)

plots <- list()
for (i in 1:length(program_cols)) {
  plots[[i]] <- plot_program_center_radial(
    heat_data,
    program_name = program_cols[i],
    celltypes = c("Tumor.cells","Fibroblasts","Smooth.muscle.cells","Pericytes","Endothelial","B.cells","Plasma.cells","T.cells","Mono.macro","DC"
                  ,"Mast.cells")
  )
}

p_all <- wrap_plots(plots, ncol = 4, nrow = 2)
p_all


