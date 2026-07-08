library(ggraph)
library(igraph)
library(tidygraph)
library(network)
library(ggnet)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
library(dplyr)
library(ggnewscale)

table_merged_proportion_raw <- readRDS("/cluster3/yflu/STS/cpdb/table_merged_proportion_raw.rds")
table_merged_proportion_adj <- scale(table_merged_proportion_raw)
table_merged_proportion_adj <- na.omit(t(table_merged_proportion_adj))
table_merged_proportion_adj <- as.matrix(t(table_merged_proportion_adj))

anno_immune <- readRDS("/cluster3/yflu/STS/microenvironmnet/anno_immune.rds")
subtype_group_map <- anno_immune$Group
names(subtype_group_map) <- anno_immune$celltype
group_colors <- RColorBrewer::brewer.pal(length(unique(anno_immune$Group)), "Set2")
names(group_colors) <- unique(anno_immune$Group)

# 读取一次注释
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

# 储存所有Weight
all_weights <- c()

# ======= 第一次循环：收集所有 Weight ========
Jaccard_df_list <- list()
subtype_list <- list()

for (i in 1:5) {
  louvain <- subset(anno_sample_cluster_extended, Louvain == i)
  louvain <- rownames(louvain)
  
  Freq_Jac <- table_merged_proportion_adj[louvain, 1:38]
  Freq_Jac <- round(Freq_Jac, 2)
  
  Jaccard_df <- data.frame()
  for (m in 1:37) {
    for (n in (m+1):38) {
      sub <- Freq_Jac[, c(m, n)]
      sub[sub >= 0.5] <- 1
      sub[sub < 0.5] <- 0
      intersection <- sum(rowSums(sub) == 2)
      union <- sum(rowSums(sub) >= 1)
      temp <- data.frame(Subtype1 = colnames(Freq_Jac)[m],
                         Subtype2 = colnames(Freq_Jac)[n],
                         Weight = round(intersection / union, 4))
      Jaccard_df <- rbind(Jaccard_df, temp)
    }
  }
  
  Jaccard_df <- na.omit(Jaccard_df)
  Jaccard_df <- subset(Jaccard_df, Weight > 0)
  
  Jaccard_df_list[[i]] <- Jaccard_df
  all_weights <- c(all_weights, Jaccard_df$Weight)
}

# 计算全局的最小和最大 Weight
weight_min <- min(all_weights)
weight_max <- max(all_weights)

# ======= 第二次循环：生成图并绘图 ========
plot_list <- list()

for (i in 1:5) {
  Jaccard_df <- Jaccard_df_list[[i]]
  Jaccard_df_sorted <- Jaccard_df[order(-Jaccard_df$Weight), ]
  
  selected_subtypes <- character()
  selected_rows <- data.frame()
  
  for (j in 1:nrow(Jaccard_df_sorted)) {
    row <- Jaccard_df_sorted[j, ]
    temp_subtypes <- unique(c(selected_subtypes, row$Subtype1, row$Subtype2))
    if (length(temp_subtypes) <= 20) {
      selected_subtypes <- temp_subtypes
      selected_rows <- rbind(selected_rows, row)
    }
    if (length(selected_subtypes) == 20) break
  }
  
  final_df <- subset(Jaccard_df,
                     Subtype1 %in% selected_subtypes & Subtype2 %in% selected_subtypes)
  
  # 构建 igraph 对象
  Jaccard_net <- graph_from_data_frame(final_df, directed = FALSE)
  
  # 添加 Group 属性
  V(Jaccard_net)$Subtype <- V(Jaccard_net)$name
  V(Jaccard_net)$Group <- subtype_group_map[V(Jaccard_net)$Subtype]
  
  # 绘图并保存到列表
  g <- ggraph(Jaccard_net, layout = "circle") +
    geom_edge_link(aes(width = Weight), color = "gray50") +
    geom_node_point(aes(color = Group), size = 15) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_edge_width(limits = c(weight_min, weight_max),
                     range = c(0.2, 2),
                     name = "Edge Weight") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_void() +
    ggtitle(paste("Louvain Cluster", i))
  
  plot_list[[i]] <- g
}

plot_list <- lapply(plot_list, function(p) {
  p + theme(plot.margin = margin(10, 10, 10, 10))  # 上、右、下、左
})

wrap_plots(plot_list, ncol = 5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

ann_colors = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
               colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
               colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])
names(ann_colors) <- labels
ann_colors <- list(Disease = ann_colors)
malignancy <- c("lightblue","pink")
names(malignancy) <- c("Benign","Malignant")
ann_colors$Malignancy <- malignancy

anno_sample_cluster_extended

df <- anno_sample_cluster_extended %>%
  select(Malignancy, Louvain, Disease) %>%
  mutate(across(everything(), as.factor)) %>%
  group_by(Malignancy, Louvain, Disease) %>%
  dplyr::summarise(Freq = n(), .groups = "drop")

# 分为两部分
df_left <- df %>%
  select(Malignancy, Louvain, Freq) %>%
  rename(axis1 = Malignancy, axis2 = Louvain)

df_right <- df %>%
  select(Louvain, Disease, Freq) %>%
  rename(axis1 = Louvain, axis2 = Disease)
df_right$axis2 <- factor(df_right$axis2, levels = names(ann_colors$Disease))

louvain_colors <- colorRampPalette(c("#d9f0d3", "#74c476"))(5)
names(louvain_colors) <- as.character(1:5)  # Louvain 是 factor/character

p1 <- ggplot(df_left, aes(axis1 = axis1, axis2 = axis2, y = Freq)) +
  geom_alluvium(aes(fill = axis1), width = 1/12, alpha = 0.6) +
  
  # stratum 分开着色：Malignancy 用 ann_colors，Louvain 用 green
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "black", size = 3) +
  
  scale_x_discrete(limits = c("Malignancy", "Louvain"), expand = c(.1, .05)) +
  scale_fill_manual(values = c(ann_colors$Malignancy, louvain_colors)) +
  theme_minimal() +
  labs(title = "Malignancy → Louvain", fill = "Stratum") +
  theme(legend.position = "bottom")

# 右图（Louvain → Disease）
p2 <- ggplot(df_right, aes(axis1 = axis1, axis2 = axis2, y = Freq)) +
  geom_alluvium(aes(fill = axis2), width = 1/12, alpha = 0.6) +
  
  # stratum 颜色：Louvain 用绿色，Disease 用自定义配色
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "black", size = 3) +
  
  scale_x_discrete(limits = c("Louvain", "Disease"), expand = c(.05, .1)) +
  scale_fill_manual(values = c(louvain_colors, ann_colors$Disease)) +
  theme_minimal() +
  labs(title = "Louvain → Disease", fill = "Stratum") +
  theme(legend.position = "bottom")

p1 + p2 + plot_layout(ncol = 2)

