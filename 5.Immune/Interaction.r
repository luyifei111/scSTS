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

means <- read.delim("/cluster3/yflu/STS/cpdb/results_sample/statistical_analysis_means_11_12_2025_224351.txt",
                    header = TRUE, fill = TRUE, quote = "", comment.char = "", check.names = FALSE)
means_sig <- read.delim("/cluster3/yflu/STS/cpdb/results_sample/statistical_analysis_significant_means_11_12_2025_224351.txt",
                        header = TRUE, fill = TRUE, quote = "", comment.char = "", check.names = FALSE)
max(means_sig$`CH1_B cells|CH1_DC`)

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# 假设 means_sig 已经存在
# 找到第一个样本列的索引
sample_start_col <- 17

# 将长格式转换
means_long <- means_sig %>%
  pivot_longer(
    cols = sample_start_col:ncol(means_sig),
    names_to = "column_name",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%  # 去掉 NA
  mutate(
    sample = str_extract(column_name, "^[^_]+")  # 提取第一个 "_" 前的部分
  ) %>%
  select(sample, column_name, interacting_pair, value)

# 查看前几行结果
head(means_long)

means_long <- means_long %>%
  mutate(
    pair_clean = str_split_fixed(column_name, "\\|", 2),  # 先分成两列
    pair_clean = paste0(
      str_remove(pair_clean[,1], "^.*?_"), "|",           # 去掉左边部分的 "_" 前内容
      str_remove(pair_clean[,2], "^.*?_")                 # 去掉右边部分的 "_" 前内容
    )
  )

# 假设 means_sig 已经存在
# 找到第一个样本列的索引
sample_start_col <- 17

# 将长格式转换
means_long <- means_sig %>%
  pivot_longer(
    cols = sample_start_col:ncol(means_sig),
    names_to = "column_name",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%  # 去掉 NA
  mutate(
    sample = str_extract(column_name, "^[^_]+")  # 提取第一个 "_" 前的部分
  ) %>%
  select(sample, column_name, interacting_pair, value)

# 查看前几行结果
head(means_long)

means_long <- means_long %>%
  mutate(
    pair_clean = str_split_fixed(column_name, "\\|", 2),  # 先分成两列
    pair_clean = paste0(
      str_remove(pair_clean[,1], "^.*?_"), "|",           # 去掉左边部分的 "_" 前内容
      str_remove(pair_clean[,2], "^.*?_")                 # 去掉右边部分的 "_" 前内容
    )
  )
means_long$sample <- factor(means_long$sample)
means_long <- as.data.frame(means_long)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

# 每个样本的总 counts
sample_total_values <- means_long %>%
  group_by(sample) %>%
  summarise(total_counts = sum(value, na.rm = TRUE))

# 每个样本、每个 pair_clean 的 value 总和
pair_value_sum <- means_long %>%
  group_by(sample, pair_clean) %>%
  summarise(pair_value = sum(value, na.rm = TRUE), .groups = "drop")

# 转宽格式
pair_counts_wide <- pair_value_sum %>%
  pivot_wider(
    names_from = pair_clean,
    values_from = pair_value,
    values_fill = 0
  )

# 合并总和
pair_counts_wide <- left_join(sample_total_values, pair_counts_wide, by = "sample")

pair_counts_wide <- as.data.frame(pair_counts_wide)
rownames(pair_counts_wide) <- pair_counts_wide$sample

pair_counts_wide <- cbind(anno_sample_cluster_extended[rownames(pair_counts_wide),],pair_counts_wide)

pair_cols <- colnames(pair_counts_wide)[5:ncol(pair_counts_wide)]

# 对每列进行分组差异检验
pval_df <- map_dfr(pair_cols, function(colname) {
  df_sub <- pair_counts_wide[, c("Malignancy", colname)]
  names(df_sub) <- c("Malignancy", "value")
  
  # 若某列在某组全为0，检验可能报错，需 tryCatch 保护
  pval <- tryCatch({
    wilcox.test(value ~ Malignancy, data = df_sub)$p.value
  }, error = function(e) NA)
  
  data.frame(pair_clean = colname, p_value = pval)
})

pair_cols <- setdiff(colnames(pair_counts_wide), 
                     c("Disease", "Malignancy", "Louvain", "sample"))

# 计算差异倍数
fold_changes <- pair_counts_wide %>%
  group_by(Malignancy) %>%
  summarise(across(all_of(pair_cols), mean, na.rm = TRUE)) %>%
  pivot_longer(-Malignancy, names_to = "pair_clean", values_to = "mean_value") %>%
  pivot_wider(names_from = Malignancy, values_from = mean_value) %>%
  mutate(
    fold_change = Malignant / Benign,
    log2FC = log2(fold_change)
  )

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# 拆分 sender 和 receiver
pval_matrix <- pval_df[-1,] %>%
  separate(pair_clean, into = c("sender", "receiver"), sep = "\\|") %>%
  pivot_wider(
    names_from = receiver,
    values_from = p_value
  ) %>%
  column_to_rownames("sender")  # sender 作为行名

# 转成矩阵
pval_matrix <- as.matrix(pval_matrix)

# 假设 fold_changes 已包含 pair_clean 和 log2FC 列
fold_matrix <- fold_changes[-1,] %>%
  separate(pair_clean, into = c("sender", "receiver"), sep = "\\|") %>%
  select(sender, receiver, log2FC) %>%
  pivot_wider(
    names_from = receiver,
    values_from = log2FC
  ) %>%
  column_to_rownames("sender")  # sender 作为行名

# 转为矩阵
fold_matrix <- as.matrix(fold_matrix)

order <- rownames(pval_matrix)[c(11,4,12,8,3,1,9,10,6,2,5,7)]
pval_matrix <- pval_matrix[order,order]
fold_matrix <- fold_matrix[order,order]

# 创建显著性标注矩阵
sig_matrix <- pval_matrix %>%
  apply(., c(1,2), function(p) {
    if (is.na(p)) return("")
    else if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("")
  })

# 绘制热图
pheatmap(
  fold_matrix[order[1:11],order[1:11]],
  cluster_rows = F,
  cluster_cols = F,
  color = c(viridis(120)[1:40],viridis(50)[40:50]),
  display_numbers = sig_matrix[order[1:11],order[1:11]],  # 显著性标注
  fontsize_number = 20,
  number_color = "white",         # 星号颜色白色
  main = "log2FC of Malignant vs Benign"
)

coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", 
                     means_long$interacting_pair, value = FALSE)
means_long_coinhibitory <- means_long[coinhibitory, ]
means_long_coinhibitory <- means_long_coinhibitory %>%
  separate(pair_clean, into = c("sender", "receiver"), sep = "\\|", remove = FALSE)

means_long_coinhibitory <- subset(means_long_coinhibitory,receiver %in% unique(means_long_coinhibitory$receiver)[c(2,5,6,7,8,9)])
means_long_coinhibitory$interacting_pair <- as.factor(means_long_coinhibitory$interacting_pair)

library(dplyr)

coinhibitory_counts <- means_long_coinhibitory %>%
  group_by(sample, interacting_pair) %>%
  summarise(total_value = sum(value, na.rm = TRUE), .groups = "drop")

coinhibitory_counts_wide <- coinhibitory_counts %>%
  pivot_wider(
    names_from = interacting_pair,
    values_from = total_value,
    values_fill = 0
  )
coinhibitory_counts_wide <- as.data.frame(coinhibitory_counts_wide)
rownames(coinhibitory_counts_wide) <- coinhibitory_counts_wide$sample
coinhibitory_counts_wide <- cbind(anno_sample_cluster_extended[rownames(coinhibitory_counts_wide),],coinhibitory_counts_wide)

pair_cols <- colnames(coinhibitory_counts_wide)[5:ncol(coinhibitory_counts_wide)]

# 对每列进行分组差异检验
pval_df_coinhibitory <- map_dfr(pair_cols, function(colname) {
  df_sub <- coinhibitory_counts_wide[, c("Malignancy", colname)]
  names(df_sub) <- c("Malignancy", "value")
  
  # 若某列在某组全为0，检验可能报错，需 tryCatch 保护
  pval <- tryCatch({
    wilcox.test(value ~ Malignancy, data = df_sub)$p.value
  }, error = function(e) NA)
  
  data.frame(pair_clean = colname, p_value = pval)
})

fold_changes_coinhibitory <- coinhibitory_counts_wide %>%
  group_by(Malignancy) %>%
  summarise(across(all_of(pair_cols), mean, na.rm = TRUE)) %>%
  pivot_longer(-Malignancy, names_to = "pair_clean", values_to = "mean_value") %>%
  pivot_wider(names_from = Malignancy, values_from = mean_value) %>%
  mutate(
    fold_change = Malignant / Benign,
    log2FC = log2(fold_change)
  )
pair_order <- fold_changes_coinhibitory$pair_clean[c(7:9,13:14,12,16,11,1,2,4,5,10,3,6,15)]
fold_changes_coinhibitory <- as.data.frame(fold_changes_coinhibitory)
rownames(fold_changes_coinhibitory) <- fold_changes_coinhibitory$pair_clean
fold_changes_coinhibitory <- fold_changes_coinhibitory[pair_order,]

pval_df_coinhibitory <- as.data.frame(pval_df_coinhibitory)
rownames(pval_df_coinhibitory) <- pval_df_coinhibitory$pair_clean
pval_df_coinhibitory <- pval_df_coinhibitory[pair_order,]

coinhibitory_dot <- cbind(pval_df_coinhibitory,fold_changes_coinhibitory)
coinhibitory_dot$neg_log10_p <- -log10(coinhibitory_dot$p_value)
coinhibitory_dot <- coinhibitory_dot[,-1]

max_values <- coinhibitory_counts_wide %>%
  select(-Disease, -Malignancy, -Louvain, -sample) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "pair_clean", values_to = "max_value")

# 2️⃣ 将最大值表合并到 coinhibitory_dot
coinhibitory_dot_norm <- coinhibitory_dot %>%
  left_join(max_values, by = "pair_clean") %>%
  mutate(
    Benign_rel = Benign / max_value,
    Malignant_rel = Malignant / max_value
  )

# 3️⃣ 查看结果
coinhibitory_dot_norm %>%
  select(pair_clean, Benign, Malignant, max_value, Benign_rel, Malignant_rel)
rownames(coinhibitory_dot_norm) <- coinhibitory_dot_norm$pair_clean

coinhibitory_dot_norm$pair_clean <- factor(coinhibitory_dot_norm$pair_clean,levels = rev(rownames(coinhibitory_dot_norm)))

# pivot_longer
coinhibitory_dot_plot <- coinhibitory_dot_norm[-16,] %>%
  pivot_longer(cols = c("Benign_rel", "Malignant_rel", "log2FC"),
               names_to = "Variable", values_to = "Value") %>%
  mutate(highlight = p_value < 0.05)

coinhibitory_dot_plot <- coinhibitory_dot_plot %>%
  mutate(Variable = factor(Variable, levels = c("Malignant_rel", "Benign_rel", "log2FC")))

# 绘图
ggplot(coinhibitory_dot_plot, aes(x = Variable, y = pair_clean)) +
  # 基础点：颜色表示 Value，大小表示 neg_log10_p
  geom_point(aes(size = neg_log10_p, color = Value), shape = 16, alpha = 0.9) +
  # 红圈标记 p_value < 0.05，大小与底层点一致
  geom_point(
    data = subset(coinhibitory_dot_plot, highlight),
    aes(x = Variable, y = pair_clean, size = neg_log10_p),
    shape = 21, fill = NA, color = "red", stroke = 1.2
  ) +
  # 配色：反向 inferno
  scale_color_viridis_c(option = "viridis", direction = 1, name = "Value") +
  scale_size_continuous(name = "-log10(p)", range = c(2, 8)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Co-inhibitory interactions",
    subtitle = "Color: relative/log2FC value | Size: -log10(p) | Red circle: p-value < 0.05"
  )

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
sample_anno <- as.data.frame(table(rownames(anno_sample_cluster_extended),anno_sample_cluster_extended$Disease))
sample_anno <- subset(sample_anno,Freq > 0)
sample_anno <- sample_anno[,-3]
colnames(sample_anno) <- c("Samples","Diseases")
rownames(sample_anno) <- sample_anno$Samples

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]

p = pheatmap::pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
disease_order <- substr(rownames(aurocs_disease)[p$tree_row$order],5,nchar(rownames(aurocs_disease)[p$tree_row$order]))

disease_colors = c(
  "#8DD3C7","#97D7C5","#A1DBC3","#ACDFC1","#B6E3BF","#C0E7BD","#CBEBBC","#D5EFBA",
  "#FFFFB3","#F5F5B8","#FCCDE5","#F0D1E1","#E4D5DD","#B3DE69","#80B1D3","#FDB462",
  "#ADD8E6","#BC80BD","#BE8FBE","#C09EBF"
)

means_long_coinhibitory_box <- means_long_coinhibitory
colnames(means_long_coinhibitory_box)[1] <- "Samples"
means_long_coinhibitory_box_1 <- subset(means_long_coinhibitory_box,interacting_pair == "NECTIN2_TIGIT")

sample_levels <- levels(means_long_coinhibitory_box$Samples)

means_long_coinhibitory_sum <- means_long_coinhibitory_box_1 %>%
  group_by(Samples) %>%
  summarise(total_value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(Samples = sample_levels, fill = list(total_value = 0))
means_long_coinhibitory_sum <- as.data.frame(means_long_coinhibitory_sum)
rownames(means_long_coinhibitory_sum) <- means_long_coinhibitory_sum$Samples

means_long_coinhibitory_box_1 <- subset(means_long_coinhibitory_box,interacting_pair == "NECTIN3_TIGIT")

means_long_coinhibitory_sum_1 <- means_long_coinhibitory_box_1 %>%
  group_by(Samples) %>%
  summarise(total_value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(Samples = sample_levels, fill = list(total_value = 0))
means_long_coinhibitory_sum_1 <- as.data.frame(means_long_coinhibitory_sum_1)
rownames(means_long_coinhibitory_sum_1) <- means_long_coinhibitory_sum_1$Samples

means_long_coinhibitory_sum <- cbind(means_long_coinhibitory_sum,means_long_coinhibitory_sum_1)
means_long_coinhibitory_sum <- means_long_coinhibitory_sum[,-c(1,3)]
colnames(means_long_coinhibitory_sum) <- c("NECTIN2_TIGIT","NECTIN3_TIGIT")

p1 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum),
                                sample_anno, 
                                target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                                genes = c("NECTIN2_TIGIT"),
                                disease_order = disease_order,
                                disease_colors = disease_colors,
                                test="wilcox")
p2 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum),
                                sample_anno, 
                                target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                                genes = c("NECTIN3_TIGIT"),
                                disease_order = disease_order,
                                disease_colors = disease_colors,
                                test="wilcox")
p1+p2

sample_anno_1 <- cbind(sample_anno,anno_sample_cluster_extended[rownames(sample_anno),])
sample_anno_1 <- sample_anno_1[,-2]
sample_anno_1$Diseases <- as.character(sample_anno_1$Louvain)
p1 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum[rownames(subset(sample_anno_1,Malignancy == "Malignant"&Louvain %in% c("4","5"))),]),
                                subset(sample_anno_1,Malignancy == "Malignant"&Louvain %in% c("4","5")), 
                                target_diseases  = c("5"), 
                                genes = c("NECTIN2_TIGIT"),
                                test="wilcox",disease_order = c("4","5"),
                                disease_colors = disease_colors[c(14,15)])

p2 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum[rownames(subset(sample_anno_1,Malignancy == "Malignant"&Louvain %in% c("4","5"))),]),
                                subset(sample_anno_1,Malignancy == "Malignant"&Louvain %in% c("4","5")), 
                                target_diseases  = c("5"), 
                                genes = c("NECTIN3_TIGIT"),
                                test="wilcox",disease_order = c("4","5"),
                                disease_colors = disease_colors[c(14,15)])
p1+p2

means_long_coinhibitory_tigit <- subset(means_long_coinhibitory,interacting_pair %in% c("NECTIN2_TIGIT","NECTIN3_TIGIT"))

pair_counts <- as.data.frame(table(means_long_coinhibitory_tigit$pair_clean))
pair_counts <- pair_counts[base::order(pair_counts$Freq,decreasing = T),]

pair_counts <- pair_counts %>%
  separate(Var1, into = c("sender", "receiver"), sep = "\\|", remove = FALSE)

pair_counts$sender <- factor(pair_counts$sender,levels = rev(order))
pair_counts$receiver <- factor(pair_counts$receiver,levels = order[c(1:5,8,6,7,9:12)])

all_pairs <- expand.grid(
  sender = unique(pair_counts$sender),
  receiver = unique(pair_counts$receiver)
)

# 与原数据合并，缺失的组合 Freq = NA
pair_counts_full <- all_pairs %>%
  left_join(pair_counts, by = c("sender", "receiver"))

# 绘图
library(viridis)  # 如果还没加载

ggplot(pair_counts_full, aes(x = receiver, y = sender, fill = Freq)) +
  geom_tile(color = NA) +
  scale_fill_viridis_c(
    option = "inferno",
    begin = 0.2,    # 避开最暗的黑色
    end = 1,
    na.value = "grey",
    name = "Freq"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Receiver cell type",
    y = "Sender cell type",
    title = "Interaction frequency heatmap"
  )

means_long_coinhibitory_tigit_Tu_TNK <- subset(means_long_coinhibitory_tigit,pair_clean == "Tumor cells|T/NK cells")

means_long_coinhibitory_box_Tu_TNK <- means_long_coinhibitory_tigit_Tu_TNK
colnames(means_long_coinhibitory_box_Tu_TNK)[1] <- "Samples"
means_long_coinhibitory_box_Tu_TNK_1 <- subset(means_long_coinhibitory_box_Tu_TNK,interacting_pair == "NECTIN2_TIGIT")

sample_levels <- levels(means_long_coinhibitory_box_Tu_TNK$Samples)

means_long_coinhibitory_sum_Tu_TNK <- means_long_coinhibitory_box_Tu_TNK_1 %>%
  group_by(Samples) %>%
  summarise(total_value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(Samples = sample_levels, fill = list(total_value = 0))
means_long_coinhibitory_sum_Tu_TNK <- as.data.frame(means_long_coinhibitory_sum_Tu_TNK)
rownames(means_long_coinhibitory_sum_Tu_TNK) <- means_long_coinhibitory_sum_Tu_TNK$Samples

means_long_coinhibitory_box_Tu_TNK_1 <- subset(means_long_coinhibitory_box_Tu_TNK,interacting_pair == "NECTIN3_TIGIT")

means_long_coinhibitory_sum_Tu_TNK_1 <- means_long_coinhibitory_box_Tu_TNK_1 %>%
  group_by(Samples) %>%
  summarise(total_value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(Samples = sample_levels, fill = list(total_value = 0))
means_long_coinhibitory_sum_Tu_TNK_1 <- as.data.frame(means_long_coinhibitory_sum_Tu_TNK_1)
rownames(means_long_coinhibitory_sum_Tu_TNK_1) <- means_long_coinhibitory_sum_Tu_TNK_1$Samples

means_long_coinhibitory_sum_Tu_TNK <- cbind(means_long_coinhibitory_sum_Tu_TNK,means_long_coinhibitory_sum_Tu_TNK_1)
means_long_coinhibitory_sum_Tu_TNK <- means_long_coinhibitory_sum_Tu_TNK[,-c(1,3)]
colnames(means_long_coinhibitory_sum_Tu_TNK) <- c("NECTIN2_TIGIT","NECTIN3_TIGIT")

p1 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum_Tu_TNK),
                                sample_anno, 
                                target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                                genes = c("NECTIN2_TIGIT"),
                                disease_order = disease_order,
                                disease_colors = disease_colors,
                                test="wilcox")
p2 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum_Tu_TNK),
                                sample_anno, 
                                target_diseases  = disease_order[-c(1,4,6,7,9:13,18)], 
                                genes = c("NECTIN3_TIGIT"),
                                disease_order = disease_order,
                                disease_colors = disease_colors,
                                test="wilcox")
p1+p2

p1 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum_Tu_TNK[rownames(subset(sample_anno_1,Malignancy == "Malignant")),]),
                                subset(sample_anno_1,Malignancy == "Malignant"), 
                                target_diseases  = c("5"), 
                                genes = c("NECTIN2_TIGIT"),
                                test="wilcox",disease_order = c("1","2","3","4","5"),
                                disease_colors = disease_colors[c(1,9,11,14,15)])

p2 <- plot_gene_average_boxplot(t(means_long_coinhibitory_sum_Tu_TNK[rownames(subset(sample_anno_1,Malignancy == "Malignant")),]),
                                subset(sample_anno_1,Malignancy == "Malignant"), 
                                target_diseases  = c("5"), 
                                genes = c("NECTIN3_TIGIT"),
                                test="wilcox",disease_order = c("1","2","3","4","5"),
                                disease_colors = disease_colors[c(1,9,11,14,15)])
p1+p2

# 每个样本、每个 pair_clean 的 value 总和
pair_value_sum <- means_long_coinhibitory_tigit %>%
  group_by(sample, pair_clean) %>%
  summarise(pair_value = sum(value, na.rm = TRUE), .groups = "drop")

# 转宽格式
pair_counts_wide <- pair_value_sum %>%
  complete(
    sample = sample_levels,                                  # 按 sample 的 levels 补全
    pair_clean = unique(pair_value_sum$pair_clean),          # 补全所有 pair
    fill = list(pair_value = 0)
  ) %>%
  pivot_wider(
    names_from = pair_clean,
    values_from = pair_value,
    values_fill = 0
  )

# 合并总和
#pair_counts_wide <- left_join(sample_total_values, pair_counts_wide, by = "sample")

pair_counts_wide <- as.data.frame(pair_counts_wide)
rownames(pair_counts_wide) <- pair_counts_wide$sample

pair_counts_wide <- cbind(anno_sample_cluster_extended[rownames(pair_counts_wide),],pair_counts_wide)

pair_cols <- colnames(pair_counts_wide)[5:ncol(pair_counts_wide)]

# 对每列进行分组差异检验
pval_df <- map_dfr(pair_cols, function(colname) {
  df_sub <- pair_counts_wide[, c("Malignancy", colname)]
  names(df_sub) <- c("Malignancy", "value")
  
  # 若某列在某组全为0，检验可能报错，需 tryCatch 保护
  pval <- tryCatch({
    wilcox.test(value ~ Malignancy, data = df_sub)$p.value
  }, error = function(e) NA)
  
  data.frame(pair_clean = colname, p_value = pval)
})

# 计算差异倍数
fold_changes <- pair_counts_wide %>%
  group_by(Malignancy) %>%
  summarise(across(all_of(pair_cols), mean, na.rm = TRUE)) %>%
  pivot_longer(-Malignancy, names_to = "pair_clean", values_to = "mean_value") %>%
  pivot_wider(names_from = Malignancy, values_from = mean_value) %>%
  mutate(
    fold_change = Malignant - Benign,
    log2FC = log2(fold_change)
  )

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# 拆分 sender 和 receiver
pval_matrix <- pval_df %>%
  separate(pair_clean, into = c("sender", "receiver"), sep = "\\|") %>%
  pivot_wider(
    names_from = receiver,
    values_from = p_value
  ) %>%
  column_to_rownames("sender")  # sender 作为行名

# 转成矩阵
pval_matrix <- as.matrix(pval_matrix)

# 假设 fold_changes 已包含 pair_clean 和 log2FC 列
fold_matrix <- fold_changes %>%
  separate(pair_clean, into = c("sender", "receiver"), sep = "\\|") %>%
  select(sender, receiver, fold_change) %>%
  pivot_wider(
    names_from = receiver,
    values_from = fold_change
  ) %>%
  column_to_rownames("sender")  # sender 作为行名

# 转为矩阵
fold_matrix <- as.matrix(fold_matrix)

order <- unique(c(rownames(pval_matrix),colnames(pval_matrix)))[c(12,4,10,8,3,1,9,11,6,2,5,7)]
pval_matrix <- pval_matrix[intersect(order,rownames(pval_matrix)),intersect(order[c(1:5,8,6:7,9:12)],colnames(pval_matrix))]
fold_matrix <- fold_matrix[intersect(order,rownames(fold_matrix)),intersect(order[c(1:5,8,6:7,9:12)],colnames(fold_matrix))]

# 创建显著性标注矩阵
sig_matrix <- pval_matrix %>%
  apply(., c(1,2), function(p) {
    if (is.na(p)) return("")
    else if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("")
  })

fold_matrix_clean <- fold_matrix
fold_matrix_clean[is.infinite(fold_matrix_clean)] <- NA

# 绘制热图
pheatmap(
  fold_matrix_clean,
  cluster_rows = F,
  cluster_cols = F,
  color = c(viridis(40)[1:10],viridis(60)[15:60]),
  display_numbers = sig_matrix,  # 显著性标注
  fontsize_number = 20,
  number_color = "black",         # 星号颜色白色
  main = "Difference of Malignant vs Benign"
)

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
DimPlot(
  STS_pega_TNK_2,
  reduction = "umap",
  label = TRUE,
  raster = FALSE,
  group.by = "Celltype_new",
  cols = cols,
)
sample_info <- anno_sample_cluster_extended
# 按 Channel 匹配 Disease 和 Malignancy
STS_pega_TNK_2$Disease <- sample_info[STS_pega_TNK_2$Channel, "Disease"]
STS_pega_TNK_2$Malignancy <- sample_info[STS_pega_TNK_2$Channel, "Malignancy"]

table(STS_pega_TNK_2$Disease, useNA = "ifany")
table(STS_pega_TNK_2$Malignancy, useNA = "ifany")

FeaturePlot(STS_pega_TNK_2,features = c("TIGIT"),split.by = "Malignancy")
library(Nebulosa)

p_list <- lapply(unique(STS_pega_TNK_2$Malignancy), function(x) {
  cells_use <- rownames(STS_pega_TNK_2@meta.data)[STS_pega_TNK_2$Malignancy == x]
  
  p <- plot_density(
    STS_pega_TNK_2[, cells_use],
    features = "TIGIT",
    reduction = "umap"
  ) + ggtitle(x)
  
  # ---- 固定颜色 scale 到 0-0.055 ----
  p <- p +
    scale_color_gradientn(colours = viridis::turbo(256), limits = c(0, 0.055)) +
    scale_fill_gradientn(colours = viridis::turbo(256), limits = c(0, 0.055))
  
  p
})

wrap_plots(p_list, ncol = length(p_list))

p1 = VlnPlot(subset(STS_pega_TNK_2,Malignancy == "Benign"),group.by = "Celltype_new",cols = cols,features = "TIGIT",pt.size = 0)
p2 = VlnPlot(subset(STS_pega_TNK_2,Malignancy == "Malignant"),group.by = "Celltype_new",cols = cols,features = "TIGIT",pt.size = 0)
p1+p2

STS.integrated.pega <- readRDS("/cluster3/yflu/STS/pegasus/STS.integrated.pega_20240507.rds")

STS.integrated.pega$Channel <- as.character(STS.integrated.pega$Channel)
STS.integrated.pega$Disease <- sample_info[STS.integrated.pega$Channel, "Disease"]
STS.integrated.pega$Malignancy <- sample_info[STS.integrated.pega$Channel, "Malignancy"]

FeaturePlot(STS.integrated.pega,features = c("NECTIN2"),split.by = "Malignancy")

nectin_average <- AverageExpression(STS.integrated.pega,features = c("NECTIN2","NECTIN3"),group.by = "Channel")
tigit_average <- AverageExpression(subset(STS_pega_TNK_2,Celltype_new %in% c("TNFRSF+ Treg","ZNF683+ Tem","GZMK+ early Tem","GZMK+ Tem","GZMK+ Tex")),features = "TIGIT",group.by = "Channel")
nectin_average <- as.data.frame(nectin_average$RNA)
tigit_average <- as.data.frame(tigit_average$RNA)

tigit_average$T827 <- 0
nectin_tigit_average <- rbind(nectin_average,tigit_average[,colnames(nectin_average)])
nectin_tigit_average <- t(nectin_tigit_average)

colnames(nectin_tigit_average) <- c("NECTIN2","NECTIN3","TIGIT")
nectin_tigit_average <- cbind(nectin_tigit_average,means_long_coinhibitory_sum_Tu_TNK[rownames(nectin_tigit_average),])

res <- rcorr(as.matrix(nectin_tigit_average), type = "pearson")

VlnPlot(STS.integrated.pega,group.by = "Malignancy",features = "NECTIN3",pt.size = 0)
res <- rcorr(as.matrix(means_long_coinhibitory_sum_program), type = "pearson")

table_merged_proportion_raw <- readRDS("/cluster3/yflu/STS/cpdb/table_merged_proportion_raw.rds")

means_long_coinhibitory_sum_Tu_TNK$sum <- rowSums(means_long_coinhibitory_sum_Tu_TNK)
means_long_coinhibitory_sum_Tu_TNK_program <- cbind(means_long_coinhibitory_sum_Tu_TNK,table_merged_proportion_raw[rownames(means_long_coinhibitory_sum),c(41:48)])

res <- rcorr(as.matrix(means_long_coinhibitory_sum_Tu_TNK_program), type = "pearson")

cols <- c(
  colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
  colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
  colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
  colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
  colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
  colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
  "#ADD8E6",
  colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)]
)

full_names <- c("Hemangioma", "KHE", "Schwannoma", "MPNST", "Undifferentiated sarcoma",
                "RMS", "MRT", "IMT", "Angiosarcoma", "EWS/PNET",
                "NF", "Aggressive fibromatosis", "Liposarcoma", "Spindle cell tumor", "ASPS",
                "Infantile fibrosarcoma", "Synovial sarcoma", "Lipoblastoma", "Pecoma", "Lymphangioma")

# 对应缩写
abbreviations <- c("HE","KHE","SWN","MPNST","US","RMS","MRT","IMT","AS","EWS",
                   "NF","AF","LPS","SCT","ASPS","IFS","SS","LPB","PECOMA","LYM")

# 创建替换映射
name_map <- setNames(full_names,abbreviations)

disease_levels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS",
                    "NF","SWN","LYM","HE","KHE","MPNST","US",
                    "AS","RMS","PECOMA","EWS","MRT")
disease_levels <- name_map[disease_levels]
# 映射到命名向量
names(cols) <- disease_levels

means_long_coinhibitory_sum_Tu_TNK_program <- cbind(means_long_coinhibitory_sum_Tu_TNK_program,anno_sample_cluster_extended[rownames(means_long_coinhibitory_sum_Tu_TNK_program),])

means_long_coinhibitory_sum_Tu_TNK_program$Disease <- factor(
  means_long_coinhibitory_sum_Tu_TNK_program$Disease,
  levels = disease_levels
)

p <- ggscatter(
  means_long_coinhibitory_sum_Tu_TNK_program, 
  x = "sum", 
  y = "programD", 
  color = "Disease",
  size = 5,
  xlab = "NECTIN_TIGIT_sum", 
  ylab = "programD"
) +
  geom_smooth(
    method = "lm", 
    aes(x = sum, y = programD), 
    color = "black", 
    se = TRUE
  ) +
  stat_cor(
    aes(x = sum, y = programD), 
    method = "spearman"
  ) +
  scale_color_manual(values = cols) +
  theme_classic()
p

p <- ggscatter(
  means_long_coinhibitory_sum_Tu_TNK_program[rownames(subset(anno_sample_cluster_extended,Malignancy == "Benign")),], 
  x = "sum", 
  y = "programD", 
  color = "Disease",
  size = 5,
  xlab = "NECTIN_TIGIT_sum", 
  ylab = "programD"
) +
  geom_smooth(
    method = "lm", 
    aes(x = sum, y = programD), 
    color = "black", 
    se = TRUE
  ) +
  stat_cor(
    aes(x = sum, y = programD), 
    method = "spearman"
  ) +
  scale_color_manual(values = cols) +
  theme_classic()
p

f <- do.call(rbind, lapply(names(drug_hgnc_list_combined), function(nm){
  data.frame(
    drug = nm,
    gene = drug_hgnc_list_combined[[nm]],
    stringsAsFactors = FALSE
  )
}))

# 查看前几行
head(df)

# 保存为 CSV
write.csv(df, "drug_hgnc_table.csv", row.names = FALSE)

library(reshape2)
library(ggplot2)
library(dplyr)

# 提取矩阵
cor_mat <- res$r
p_mat   <- res$P

# 右上三角
cor_upper <- cor_mat
cor_upper[lower.tri(cor_upper, diag = FALSE)] <- NA

p_upper <- p_mat
p_upper[lower.tri(p_upper, diag = FALSE)] <- NA

# melt 相关矩阵
df <- melt(cor_upper, varnames = c("Var1", "Var2"))
# melt p 矩阵
df_p <- melt(p_upper, varnames = c("Var1", "Var2"))

# ⭐ 生成显著性星号
df$star <- sapply(df_p$value, function(p) {
  if (is.na(p)) return("")
  else if (p <= 0.001) return("***")
  else if (p <= 0.01)  return("**")
  else if (p <= 0.05)  return("*")
  else return("")
})

# 只显示上三角的星号
df$label <- ifelse(is.na(df$value), "", df$star)

# 恢复顺序
df$Var1 <- factor(df$Var1, levels = rownames(cor_mat))
df$Var2 <- factor(df$Var2, levels = colnames(cor_mat))

# 绘图
ggplot(subset(df,Var1 %in% c("NECTIN2_TIGIT","NECTIN3_TIGIT","sum")&Var2 %in% c("programA","programB","programC","programD","programE","programF","programG","programH")), aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  scale_fill_gradient2(
    low = "#3B9AB2", mid = "white", high = "#F21A00",
    midpoint = 0, limits = c(min(df$value), max(df$value))
  ) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "", fill = "Correlation",
       title = "Upper Triangle Correlation Heatmap (Significant Stars Only)")
