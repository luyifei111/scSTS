library(openxlsx)
library(dplyr)
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
library(irGSEA)
library(igraph)
library(stats)
library(spam)
library(muscat)
library(DESeq2)
library(apeglm)
library(irGSEA)
library(ggcorrplot)

results_1 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-1.xlsx","Sheet1")
results_2 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-2.xlsx","Sheet1")
results_3 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-3.xlsx","Sheet1")
results_4 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-4.xlsx","Sheet1")

colnames(results_2) <- colnames(results_1)
colnames(results_3) <- colnames(results_1)
colnames(results_4) <- colnames(results_1)

merged <- rbind(results_1,results_2,results_3,results_4)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

rownames(merged) <- merged$Image.Tag
merged_norm <- merged %>%
  mutate(across(-c(Image.Tag, Total.Cells), ~ . / Total.Cells))
merged_norm <- cbind(anno_sample_cluster_extended[rownames(merged_norm),],merged_norm)

plot_box_sig_global <- function(df, group_col, value_cols,
                                method_pair = "wilcox.test",
                                method_global = "kruskal.test",
                                p.adjust.method = "BH",
                                p.cutoff = 0.05,
                                point_color_col = NULL,
                                point_palette = NULL,
                                group_palette = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(rstatix)
  
  # =========================
  # 0️⃣ 检查列是否存在
  # =========================
  required_cols <- c(group_col, value_cols)
  if (!is.null(point_color_col)) {
    required_cols <- c(required_cols, point_color_col)
  }
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("缺少列:", paste(missing_cols, collapse = ", ")))
  }
  
  # =========================
  # 1️⃣ 转长表
  # =========================
  df_long <- df %>%
    dplyr::select(all_of(required_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(value_cols),
      names_to = "variable",
      values_to = "value"
    )
  
  df_long[[group_col]] <- as.factor(df_long[[group_col]])
  
  if (!is.null(point_color_col)) {
    df_long[[point_color_col]] <- as.factor(df_long[[point_color_col]])
  }
  
  plots <- list()
  
  # =========================
  # 2️⃣ 循环每个变量
  # =========================
  for (var in unique(df_long$variable)) {
    
    df_sub <- df_long %>%
      dplyr::filter(variable == var) %>%
      dplyr::filter(!is.na(value))
    
    # =========================
    # 2.1️⃣ 全局检验
    # =========================
    global_p <- tryCatch({
      df_sub %>%
        rstatix::kruskal_test(
          formula = as.formula(paste("value ~", group_col))
        ) %>%
        dplyr::pull(p)
    }, error = function(e) NA)
    
    global_label <- ifelse(
      is.na(global_p),
      "kruskal-wallis p = NA",
      paste0("kruskal-wallis p = ", signif(global_p, 3))
    )
    
    # =========================
    # 2.2️⃣ 作图
    # =========================
    if (is.null(point_color_col)) {
      
      p <- ggplot(df_sub,
                  aes_string(x = group_col, y = "value", fill = group_col)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 3, alpha = 0.7)
      
    } else {
      
      p <- ggplot(df_sub,
                  aes_string(x = group_col, y = "value", fill = group_col)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes_string(color = point_color_col),
                    width = 0.2, size = 3, alpha = 0.7)
    }
    
    p <- p +
      labs(
        title = var,
        subtitle = global_label,
        x = group_col,
        y = "Proportion"
      ) +
      theme_bw() +
      theme(
        legend.position = ifelse(is.null(point_color_col), "none", "right"),
        plot.title = element_text(hjust = 0.5)
      )
    
    # =========================
    # 🎨 配色
    # =========================
    if (!is.null(group_palette)) {
      p <- p + scale_fill_manual(values = group_palette)
    }
    
    if (!is.null(point_palette) && !is.null(point_color_col)) {
      p <- p + scale_color_manual(values = point_palette)
    }
    
    # =========================
    # 2.3️⃣ 两两检验
    # =========================
    stat_res <- tryCatch({
      df_sub %>%
        rstatix::pairwise_wilcox_test(
          formula = as.formula(paste("value ~", group_col)),
          p.adjust.method = p.adjust.method
        ) %>%
        dplyr::filter(p.adj < p.cutoff)
    }, error = function(e) NULL)
    
    # =========================
    # 2.4️⃣ 显著性标注
    # =========================
    if (!is.null(stat_res) && nrow(stat_res) > 0) {
      
      y_max <- max(df_sub$value, na.rm = TRUE)
      y_min <- min(df_sub$value, na.rm = TRUE)
      range_y <- y_max - y_min
      
      if (range_y == 0) range_y <- abs(y_max)
      if (range_y == 0) range_y <- 1
      
      step <- range_y * 0.15
      
      stat_res$y.position <- y_max + seq_len(nrow(stat_res)) * step
      
      p <- p +
        stat_pvalue_manual(
          stat_res,
          label = "p.adj.signif",
          xmin = "group1",
          xmax = "group2",
          y.position = "y.position",
          inherit.aes = FALSE
        )
    }
    
    plots[[var]] <- p
  }
  
  return(plots)
}

ann_colors = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
               colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
               colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
               "#ADD8E6",
               colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])
names(ann_colors) <- labels

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

merged_norm$Louvain <- as.character(merged_norm$Louvain)

plots <- plot_box_sig_global(
  df = merged_norm,
  group_col = "Louvain",
  point_color_col = "Disease",
  point_palette = ann_colors,
  group_palette = c("#edf8fb", "#c2e8e8", "#8dd2c3","#57ba93","#2ca25f"),
  value_cols = c("Immune.cells","B.Cells","T.Cells","Mono.macro","Neutrophil","NK","CD4..T","CD8..T","Endo","Schwann","Fibroblasts")
)
plots

close_1 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-1_close.xlsx","Summary")
close_2 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-2_close.xlsx","Summary")
close_3 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-3_close.xlsx","Summary")
close_4 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-4_close.xlsx","Sheet1")

merge_close <- rbind(close_1[,-c(2,3)],close_2[,-c(2,3)],close_3,close_4)
rownames(merge_close) <- merge_close$Image.Tag
merge_close <- cbind(anno_sample_cluster_extended[merge_close$Image.Tag,],merge_close)

neighbor_1 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-1_neighbor.xlsx","Sheet1")
neighbor_2 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-2_neighbor.xlsx","Summary")
neighbor_3 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-3_neighbor.xlsx","Summary")
neighbor_4 <- read.xlsx("/cluster3/yflu/STS/mIHC/results/4-4_neighbor.xlsx","Sheet1")
merge_neighbor <- rbind(neighbor_1,neighbor_2,neighbor_3,neighbor_4)

rownames(merge_neighbor) <- merge_neighbor$Image.Tag
merge_neighbor <- cbind(anno_sample_cluster_extended[merge_neighbor$Image.Tag,],merge_neighbor)

merge_neighbor_norm <- merge_neighbor %>%
  dplyr::mutate(
    dplyr::across(
      c(
        Target.Within.0.50µm.of.Reference,
        Target.Within.50.100µm.of.Reference,
        Target.Within.100.150µm.of.Reference,
        Target.Within.150.200µm.of.Reference
      ),
      ~ .x / Target.Cells.Within.200µm.of.Reference
    )
  )

dist_centers <- c(25, 75, 125, 175)

library(dplyr)

# 距离中心点（单位 µm）
dist_centers <- c(25, 75, 125, 175)

merge_neighbor_log_slope <- merge_neighbor %>%
  
  # 先转为比例（如果你已经有 ratio 可以删掉这一步）
  dplyr::mutate(
    dplyr::across(
      c(
        Target.Within.0.50µm.of.Reference,
        Target.Within.50.100µm.of.Reference,
        Target.Within.100.150µm.of.Reference,
        Target.Within.150.200µm.of.Reference
      ),
      ~ .x / Target.Cells.Within.200µm.of.Reference
    )
  ) %>%
  
  dplyr::rowwise() %>%
  dplyr::mutate(
    log_slope = {
      y <- c(
        Target.Within.0.50µm.of.Reference,
        Target.Within.50.100µm.of.Reference,
        Target.Within.100.150µm.of.Reference,
        Target.Within.150.200µm.of.Reference
      )
      
      # 👉 处理0值（关键！）
      y <- y + 1e-6
      
      # 👉 如果全是NA或异常
      if(all(is.na(y))) {
        NA
      } else {
        coef(lm(log(y) ~ dist_centers))[2]
      }
    }
  ) %>%
  dplyr::ungroup()

library(dplyr)
library(ggplot2)
library(emmeans)
library(broom)
library(ggpubr)
library(patchwork)

# =========================
# 1. 数据 + reference
# =========================
df <- merge_neighbor_log_slope %>%
  dplyr::mutate(Louvain = factor(Louvain, levels = c("1","2","3","4","5")))

df$Louvain <- relevel(df$Louvain, ref = "2")

# =========================
# 2. weighted lm
# =========================
fit <- lm(
  log_slope ~ Louvain,
  data = df,
  weights = log1p(Target.Cells.Within.200µm.of.Reference)
)

anova_res <- anova(fit)
anova_p <- signif(anova_res$`Pr(>F)`[1], 3)

emm <- emmeans(fit, ~ Louvain)
emm_df <- as.data.frame(emm)

# =========================
# 3. forest plot数据
# =========================
coef_df <- broom::tidy(fit, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    group = gsub("Louvain", "", term)
  )

coef_df$group <- factor(coef_df$group, levels = c("1","3","4","5"))

# =========================
# 4. 显著性标记函数
# =========================
coef_df$signif <- case_when(
  coef_df$p.value < 0.001 ~ "***",
  coef_df$p.value < 0.01 ~ "**",
  coef_df$p.value < 0.05 ~ "*",
  TRUE ~ "ns"
)

df$Louvain <- factor(df$Louvain, levels = c("1","2","3","4","5"))

# =========================
# 5. boxplot（带ANOVA）
# =========================
p_box <- ggplot(df, aes(x = Louvain, y = log_slope)) +
  
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  
  geom_jitter(
    aes(size = Target.Cells.Within.200µm.of.Reference),
    width = 0.2, alpha = 0.6
  ) +
  
  geom_point(
    data = emm_df,
    aes(x = Louvain, y = emmean),
    color = "red",
    size = 3
  ) +
  
  # =========================
# ⭐ ANOVA标注（核心）
# =========================
annotate(
  "text",
  x = Inf,
  y = Inf,
  label = paste0("ANOVA p = ", anova_p),
  hjust = 1.1,
  vjust = 1.5,
  size = 4
) +
  
  theme_classic() +
  
  labs(
    title = "log_slope distribution",
    x = "Louvain",
    y = "log_slope"
  )

# =========================
# 6. forest plot（带显著性）
# =========================
p_forest <- ggplot(coef_df, aes(x = group, y = estimate)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.2
  ) +
  
  geom_text(
    aes(label = signif, y = conf.high + 0.005),
    size = 5
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    title = "Effect vs Louvain2 (weighted lm)",
    x = NULL,
    y = "β (log_slope change)"
  )

# =========================
# 7. 拼图
# =========================
p_all <- p_box / p_forest +
  plot_layout(heights = c(2, 1))
p_all

