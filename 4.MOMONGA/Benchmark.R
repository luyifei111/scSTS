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
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(reshape2) 
library(ggplot2)
library(ggpubr)
library(Seurat)
library(ggtree)
library(aplot)
library(stringr)

CNV_disease_durgs_sub <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_CNV_disease_durgs_sub_251203.rds")
aucell_druggable_genesets <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets.rds")

aucell_druggable_genesets <- t(aucell_druggable_genesets)

colnames(aucell_druggable_genesets) <- gsub("[|]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub(" ",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("-",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[(]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[)]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub(",",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[/]",".",  colnames(aucell_druggable_genesets))

samplenames <- unique(rownames(aucell_druggable_genesets))

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

CNV_disease_durgs_sub$Drugs <- as.character(CNV_disease_durgs_sub$Drugs)

## =========================
## 1. 提取目标药物
## =========================
target_drugs <- unique(CNV_disease_durgs_sub$Drugs)

## =========================
## 2. 构建 sample × drug 矩阵
## =========================
killing_subset_mat <- aucell_druggable_genesets[samplenames,target_drugs]

library(dplyr)
library(reshape2)
library(ggplot2)

## =========================
## 1. 构建 killing_long（sample × drug × score）并合并 subtype
## =========================
killing_long <- reshape2::melt(killing_subset_mat)
colnames(killing_long) <- c("Sample", "Drug", "Score")

killing_long <- merge(
  killing_long,
  anno_sample_cluster_extended,
  by.x = "Sample",
  by.y = "row.names"
)

## =========================
## 2. 先对 Score 做 min-max 归一化（0~1） per Drug
## =========================
killing_long <- killing_long %>%
  group_by(Drug) %>%
  mutate(
    Score_scaled = (Score - min(Score, na.rm = TRUE)) / (max(Score, na.rm = TRUE) - min(Score, na.rm = TRUE))
  ) %>%
  ungroup()

## =========================
## 3. 计算 subtype 内一致性（跳过单样本）
## =========================
drug_subtype_consistency <- killing_long %>%
  group_by(Disease, Drug) %>%
  summarise(
    mean_score = mean(Score_scaled, na.rm = TRUE),
    sd_score   = ifelse(n() > 1, sd(Score_scaled, na.rm = TRUE), NA),
    n_samples  = n(),
    .groups    = "drop"
  ) %>%
  mutate(
    cv = ifelse(!is.na(sd_score) & mean_score != 0, sd_score / mean_score, NA)
  ) %>%
  filter(n_samples > 1)

## =========================
## 4. 计算 overall 一致性
## =========================
drug_overall_consistency <- killing_long %>%
  group_by(Drug) %>%
  summarise(
    mean_score = mean(Score_scaled, na.rm = TRUE),
    sd_score   = ifelse(n() > 1, sd(Score_scaled, na.rm = TRUE), 0),
    cv         = ifelse(mean_score != 0, sd_score / mean_score, 0),
    n_samples  = n(),
    .groups    = "drop"
  )

## =========================
## 5. 计算 Shannon entropy
## =========================
drug_disease_mean <- killing_long %>%
  group_by(Drug, Disease) %>%
  summarise(mean_score = mean(Score_scaled, na.rm = TRUE), .groups = "drop")

drug_prob <- drug_disease_mean %>%
  group_by(Drug) %>%
  mutate(
    total = sum(mean_score),
    p = mean_score / (total + 1e-10)
  ) %>%
  ungroup()

drug_entropy <- drug_prob %>%
  group_by(Drug) %>%
  summarise(
    entropy = -sum(p * log(p + 1e-10)),
    .groups = "drop"
  ) %>%
  mutate(
    n_types = length(unique(drug_prob$Disease)),
    entropy_norm = entropy / log(n_types)
  )

## =========================
## 6. 整合为综合表格
## =========================
drug_summary <- drug_overall_consistency %>%
  select(Drug, mean_score, sd_score, cv) %>%
  left_join(drug_entropy %>% select(Drug, entropy_norm), by = "Drug") %>%
  rename(
    overall_mean_scaled = mean_score,
    overall_sd_scaled   = sd_score,
    overall_cv_scaled   = cv,
    shannon_entropy     = entropy_norm
  )

head(drug_summary)

library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

## =========================
## 1. 准备 subtype consistency 平均值
## =========================
drug_subtype_mean <- drug_subtype_consistency %>%
  group_by(Drug) %>%
  summarise(
    subtype_mean_score = mean(mean_score, na.rm = TRUE),
    subtype_sd        = mean(sd_score, na.rm = TRUE),
    .groups = "drop"
  )

## =========================
## 2. 整合指标
## =========================
drug_plot_df <- drug_summary %>%
  left_join(drug_subtype_mean, by = "Drug") %>%
  select(Drug, overall_mean_scaled, overall_sd_scaled, overall_cv_scaled, shannon_entropy, subtype_mean_score, subtype_sd)

# 转长格式
drug_plot_long <- drug_plot_df %>%
  pivot_longer(-Drug, names_to = "Metric", values_to = "Value")

## =========================
## 3. 设置顺序，让 Drug 保持原顺序
## =========================
drug_plot_long$Drug <- factor(drug_plot_long$Drug, levels = unique(drug_plot_df$Drug))

## =========================
## 4. 绘图
## =========================
ggplot(drug_plot_long, aes(x = 1, y = Drug, fill = Value)) +
  geom_tile(width = 0.8) +
  facet_wrap(~Metric, ncol = 1, scales = "free_x") +
  scale_fill_gradient(low = "red", high = "blue") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(fill = "Scaled Value")

# 假设已有数据
# drug_subtype_consistency, drug_overall_consistency, drug_entropy

library(dplyr)

## =========================
## 1. 对 subtype 一致性做 min-max
## =========================
drug_subtype_cv <- drug_subtype_consistency %>%
  filter(!is.na(cv)) %>%   # 跳过 NA
  group_by(Drug) %>%
  summarise(subtype_cv = mean(cv, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    subtype_cv_norm = (subtype_cv - min(subtype_cv)) / (max(subtype_cv) - min(subtype_cv))
  )

## =========================
## 2. 对 overall 一致性做 min-max
## =========================
drug_overall_cv <- drug_overall_consistency %>%
  filter(!is.na(cv)) %>%
  mutate(overall_cv_norm = (cv - min(cv)) / (max(cv) - min(cv))) %>%
  select(Drug, overall_cv_norm)

## =========================
## 3. 确保 drug_entropy 已有 entropy_norm
## =========================
# drug_entropy <- drug_entropy %>% mutate(entropy_norm = entropy_norm)  # 已经算过

## =========================
## 4. 合并三个指标
## =========================
drug_summary <- drug_subtype_cv %>%
  left_join(drug_overall_cv, by = "Drug") %>%
  left_join(drug_entropy %>% select(Drug, entropy_norm), by = "Drug")

## =========================
## 5. 计算复合指标
## =========================
# 因为 CV 和 entropy 越低越好，所以用 1 - norm
drug_summary <- drug_summary %>%
  mutate(
    composite_score = (1 - subtype_cv_norm + 1 - 1 + 1 - entropy_norm)/3
  ) %>%
  arrange(desc(composite_score))

## =========================
## 6. 查看 top 药物
## =========================
head(drug_summary)

library(ggplot2)
library(dplyr)

# 按 composite_score 排序
drug_summary_plot <- drug_summary

# 条形图
ggplot(drug_summary_plot, aes(x = Drug, y = composite_score, fill = composite_score)) +
  geom_col() +
  coord_flip() +  # 横向显示
  scale_fill_viridis_c(option = "plasma", direction = 1, limits = c(0,1)) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Composite Score (0-1)",
    fill = "Composite Score"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
