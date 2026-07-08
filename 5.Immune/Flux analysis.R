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
library(tidyr)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
samplenames <- rownames(anno_sample_cluster_extended)

scFEA_path <- paste("/cluster3/yflu/STS/scFEA/output/",samplenames,"_flux.csv",sep = "")

celltype <- readRDS("/cluster3/yflu/STS/scFEA/celltype.rds")

i=1
scFEA <- read.csv(scFEA_path[i])

rownames(scFEA) <- scFEA$X
scFEA <- scFEA[,-1]
rownames(scFEA) <- gsub("\\.", "-", rownames(scFEA))

scFEA_celltype <- celltype[rownames(scFEA)]
df <- as.data.frame(scFEA)
df$celltype <- scFEA_celltype[rownames(df)]

celltype_levels <- levels(df$celltype)

scFEA_celltype_mean <- df %>%
  group_by(celltype) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            .groups = "drop") %>%
  complete(celltype = factor(celltype_levels, levels = celltype_levels))

scFEA_celltype_mean_list <- list()
scFEA_celltype_mean_list[[i]] <- scFEA_celltype_mean
names(scFEA_celltype_mean_list)[i] <- samplenames[i]

for (i in 2:length(samplenames)) {
  scFEA <- read.csv(scFEA_path[i])
  
  rownames(scFEA) <- scFEA$X
  scFEA <- scFEA[,-1]
  rownames(scFEA) <- gsub("\\.", "-", rownames(scFEA))
  
  scFEA_celltype <- celltype[rownames(scFEA)]
  df <- as.data.frame(scFEA)
  df$celltype <- scFEA_celltype[rownames(df)]
  
  celltype_levels <- levels(df$celltype)
  
  scFEA_celltype_mean <- df %>%
    group_by(celltype) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE),
              .groups = "drop") %>%
    complete(celltype = factor(celltype_levels, levels = celltype_levels))
  
  scFEA_celltype_mean_list[[i]] <- scFEA_celltype_mean
  names(scFEA_celltype_mean_list)[i] <- samplenames[i]
}

library(dplyr)
library(tidyr)
library(purrr)

# 合并 scFEA 数据
scFEA_all <- map2_df(scFEA_celltype_mean_list, names(scFEA_celltype_mean_list), 
                     ~ mutate(.x, Sample = .y))

# 加入 Louvain 信息
scFEA_all <- scFEA_all %>%
  left_join(anno_sample_cluster_extended %>% tibble::rownames_to_column("Sample"), by = "Sample")

# pivot_longer
scFEA_long <- scFEA_all %>%
  pivot_longer(cols = starts_with("M_"), names_to = "Module", values_to = "flux")

library(dplyr)

# 去掉 flux NA
scFEA_clean <- scFEA_long %>% filter(!is.na(flux))

celltypes <- unique(scFEA_clean$celltype)
louvains  <- unique(scFEA_clean$Louvain)
modules   <- unique(scFEA_clean$Module)

result_list <- list()

for(ct in celltypes){
  df_ct <- scFEA_clean %>% filter(celltype == ct)
  
  for(mod in modules){
    df_mod <- df_ct %>% filter(Module == mod)
    
    for(lv in louvains){
      flux_group <- df_mod %>% filter(Louvain == lv) %>% pull(flux)
      flux_other <- df_mod %>% filter(Louvain != lv) %>% pull(flux)
      
      # 如果任意一组为空就跳过
      if(length(flux_group) == 0 || length(flux_other) == 0){
        mean_target <- NA
        mean_other  <- NA
        median_diff <- NA
        wilcox_p   <- NA
      } else {
        mean_target <- mean(flux_group)
        mean_other  <- mean(flux_other)
        median_diff <- median(flux_group) - median(flux_other)
        mean_diff <-  mean(flux_group) - mean(flux_other)
        wilcox_p    <- wilcox.test(flux_group, flux_other, exact = FALSE)$p.value
      }
      
      result_list[[length(result_list)+1]] <- data.frame(
        celltype    = ct,
        Module      = mod,
        Louvain     = lv,
        mean_target = mean_target,
        mean_other  = mean_other,
        median_diff = median_diff,
        mean_diff = mean_diff,
        wilcox_p    = wilcox_p
      )
    }
  }
  print(ct)
}


# 合并为 data.frame
result_df <- bind_rows(result_list)

# 按 celltype + Module 做 FDR 校正
result_df <- result_df %>%
  group_by(celltype, Module) %>%
  mutate(wilcox_padj = p.adjust(wilcox_p, method = "BH")) %>%
  ungroup()

head(result_df)
result_df_sig <- subset(result_df,wilcox_padj < 0.05)

library(dplyr)
library(tidyr)

df_tc <- result_df_sig %>%
  filter(celltype == "T cells", Louvain %in% c(2,4)) %>%
  select(Module, Louvain, mean_diff,mean_target, mean_other) %>%
  pivot_wider(names_from = Louvain, values_from = c(mean_diff,mean_target, mean_other), names_prefix = "L")

# 找方向相反的模块
df_tc_opposite <- df_tc %>%
  filter(!is.na(mean_diff_L2) & !is.na(mean_diff_L4)) %>%
  filter(sign(mean_diff_L2) != sign(mean_diff_L4))

df_tc_opposite

library(dplyr)
library(tidyr)

# 先筛选 T cells
df_tc <- result_df_sig %>% filter(celltype == "T cells")

`%notin%` <- Negate(`%in%`)

# 遍历每个 Module
modules_only_2_4 <- df_tc %>%
  group_by(Module) %>%
  # 确保 L2 和 L4 都存在
  filter(all(c(2,4) %in% Louvain)) %>%
  # L2 和 L4 都显著
  filter(all(Louvain %in% c(2,4) & wilcox_padj < 0.05)) %>%
  # 方向相反
  filter(sign(mean_diff[Louvain == 2]) != sign(mean_diff[Louvain == 4])) %>%
  # 其他 Louvain 不显著
  #filter(all(!(Louvain %notin% c(2,4) & wilcox_padj < 0.05))) %>%
  ungroup() %>%
  pull(Module) %>% unique()

df_tc_opposite_only <- df_tc %>%
  filter(Module %in% modules_only_2_4, Louvain %in% c(2,4)) %>%
  select(Module, Louvain, mean_diff, mean_target, mean_other) %>%
  pivot_wider(
    names_from = Louvain,
    values_from = c(mean_diff, mean_target, mean_other),
    names_prefix = "L"
  ) %>%
  arrange(Module)

df_m <- result_df_sig %>% filter(celltype == "Monocytes/macrophages")

modules_only_3_4 <- df_m %>%
  group_by(Module) %>%
  # 确保 L3 和 L4 都存在
  filter(all(c(3,4) %in% Louvain)) %>%
  # L2 和 L4 都显著
  filter(all(Louvain %in% c(3,4) & wilcox_padj < 0.05)) %>%
  # 方向相反
  filter(sign(mean_diff[Louvain == 3]) != sign(mean_diff[Louvain == 4])) %>%
  # 其他 Louvain 不显著
  filter(all(!(Louvain %notin% c(3,4) & wilcox_padj < 0.05))) %>%
  ungroup() %>%
  pull(Module) %>% unique()

# 构建最终表格
df_tm_opposite_only <- df_m %>%
  filter(Module %in% modules_only_3_4, Louvain %in% c(3,4)) %>%
  select(Module, Louvain, mean_diff, mean_target, mean_other) %>%
  pivot_wider(
    names_from = Louvain,
    values_from = c(mean_diff, mean_target, mean_other),
    names_prefix = "L"
  ) %>%
  arrange(Module)

result_df_sig_4_tumor <- subset(result_df_sig,Louvain == 4&celltype == 'Tumor')
library(dplyr)
library(tidyr)

df_plot <- result_df_sig_4_tumor %>%
  select(Module, mean_target, mean_other, median_diff, wilcox_padj) %>%
  pivot_longer(cols = c(mean_target, mean_other, median_diff),
               names_to = "metric", values_to = "value")
df_plot <- df_plot %>%
  mutate(sig = case_when(
    wilcox_padj < 0.001 ~ "***",
    wilcox_padj < 0.01  ~ "**",
    wilcox_padj < 0.05  ~ "*",
    TRUE ~ ""
  ))
library(ggplot2)

df_plot <- df_plot %>%
  mutate(Module_num = as.numeric(gsub("M_", "", Module))) %>%
  arrange(Module_num) %>%
  mutate(Module = factor(Module, levels = rev(unique(Module))))

ggplot(df_plot, aes(x = metric, y = Module, fill = value)) +
  geom_tile(color = "grey80") +
  
  # p value符号只标在median_diff列
  geom_text(data = df_plot %>% filter(metric == "median_diff"),
            aes(label = sig), color = "black", size = 5) +
  
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "", fill = "Flux")

library(dplyr)
library(tidyr)

filter_opposite_mean_diff_multi <- function(df, celltype_input, group1, group2) {
  
  df %>%
    filter(celltype == celltype_input) %>%
    # 计算每组 Louvain 的 mean_diff 平均值（可以理解为组内方向）
    mutate(group = case_when(
      Louvain %in% group1 ~ "group1",
      Louvain %in% group2 ~ "group2",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    group_by(Module, group) %>%
    summarise(mean_diff_group = mean(mean_diff), .groups = "drop") %>%
    # 拉宽，每个 Module 两组并排
    pivot_wider(names_from = group, values_from = mean_diff_group) %>%
    # 筛选两组都有数据且符号相反
    filter(!is.na(group1), !is.na(group2)) %>%
    filter(group1 * group2 < 0)
}

# 用法：
result <- filter_opposite_mean_diff_multi(result_df_sig, "T cells", c(2:3), c(4))

M168 <- read.csv("/cluster3/yflu/STS/scFEA-master/data/Human_M168_information.symbols.csv")

result_df_sig_4_T <- subset(result_df_sig,celltype == 'T cells'&Louvain == 4)
df_opposite <- result_df_sig_4_T %>%
  dplyr::filter(
    !is.na(mean_target),
    !is.na(mean_other),
    sign(mean_target) != sign(mean_other)
  )

library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(ggraph)
library(grid)

library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(ggraph)
library(grid)

plot_scFEA_network <- function(result_df, module_anno, super_id_vec,
                               show_reaction_node = TRUE,
                               seed = 123,
                               title = NULL,
                               sig_alpha = 0.2) {   # 新增 title 和 sig_alpha 参数
  # ================= 1. 筛选模块 =================
  df <- result_df %>%
    left_join(module_anno, by = c("Module" = "X")) %>%
    filter(Supermodule_id %in% super_id_vec) %>%
    filter(!is.na(Compound_IN_name)) %>%
    filter(abs(mean_diff) > 0)
  
  # ================= 2. 拆分多代谢物 =================
  df_long <- df %>%
    mutate(
      Compound_IN_name  = str_split(Compound_IN_name, "\\+"),
      Compound_OUT_name = str_split(Compound_OUT_name, "\\+")
    ) %>%
    unnest(Compound_IN_name) %>%
    unnest(Compound_OUT_name) %>%
    mutate(
      Compound_IN_name  = str_trim(Compound_IN_name),
      Compound_OUT_name = str_trim(Compound_OUT_name)
    )
  
  # ================= 3. 构建 edges =================
  if(show_reaction_node){
    e1 <- df %>%
      mutate(reaction = Module,
             Compound_IN_name = str_split(Compound_IN_name, "\\+")) %>%
      unnest(Compound_IN_name) %>%
      mutate(Compound_IN_name = str_trim(Compound_IN_name)) %>%
      transmute(from = Compound_IN_name, to = reaction,
                mean_diff, mean_target, wilcox_p)
    
    e2 <- df %>%
      mutate(reaction = Module,
             Compound_OUT_name = str_split(Compound_OUT_name, "\\+")) %>%
      unnest(Compound_OUT_name) %>%
      mutate(Compound_OUT_name = str_trim(Compound_OUT_name)) %>%
      transmute(from = reaction, to = Compound_OUT_name,
                mean_diff, mean_target, wilcox_p)
    
    edges <- bind_rows(e1, e2)
  } else {
    edges <- df_long %>%
      transmute(from = Compound_IN_name, to = Compound_OUT_name,
                mean_diff, mean_target, wilcox_p)
  }
  
  # ================= 4. 补充绘图属性 =================
  edges <- edges %>%
    mutate(
      weight = abs(mean_diff),              # 控制线宽
      length_value = abs(mean_target),      # 控制布局长度
      sign = ifelse(mean_diff > 0, "pos", "neg"),
      line_type = ifelse(mean_target > 0, "solid", "dashed"),
      sig = wilcox_p < 0.05
    )
  
  # ================= 5. 构建 graph =================
  nodes <- tibble(name = unique(c(edges$from, edges$to)))
  g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  
  # ================= 6. 生成布局，让边长度反映 mean_target =================
  edge_lengths <- edges$length_value
  edge_lengths[edge_lengths == 0] <- min(edge_lengths[edge_lengths > 0])/10  # 防止0
  set.seed(seed)  # 保证可复现
  layout <- layout_with_fr(g, weights = 1/edge_lengths)
  
  # ================= 7. 绘图 =================
  p <- ggraph(g, layout = layout) +
    geom_edge_link(aes(
      edge_width = weight,
      edge_color = sign,
      linetype = line_type,
      alpha = sig
    ),
    arrow = arrow(length = unit(2, "mm")),
    end_cap = circle(2, "mm")) +
    geom_node_point(size = 4, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c("pos" = "red", "neg" = "#7EAEDB")) +
    scale_edge_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
    scale_edge_alpha_manual(values = c("TRUE" = 1, "FALSE" = sig_alpha)) +  # 显著性透明度
    theme_void() +
    labs(title = ifelse(is.null(title),
                        paste0("Putrescine → Spermine network for supermodules: ",
                               paste(super_id_vec, collapse = ",")),
                        title))
  
  return(p)
}

result_df_4_Tumor <- subset(result_df,celltype == 'Tumor'&Louvain == 4)

plot_scFEA_network(
  result_df = result_df_4_Tumor,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 123,
  title = "Louvain 4 Tumor cells"
)

result_df_4_T <- subset(result_df,celltype == 'T cells'&Louvain == 4)
plot_scFEA_network(
  result_df = result_df_4_T,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 119,
  title = "Louvain 4 T cells"
)

result_df_5_T <- subset(result_df,celltype == 'T cells'&Louvain == 5)
plot_scFEA_network(
  result_df = result_df_5_T,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 140,
  title = "Louvain 5 T cells"
)

result_df_2_T <- subset(result_df,celltype == 'T cells'&Louvain == 2)
plot_scFEA_network(
  result_df = result_df_2_T,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 119,
  title = "Louvain 2 T cells"
)

result_df_4_B <- subset(result_df,celltype == 'B cells'&Louvain == 4)
plot_scFEA_network(
  result_df = result_df_4_B,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 119,
  title = "Louvain 4 B cells"
)

result_df_5_B <- subset(result_df,celltype == 'B cells'&Louvain == 5)
plot_scFEA_network(
  result_df = result_df_5_B,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 126,
  title = "Louvain 5 B cells"
)

result_df_2_B <- subset(result_df,celltype == 'B cells'&Louvain == 2)
plot_scFEA_network(
  result_df = result_df_2_B,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 119,
  title = "Louvain 2 B cells"
)

result_df_4_M <- subset(result_df,celltype == 'Monocytes/macrophages'&Louvain == 4)
plot_scFEA_network(
  result_df = result_df_4_M,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 123,
  title = "Louvain 4 Monocytes/macrophages"
)

result_df_5_M <- subset(result_df,celltype == 'Monocytes/macrophages'&Louvain == 5)
plot_scFEA_network(
  result_df = result_df_5_M,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 123,
  title = "Louvain 5 Monocytes/macrophages"
)

result_df_3_M <- subset(result_df,celltype == 'Monocytes/macrophages'&Louvain == 3)
plot_scFEA_network(
  result_df = result_df_3_M,
  module_anno = M168,
  super_id = c(10,11,6),
  show_reaction_node = F,
  seed = 120,
  title = "Louvain 3 Monocytes/macrophages"
)
