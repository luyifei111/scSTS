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

STS.integrated.pega <- readRDS("/cluster3/yflu/STS/pegasus/STS.integrated.pega_20240507.rds")
samplenames <- unique(as.character(STS.integrated.pega$Channel))

cnv_genes <- readRDS("/cluster3/yflu/STS/development/target/cnv_genes.rds")
drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))

cnv_genes_top_list <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_cnv_genes_top_list_251202.rds")
CNV_correlated_drugs_list <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_CNV_correlated_drugs_list_251202.rds")

aucell_druggable_genesets <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets.rds")
aucell_druggable_genesets <- t(aucell_druggable_genesets)
aucell_druggable_genesets_scaled <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets_scaled.rds")

colnames(aucell_druggable_genesets) <- gsub("[|]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub(" ",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("-",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[(]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[)]",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub(",",".",  colnames(aucell_druggable_genesets))
colnames(aucell_druggable_genesets) <- gsub("[/]",".",  colnames(aucell_druggable_genesets))

colnames(aucell_druggable_genesets_scaled) <- gsub("[|]",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub(" ",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub("-",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub("[(]",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub("[)]",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub(",",".",  colnames(aucell_druggable_genesets_scaled))
colnames(aucell_druggable_genesets_scaled) <- gsub("[/]",".",  colnames(aucell_druggable_genesets_scaled))
#degs_tumor_normal_sample_list <- readRDS("/cluster3/yflu/STS/Drug_screen/degs_tumor_normal_sample_list.rds")

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

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

correlated_drugs_counts <- sapply(CNV_correlated_drugs_list, length)
correlated_drugs_counts
correlated_drugs_counts <- cbind(as.data.frame(correlated_drugs_counts),anno_sample_cluster_extended[names(correlated_drugs_counts),])
colnames(correlated_drugs_counts)[1] <- c("CNV_score")
correlated_drugs_counts$Disease <- factor(correlated_drugs_counts$Disease,levels = labels)
plot_cnv_box_by_malignancy(correlated_drugs_counts, plot_title = "CNV_correlated_druggable_genesets_counts")

samplenames <- rownames(correlated_drugs_counts)
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

cnv_score_DSGs <- cbind(cnv_score_sample,correlated_drugs_counts)

colnames(cnv_score_DSGs)[1] <- "CNV_score"
colnames(cnv_score_DSGs)[2] <- "DSG_counts"

full_names <- c("Hemangioma", "KHE", "Schwannoma", "MPNST", "Undifferentiated sarcoma",
                "RMS", "MRT", "IMT", "Angiosarcoma", "EWS/PNET",
                "NF", "Aggressive fibromatosis", "Liposarcoma", "Spindle cell tumor", "ASPS",
                "Infantile fibrosarcoma", "Synovial sarcoma", "Lipoblastoma", "Pecoma", "Lymphangioma")

# 对应缩写
abbreviations <- c("HE","KHE","SWN","MPNST","US","RMS","MRT","IMT","AS","EWS",
                   "NF","AF","LPS","SCT","ASPS","IFS","SS","LPB","PECOMA","LYM")

# 创建替换映射
name_map <- setNames(abbreviations, full_names)
cnv_score_DSGs$Disease <- as.character(cnv_score_DSGs$Disease)
cnv_score_DSGs$Disease <- name_map[cnv_score_DSGs$Disease]

# 疾病顺序
disease_levels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS",
                    "NF","SWN","LYM","HE","KHE","MPNST","US",
                    "AS","RMS","PECOMA","EWS","MRT")

# 对应颜色
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

# 映射到命名向量
names(cols) <- disease_levels

cnv_score_DSGs$Disease <- factor(
  cnv_score_DSGs$Disease,
  levels = disease_levels
)

p <- ggscatter(
  cnv_score_DSGs, 
  x = "CNV_score", 
  y = "DSG_counts", 
  color = "Disease",
  size = 5,
  xlab = "CNV_score", 
  ylab = "DSG_counts"
) +
  geom_smooth(
    method = "lm", 
    aes(x = CNV_score, y = DSG_counts), 
    color = "black", 
    se = TRUE
  ) +
  stat_cor(
    aes(x = CNV_score, y = DSG_counts), 
    method = "spearman"
  ) +
  scale_color_manual(values = cols) +
  theme_classic()
p

library(dplyr)

combine_aucell_drugs_keep_rownames_full <- function(drug_list) {
  sample_names <- names(drug_list)
  df_list <- list()
  
  # 先收集所有可能的 drug 名称
  all_drugs <- character(0)
  for (name in sample_names) {
    df <- drug_list[[name]]
    if (is.null(df) || (is.numeric(df) && length(df) == 1 && df == 0)) next
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    all_drugs <- unique(c(all_drugs, rownames(df)))
  }
  
  # 遍历每个样本
  for (name in sample_names) {
    df <- drug_list[[name]]
    
    # 如果是 NULL 或 0，生成全 0 列
    if (is.null(df) || (is.numeric(df) && length(df) == 1 && df == 0)) {
      df_list[[name]] <- tibble(drug = all_drugs, !!name := 0)
      next
    }
    
    # 确保 df 是 data.frame 并保留行名
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df$drugnames <- rownames(df)
    
    # 检查 aucell_druggable_genesets 列
    if (!"aucell_druggable_genesets" %in% colnames(df)) {
      df_list[[name]] <- tibble(drug = all_drugs, !!name := 0)
      next
    }
    
    # 处理重复 drugnames，取最大值
    df_clean <- df %>%
      group_by(drugnames) %>%
      summarise(val = max(aucell_druggable_genesets, na.rm = TRUE), .groups = "drop") %>%
      mutate(drug = drugnames) %>%
      mutate(!!name := val)
    
    df_clean$val <- NULL
    df_clean$drugnames <- NULL
    
    # 补充缺失 drug 行
    missing_drugs <- setdiff(all_drugs, df_clean$drug)
    if (length(missing_drugs) > 0) {
      df_clean <- bind_rows(df_clean, tibble(drug = missing_drugs, !!name := 0))
    }
    
    # 按 drug 排序
    df_clean <- df_clean[order(df_clean$drug), ]
    
    df_list[[name]] <- df_clean
  }
  
  # 全合并
  merged_df <- Reduce(function(x, y) full_join(x, y, by = "drug"), df_list)
  
  # 确保所有原始样本列都存在
  for (col in sample_names) {
    if (!col %in% colnames(merged_df)) merged_df[[col]] <- 0
  }
  
  # 缺失值补 0
  merged_df <- merged_df %>% replace(is.na(.), 0)
  
  # 排序
  merged_df <- merged_df[order(merged_df$drug), ]
  
  return(merged_df)
}

CNV_aucell_drugs_top_list <- readRDS("/cluster3/yflu/STS/Drug_screen/CNV_aucell_drugs_top_list_251203.rds")

# 使用
CNV_aucell_drugs_merged <- combine_aucell_drugs_keep_rownames_full(CNV_aucell_drugs_top_list)
head(CNV_aucell_drugs_merged)

CNV_aucell_drugs_merged$nonzero_count <- apply(
  CNV_aucell_drugs_merged[ , setdiff(colnames(CNV_aucell_drugs_merged), "drug")],
  1,
  function(x) sum(x != 0 & !is.na(x))
)
CNV_aucell_drugs_merged <- as.data.frame(CNV_aucell_drugs_merged)
rownames(CNV_aucell_drugs_merged) <- CNV_aucell_drugs_merged$drug

CNV_aucell_drugs_merged_sub <- subset(CNV_aucell_drugs_merged,nonzero_count > 10)

orig_names <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets.rds")
orig_names <- rownames(orig_names)
orig_names <- cbind(as.data.frame(orig_names),as.data.frame(orig_names))
colnames(orig_names) <- c("orig","changed")

orig_names$changed <- gsub("[|]",".",  orig_names$changed)
orig_names$changed <- gsub(" ",".",  orig_names$changed)
orig_names$changed <- gsub("-",".",  orig_names$changed)
orig_names$changed <- gsub("[(]",".",  orig_names$changed)
orig_names$changed <- gsub("[)]",".",  orig_names$changed)
orig_names$changed <- gsub(",",".",  orig_names$changed)
orig_names$changed <- gsub("[/]",".",  orig_names$changed)

rownames(orig_names) <- orig_names$changed
rownames(CNV_aucell_drugs_merged_sub) <- orig_names[rownames(CNV_aucell_drugs_merged_sub),]$orig
names_new <- sub("-.*$", "", rownames(CNV_aucell_drugs_merged_sub))

# 首字母大写，其余小写
rownames(CNV_aucell_drugs_merged_sub) <- paste0(toupper(substr(names_new, 1, 1)),
                                                tolower(substr(names_new, 2, nchar(names_new))))
CNV_aucell_drugs_merged_sub_1 <- CNV_aucell_drugs_merged_sub[,-c(1,80)]

pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_1,annotation_col = anno_sample_cluster_extended)

cols_to_check <- setdiff(colnames(CNV_aucell_drugs_merged_sub), c("drug", "nonzero_count"))
# 每行非零列名保存为 list
row_nonzero_cols <- apply(CNV_aucell_drugs_merged_sub[, cols_to_check], 1, function(x) {
  names(x)[x != 0 & !is.na(x)]
})
# 转成 list（apply 返回的是 matrix 时可能自动转成 vector，所以加 as.list 保证 list）
row_nonzero_cols <- as.list(row_nonzero_cols)
names(row_nonzero_cols) <- CNV_aucell_drugs_merged_sub$drug
row_nonzero_cols
druggenes_list <- drug_hgnc_list_combined[CNV_aucell_drugs_merged_sub$drug]

drug_sample_gene_intersect <- list()

for (drug in names(row_nonzero_cols)) {
  
  # 有值的样本
  samples_with_value <- row_nonzero_cols[[drug]]
  
  # drug 的 target genes
  drug_targets <- druggenes_list[[drug]]
  
  # 初始化该药物的结果
  intersect_list <- list()
  
  for (sample in samples_with_value) {
    # 检查 cnv_genes_top_list 中是否有该 sample
    if (!is.null(cnv_genes_top_list[[sample]])) {
      sample_genes <- cnv_genes_top_list[[sample]]$genes
      # 取交集
      intersect_genes <- intersect(drug_targets, sample_genes)
      intersect_list[[sample]] <- intersect_genes
    } else {
      intersect_list[[sample]] <- character(0)
    }
  }
  
  # 保存该药物的样本-基因交集
  drug_sample_gene_intersect[[drug]] <- intersect_list
}

cnv_genes_top_list <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_cnv_genes_top_list_251202.rds")

library(dplyr)
library(stringr)

replace_summary <- list()

drug_sample_gene_intersect_updated <- lapply(names(drug_sample_gene_intersect), function(drug) {
  sample_list <- drug_sample_gene_intersect[[drug]]
  
  lapply(names(sample_list), function(sample_id) {
    genes <- sample_list[[sample_id]]
    
    # 如果该样本存在于 cnv_genes_top_list
    if (sample_id %in% names(cnv_genes_top_list)) {
      cnv_df <- cnv_genes_top_list[[sample_id]]
      
      # 对每个基因依次判断是否需要替换
      updated_genes <- unlist(lapply(genes, function(g) {
        matched <- cnv_df %>% filter(genes == g)
        
        if (nrow(matched) > 0 && all(matched$status == "TF target")) {
          cnv_prefixes <- unique(str_extract(matched$CNV, "^[^_]+"))
          
          # 记录替换信息
          replace_summary[[length(replace_summary) + 1]] <<- data.frame(
            Drug = drug,
            Sample = sample_id,
            Gene_original = g,
            Gene_replaced = paste(cnv_prefixes, collapse = ","),
            stringsAsFactors = FALSE
          )
          
          # ✅ 如果有多个 CNV → 返回字符向量
          return(cnv_prefixes)
          
        } else {
          # 未替换
          replace_summary[[length(replace_summary) + 1]] <<- data.frame(
            Drug = drug,
            Sample = sample_id,
            Gene_original = g,
            Gene_replaced = g,
            stringsAsFactors = FALSE
          )
          return(g)
        }
      }))
      
      # ✅ 输出结构与原始一致（字符向量）
      return(updated_genes)
      
    } else {
      # 样本不在 cnv_genes_top_list → 原样返回
      replace_summary[[length(replace_summary) + 1]] <<- data.frame(
        Drug = drug,
        Sample = sample_id,
        Gene_original = genes,
        Gene_replaced = genes,
        stringsAsFactors = FALSE
      )
      return(genes)
    }
  }) |> rlang::set_names(names(sample_list))
}) |> rlang::set_names(names(drug_sample_gene_intersect))

# 汇总表
replace_summary_df <- bind_rows(replace_summary)

drug_sample_gene_intersect <- drug_sample_gene_intersect_updated

drug_gene_counts <- lapply(drug_sample_gene_intersect, function(sample_list) {
  # 提取所有样本交集基因
  all_genes <- unlist(sample_list, use.names = FALSE)
  
  # 统计频数
  gene_counts <- sort(table(all_genes), decreasing = TRUE)
  
  # 如果只有一个基因，手动构建 data.frame
  if (length(gene_counts) == 1) {
    df <- data.frame(
      all_genes = names(gene_counts),
      Freq = as.numeric(gene_counts),
      stringsAsFactors = FALSE
    )
  } else {
    df <- as.data.frame(gene_counts)
    colnames(df) <- c("all_genes", "Freq")
  }
  
  return(df)
})

library(dplyr)

extract_top_gene_per_drug <- function(drug_gene_counts) {
  out <- lapply(names(drug_gene_counts), function(drug) {
    x <- drug_gene_counts[[drug]]
    # empty
    if (is.null(x) || length(x) == 0) {
      return(data.frame(drug = drug, top_gene = NA_character_, count = 0, stringsAsFactors = FALSE))
    }
    
    # case: table
    if (is.table(x)) {
      df <- data.frame(gene = names(x), count = as.numeric(x), stringsAsFactors = FALSE)
      
      # atomic vector (named or not)
    } else if (is.atomic(x) && !is.list(x)) {
      nm <- names(x)
      if (!is.null(nm)) {
        # named numeric/character vector: names are genes, values are counts (or 1)
        df <- data.frame(gene = nm, count = as.numeric(x), stringsAsFactors = FALSE)
      } else if (is.character(x)) {
        # un-named character vector: entries are gene names, freq = 1 each
        df <- data.frame(gene = as.character(x), count = rep(1, length(x)), stringsAsFactors = FALSE)
      } else {
        # un-named numeric vector (unlikely): can't infer gene names -> skip
        return(data.frame(drug = drug, top_gene = NA_character_, count = 0, stringsAsFactors = FALSE))
      }
      
      # data.frame
    } else if (is.data.frame(x)) {
      cols <- colnames(x)
      
      # case: single column data.frame
      if (ncol(x) == 1) {
        # if rownames look informative (not "1","2",...), use them as genes
        rn <- rownames(x)
        if (!is.null(rn) && any(rn != as.character(seq_along(rn)))) {
          df <- data.frame(gene = rn, count = as.numeric(x[[1]]), stringsAsFactors = FALSE)
        } else {
          # otherwise assume the single column contains gene names
          if (is.character(x[[1]]) || is.factor(x[[1]])) {
            df <- data.frame(gene = as.character(x[[1]]), count = rep(1, nrow(x)), stringsAsFactors = FALSE)
          } else {
            # fallback: use values as counts but no gene names -> skip
            return(data.frame(drug = drug, top_gene = NA_character_, count = 0, stringsAsFactors = FALSE))
          }
        }
        
      } else {
        # multi-column data.frame: try to detect gene column (character/factor) and count column (numeric)
        gene_col_idx <- which(sapply(x, function(col) is.character(col) || is.factor(col)))[1]
        count_col_idx <- which(sapply(x, is.numeric))[1]
        
        # also try heuristic by name
        if (is.na(gene_col_idx) && any(grepl("gene", cols, ignore.case = TRUE))) {
          gene_col_idx <- which(grepl("gene", cols, ignore.case = TRUE))[1]
        }
        if (is.na(count_col_idx) && any(grepl("count|freq|n$", cols, ignore.case = TRUE))) {
          count_col_idx <- which(grepl("count|freq|n$", cols, ignore.case = TRUE))[1]
        }
        
        # fallback: first col = gene, second col = count
        if (is.na(gene_col_idx) || is.na(count_col_idx)) {
          gene_col_idx <- ifelse(is.na(gene_col_idx), 1, gene_col_idx)
          count_col_idx <- ifelse(is.na(count_col_idx), min(setdiff(seq_along(cols), gene_col_idx)), count_col_idx)
        }
        
        df <- data.frame(
          gene = as.character(x[[gene_col_idx]]),
          count = as.numeric(x[[count_col_idx]]),
          stringsAsFactors = FALSE
        )
      }
      
    } else {
      # other types: skip
      return(data.frame(drug = drug, top_gene = NA_character_, count = 0, stringsAsFactors = FALSE))
    }
    
    # clean and choose top
    df <- df[!is.na(df$gene) & df$gene != "", , drop = FALSE]
    if (nrow(df) == 0) return(data.frame(drug = drug, top_gene = NA_character_, count = 0, stringsAsFactors = FALSE))
    
    # if any NA counts, set to 1 (conservative)
    df$count[is.na(df$count)] <- 1
    
    top_idx <- which.max(df$count)
    data.frame(drug = drug, top_gene = as.character(df$gene[top_idx]), count = as.numeric(df$count[top_idx]), stringsAsFactors = FALSE)
  })
  
  res <- bind_rows(out)
  rownames(res) <- NULL
  res
}

# 调用：
top_gene_per_drug <- extract_top_gene_per_drug(drug_gene_counts)
top_gene_per_drug

library(S4Vectors)
library(tidyr)

infer.gistic <- readRDS("/cluster3/yflu/STS/WES_CNV/infer.gistic.rds")

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
g_score <- infer.gistic@gis.scores

library(data.table)

## 确保都是 data.table
setDT(g_score)
setDT(recurrent_cnv)

## 统一染色体格式（非常重要）
g_score[, Chromosome := paste0("chr", Chromosome)]

## 设置区间 key
setkey(g_score, Chromosome, Start_Position, End_Position)
setkey(recurrent_cnv, Chromosome, Start, End)

## 区间重叠
ov <- foverlaps(
  g_score,
  recurrent_cnv,
  by.x = c("Chromosome", "Start_Position", "End_Position"),
  by.y = c("Chromosome", "Start", "End"),
  type = "any",
  nomatch = 0L
)

## 对每个 CNV 取 G_Score 最大值
cnv_max_gscore <- ov[
  ,
  .(max_G_Score = max(G_Score, na.rm = TRUE)),
  by = Unique_Name
]

## 合并回 recurrent_cnv
recurrent_cnv <- merge(
  recurrent_cnv,
  cnv_max_gscore,
  by = "Unique_Name",
  all.x = TRUE
)

write.csv(recurrent_cnv,"recurrent_cnv.csv")

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

gene_list <- List()
gene_data <- infer.gistic@data
for (i in 1:nrow(recurrent_cnv)) {
  gene_list <- append(gene_list,list(unique(subset(gene_data,Cytoband == recurrent_cnv$Unique_Name[i])$Hugo_Symbol)))
  names(gene_list)[i+1] <- recurrent_cnv$Cytoband[i]
}
gene_list[[2]] <- gene_list[[1]]
gene_list <- gene_list[-1]

map_top_gene_to_gene_list <- function(top_gene_per_drug, gene_list) {
  # 确保输入正确
  if (!("top_gene" %in% colnames(top_gene_per_drug))) {
    stop("`top_gene_per_drug` 必须包含列 'top_gene'")
  }
  
  # 遍历每个药物
  result <- lapply(seq_len(nrow(top_gene_per_drug)), function(i) {
    gene <- top_gene_per_drug$top_gene[i]
    drug <- top_gene_per_drug$drug[i]
    count <- top_gene_per_drug$count[i]
    
    # 找出 gene 出现在 gene_list 中的样本
    matched <- names(gene_list)[vapply(gene_list, function(gset) gene %in% gset, logical(1))]
    
    # 如果匹配为空则设为 NA
    if (length(matched) == 0) matched <- NA_character_
    
    data.frame(
      drug = drug,
      top_gene = gene,
      count = count,
      matched_samples = paste(matched, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  # 合并结果
  result_df <- do.call(rbind, result)
  return(result_df)
}
mapped_df <- map_top_gene_to_gene_list(top_gene_per_drug, gene_list)

drug_gene_matched <- lapply(names(drug_gene_counts), function(drug) {
  df <- drug_gene_counts[[drug]]
  
  # 如果 df 为空，直接返回 NA
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(all_genes = NA_character_, Freq = NA_real_, matched_samples = NA_character_, stringsAsFactors = FALSE))
  }
  
  # 遍历每个基因，检查属于 gene_list 的哪些元素
  matched_samples <- sapply(df$all_genes, function(gene) {
    matched <- names(gene_list)[vapply(gene_list, function(gset) gene %in% gset, logical(1))]
    if (length(matched) == 0) NA_character_ else paste(matched, collapse = ";")
  }, USE.NAMES = FALSE)
  
  # 返回 data.frame
  data.frame(
    all_genes = df$all_genes,
    Freq = df$Freq,
    matched_samples = matched_samples,
    stringsAsFactors = FALSE
  )
})

# 保留药物名字
names(drug_gene_matched) <- names(drug_gene_counts)

extract_top_rows <- function(drug_gene_matched) {
  result_list <- lapply(names(drug_gene_matched), function(drug) {
    df <- drug_gene_matched[[drug]]
    
    # 移除空数据框
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # 1️⃣ Freq 最大的一行
    top_all <- df[which.max(df$Freq), , drop = FALSE]
    top_all$drug <- drug
    top_all$type <- "all"
    
    # 2️⃣ matched_samples 不为 NA 且 Freq 最大的一行
    df_non_na <- df[!is.na(df$matched_samples), , drop = FALSE]
    if (nrow(df_non_na) > 0) {
      top_non_na <- df_non_na[which.max(df_non_na$Freq), , drop = FALSE]
      top_non_na$drug <- drug
      top_non_na$type <- "non_na"
    } else {
      top_non_na <- NULL
    }
    
    rbind(top_all, top_non_na)
  })
  
  # 合并所有药物结果
  result_df <- do.call(rbind, result_list)
  rownames(result_df) <- NULL
  return(result_df)
}

# 使用
top_drug_gene_df <- extract_top_rows(drug_gene_matched)

top_matched_per_drug <- lapply(names(drug_gene_matched), function(drug) {
  df <- drug_gene_matched[[drug]]
  
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      drug = drug,
      top_gene = NA_character_,
      Freq = NA_real_,
      matched_samples = NA_character_,
      matched_count = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  # 确保 matched_samples 是字符向量
  df$matched_samples <- as.character(df$matched_samples)
  
  # 计算每行 matched_samples 的数量
  df$matched_count <- sapply(df$matched_samples, function(ms) {
    if (is.na(ms)) 0 else 1
  })
  df$Freq <- df$Freq * df$matched_count
  # 选出 matched_count 最大的行
  top_idx <- which.max(df$Freq)
  top_row <- df[top_idx, ]
  
  # 返回整理好的行
  data.frame(
    drug = drug,
    top_gene = top_row$all_genes,
    Freq = top_row$Freq,
    matched_samples = top_row$matched_samples,
    matched_count = top_row$matched_count,
    stringsAsFactors = FALSE
  )
})

# 合并成一个 data.frame
top_matched_per_drug_df <- do.call(rbind, top_matched_per_drug)
top_matched_per_drug_df$matched_samples_amp <- sapply(top_matched_per_drug_df$matched_samples, function(ms) {
  # 如果 NA 直接返回 NA
  if (is.na(ms)) return(NA_character_)
  
  # 拆分分号
  samples <- unlist(strsplit(ms, ";"))
  
  # 只保留以 "Amp" 结尾的
  amp_samples <- samples[grepl("Amp$", samples)]
  
  # 重新组合成字符串，如果没有匹配返回 NA
  if (length(amp_samples) == 0) return(NA_character_) else paste(amp_samples, collapse = ";")
})
rownames(top_matched_per_drug_df) <- rownames(CNV_aucell_drugs_merged_sub_1)

library(biomaRt)

# 连接 Ensembl（人类基因组）
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_to_fill <- unique(top_matched_per_drug_df$top_gene[is.na(top_matched_per_drug_df$matched_samples_amp)])

if (length(genes_to_fill) > 0) {
  gene_cytoband <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "band"),
    filters = "hgnc_symbol",
    values = genes_to_fill,
    mart = ensembl
  ) %>%
    mutate(cytoband = paste0(chromosome_name, band),
           matched_samples = paste0(cytoband, "_Amp"),
           matched_samples_amp = matched_samples)
}

# Step 3. 补全 NA
top_matched_per_drug_df <- top_matched_per_drug_df %>%
  left_join(gene_cytoband[, c("hgnc_symbol", "matched_samples", "matched_samples_amp")],
            by = c("top_gene" = "hgnc_symbol"),
            suffix = c("", "_fill")) %>%
  mutate(
    matched_samples = ifelse(is.na(matched_samples), matched_samples_fill, matched_samples),
    matched_samples_amp = ifelse(is.na(matched_samples_amp), matched_samples_amp_fill, matched_samples_amp)
  ) %>%
  dplyr::select(-ends_with("_fill"))

# Step 4. 结果查看
rownames(top_matched_per_drug_df) <- rownames(CNV_aucell_drugs_merged_sub_1)

top_matched_per_drug_df <- top_matched_per_drug_df %>%
  mutate(
    matched_samples_amp = if_else(
      is.na(matched_samples_amp) | matched_samples_amp == "",
      matched_samples,
      matched_samples_amp
    )
  )

# 先拆分 chr 和 arm+band
chr_ordered <- sapply(top_matched_per_drug_df$matched_samples_amp, function(x) {
  if (is.na(x)) return(NA_character_)
  # 去掉 _Amp
  x_clean <- sub("_Amp$", "", x)
  # 提取 chr（数字或 X/Y）和位置
  x_clean
})

# 先拆分 chr 和 arm+band
chr_ordered <- sapply(top_matched_per_drug_df$matched_samples_amp, function(x) {
  if (is.na(x)) return(NA_character_)
  # 去掉 _Amp
  x_clean <- sub("_Amp$", "", x)
  # 提取 chr（数字或 X/Y）和位置
  x_clean
})

# 自定义排序函数：先按 chr（数字优先，X=23,Y=24），再按位置
parse_chr_pos <- function(s) {
  # s 示例: "7q21.3"
  m <- regmatches(s, regexec("^([0-9XY]+)(.*)$", s))[[1]]
  chr <- m[2]
  pos <- m[3]
  
  # 数字化 chr
  chr_num <- if (chr %in% c("X","Y")) ifelse(chr=="X", 23, 24) else as.numeric(chr)
  
  # 解析位置部分为可排序数字，去掉 q/p 字母
  pos_num <- as.numeric(gsub("[^0-9.]", "", pos))
  
  c(chr_num, pos_num)
}

# 排序索引
order_idx <- order(sapply(chr_ordered, function(s) parse_chr_pos(s)[1]),
                   sapply(chr_ordered, function(s) parse_chr_pos(s)[2]))

# 排序后的向量

matched_samples_amp_sorted <- unique(top_matched_per_drug_df$matched_samples_amp[order_idx])
top_matched_per_drug_df$matched_samples_amp <- factor(top_matched_per_drug_df$matched_samples_amp,levels = matched_samples_amp_sorted)
top_matched_per_drug_df <- top_matched_per_drug_df[order(top_matched_per_drug_df$matched_samples_amp),]

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

annotation_colors <- c(
  list(Disease = disease_colors),
  list(Malignancy = malignancy_colors)
)

row_colors <- c(
  "1q21.3_Amp"   = "#FC8D62", 
  "2p25.1_Amp"   = "#66C2A5",
  "2q37.1_Amp"  = "#A6D854",
  "5q31.2_Amp" = "#ADD8E6",
  "6p22.1_Amp"   = "#D9DFED",  
  "7q21.3_Amp" = "#5E81AC",
  "11q13.2_Amp"  = "#B07CCF",
  "14q32.31_Amp" = "#1F78B4",
  "17q25.3_Amp"  = "#E377C2",
  "19p13.2_Amp" = "#FFB347",
  "20p13_Amp" = "#1F7874",
  "21q22.13_Amp" = "#8DA0CB",
  "22q11.21_Amp" = "#FFD92F"
)

row_annotation_colors <- list(matched_samples_amp = row_colors)
anno_row <- as.data.frame(top_matched_per_drug_df[,-c(1,2,3,4,5)])
rownames(anno_row) <- rownames(top_matched_per_drug_df)
colnames(anno_row) <- "matched_samples_amp"

colors_combined = colorRampPalette(brewer.pal(8,'RdBu'))(100)

pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_1[rownames(top_matched_per_drug_df),],
                   cluster_rows = F,
                   color = rev(c(colors_combined[1:40],"#FFFFFF")),
                   annotation_col = anno_sample_cluster_extended[,-3],
                   annotation_row = anno_row,
                   clustering_distance_cols = "manhattan",
                   annotation_colors = c(annotation_colors,row_annotation_colors))


library(ggplot2)
library(dplyr)

# 假设你的数据叫 CNV_disease_durgs_sub
CNV_disease_durgs_sub <- subset(CNV_disease_durgs_sub,logFC != 0)
df <- CNV_disease_durgs_sub

# 为了美观，可以让疾病顺序按 prop 或 logFC 的平均值排列
df <- df %>%
  mutate(
    Disease = factor(Disease, levels = rev(unique(Disease))),
    Drugs = factor(Drugs, levels = unique(Drugs[order(logFC, decreasing = TRUE)]))
  )

df$Drugs <- orig_names[df$Drugs,]$orig
names_new <- sub("-.*$", "", df$Drugs)

# 首字母大写，其余小写
df$Drugs <- paste0(toupper(substr(names_new, 1, 1)),
                                                tolower(substr(names_new, 2, nchar(names_new))))
disease_drugs <- as.character(unique(CNV_disease_durgs_sub$Drugs))
CNV_aucell_drugs_merged_sub <- subset(CNV_aucell_drugs_merged,drug %in% disease_drugs)

orig_names <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets.rds")
orig_names <- rownames(orig_names)
orig_names <- cbind(as.data.frame(orig_names),as.data.frame(orig_names))
colnames(orig_names) <- c("orig","changed")

orig_names$changed <- gsub("[|]",".",  orig_names$changed)
orig_names$changed <- gsub(" ",".",  orig_names$changed)
orig_names$changed <- gsub("-",".",  orig_names$changed)
orig_names$changed <- gsub("[(]",".",  orig_names$changed)
orig_names$changed <- gsub("[)]",".",  orig_names$changed)
orig_names$changed <- gsub(",",".",  orig_names$changed)
orig_names$changed <- gsub("[/]",".",  orig_names$changed)

rownames(orig_names) <- orig_names$changed
rownames(CNV_aucell_drugs_merged_sub) <- orig_names[rownames(CNV_aucell_drugs_merged_sub),]$orig
names_new <- sub("-.*$", "", rownames(CNV_aucell_drugs_merged_sub))

# 首字母大写，其余小写
rownames(CNV_aucell_drugs_merged_sub) <- paste0(toupper(substr(names_new, 1, 1)),
                                                tolower(substr(names_new, 2, nchar(names_new))))
CNV_aucell_drugs_merged_sub_1 <- CNV_aucell_drugs_merged_sub[,-c(1,80)]

pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_1,annotation_col = anno_sample_cluster_extended)

cols_to_check <- setdiff(colnames(CNV_aucell_drugs_merged_sub), c("drug", "nonzero_count"))
# 每行非零列名保存为 list
row_nonzero_cols <- apply(CNV_aucell_drugs_merged_sub[, cols_to_check], 1, function(x) {
  names(x)[x != 0 & !is.na(x)]
})
# 转成 list（apply 返回的是 matrix 时可能自动转成 vector，所以加 as.list 保证 list）
row_nonzero_cols <- as.list(row_nonzero_cols)
names(row_nonzero_cols) <- CNV_aucell_drugs_merged_sub$drug
row_nonzero_cols
druggenes_list <- drug_hgnc_list_combined[CNV_aucell_drugs_merged_sub$drug]

drug_sample_gene_intersect <- list()

for (drug in names(row_nonzero_cols)) {
  
  # 有值的样本
  samples_with_value <- row_nonzero_cols[[drug]]
  
  # drug 的 target genes
  drug_targets <- druggenes_list[[drug]]
  
  # 初始化该药物的结果
  intersect_list <- list()
  
  for (sample in samples_with_value) {
    # 检查 cnv_genes_top_list 中是否有该 sample
    if (!is.null(cnv_genes_top_list[[sample]])) {
      sample_genes <- cnv_genes_top_list[[sample]]$genes
      # 取交集
      intersect_genes <- intersect(drug_targets, sample_genes)
      intersect_list[[sample]] <- intersect_genes
    } else {
      intersect_list[[sample]] <- character(0)
    }
  }
  
  # 保存该药物的样本-基因交集
  drug_sample_gene_intersect[[drug]] <- intersect_list
}

library(dplyr)
library(stringr)

replace_summary <- list()

drug_sample_gene_intersect_updated <- lapply(names(drug_sample_gene_intersect), function(drug) {
  sample_list <- drug_sample_gene_intersect[[drug]]
  
  lapply(names(sample_list), function(sample_id) {
    genes <- sample_list[[sample_id]]
    
    # 如果该样本存在于 cnv_genes_top_list
    if (sample_id %in% names(cnv_genes_top_list)) {
      cnv_df <- cnv_genes_top_list[[sample_id]]
      
      # 对每个基因依次判断是否需要替换
      updated_genes <- unlist(lapply(genes, function(g) {
        matched <- cnv_df %>% filter(genes == g)
        
        if (nrow(matched) > 0 && all(matched$status == "TF target")) {
          cnv_prefixes <- unique(str_extract(matched$CNV, "^[^_]+"))
          
          # 记录替换信息
          replace_summary[[length(replace_summary) + 1]] <<- data.frame(
            Drug = drug,
            Sample = sample_id,
            Gene_original = g,
            Gene_replaced = paste(cnv_prefixes, collapse = ","),
            stringsAsFactors = FALSE
          )
          
          # ✅ 如果有多个 CNV → 返回字符向量
          return(cnv_prefixes)
          
        } else {
          # 未替换
          replace_summary[[length(replace_summary) + 1]] <<- data.frame(
            Drug = drug,
            Sample = sample_id,
            Gene_original = g,
            Gene_replaced = g,
            stringsAsFactors = FALSE
          )
          return(g)
        }
      }))
      
      # ✅ 输出结构与原始一致（字符向量）
      return(updated_genes)
      
    } else {
      # 样本不在 cnv_genes_top_list → 原样返回
      replace_summary[[length(replace_summary) + 1]] <<- data.frame(
        Drug = drug,
        Sample = sample_id,
        Gene_original = genes,
        Gene_replaced = genes,
        stringsAsFactors = FALSE
      )
      return(genes)
    }
  }) |> rlang::set_names(names(sample_list))
}) |> rlang::set_names(names(drug_sample_gene_intersect))

# 汇总表
replace_summary_df <- bind_rows(replace_summary)

drug_sample_gene_intersect <- drug_sample_gene_intersect_updated

drug_gene_counts <- lapply(drug_sample_gene_intersect, function(sample_list) {
  # 提取所有样本交集基因
  all_genes <- unlist(sample_list, use.names = FALSE)
  
  # 统计频数
  gene_counts <- sort(table(all_genes), decreasing = TRUE)
  
  # 如果只有一个基因，手动构建 data.frame
  if (length(gene_counts) == 1) {
    df <- data.frame(
      all_genes = names(gene_counts),
      Freq = as.numeric(gene_counts),
      stringsAsFactors = FALSE
    )
  } else {
    df <- as.data.frame(gene_counts)
    colnames(df) <- c("all_genes", "Freq")
  }
  
  return(df)
})

top_gene_per_drug <- extract_top_gene_per_drug(drug_gene_counts)

library(S4Vectors)
library(tidyr)

infer.gistic <- readRDS("/cluster3/yflu/STS/WES_CNV/infer.gistic.rds")

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

mapped_df <- map_top_gene_to_gene_list(top_gene_per_drug, gene_list)

drug_gene_matched <- lapply(names(drug_gene_counts), function(drug) {
  df <- drug_gene_counts[[drug]]
  
  # 如果 df 为空，直接返回 NA
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(all_genes = NA_character_, Freq = NA_real_, matched_samples = NA_character_, stringsAsFactors = FALSE))
  }
  
  # 遍历每个基因，检查属于 gene_list 的哪些元素
  matched_samples <- sapply(df$all_genes, function(gene) {
    matched <- names(gene_list)[vapply(gene_list, function(gset) gene %in% gset, logical(1))]
    if (length(matched) == 0) NA_character_ else paste(matched, collapse = ";")
  }, USE.NAMES = FALSE)
  
  # 返回 data.frame
  data.frame(
    all_genes = df$all_genes,
    Freq = df$Freq,
    matched_samples = matched_samples,
    stringsAsFactors = FALSE
  )
})
names(drug_gene_matched) <- names(drug_gene_counts)

top_matched_per_drug <- lapply(names(drug_gene_matched), function(drug) {
  df <- drug_gene_matched[[drug]]
  
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      drug = drug,
      top_gene = NA_character_,
      Freq = NA_real_,
      matched_samples = NA_character_,
      matched_count = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  # 确保 matched_samples 是字符向量
  df$matched_samples <- as.character(df$matched_samples)
  
  # 计算每行 matched_samples 的数量
  df$matched_count <- sapply(df$matched_samples, function(ms) {
    if (is.na(ms)) 0 else 1
  })
  df$Freq <- df$Freq * df$matched_count
  # 选出 matched_count 最大的行
  top_idx <- which.max(df$Freq)
  top_row <- df[top_idx, ]
  
  # 返回整理好的行
  data.frame(
    drug = drug,
    top_gene = top_row$all_genes,
    Freq = top_row$Freq,
    matched_samples = top_row$matched_samples,
    matched_count = top_row$matched_count,
    stringsAsFactors = FALSE
  )
})

top_matched_per_drug_df <- do.call(rbind, top_matched_per_drug)
top_matched_per_drug_df$matched_samples_amp <- sapply(top_matched_per_drug_df$matched_samples, function(ms) {
  # 如果 NA 直接返回 NA
  if (is.na(ms)) return(NA_character_)
  
  # 拆分分号
  samples <- unlist(strsplit(ms, ";"))
  
  # 只保留以 "Amp" 结尾的
  amp_samples <- samples[grepl("Amp$", samples)]
  
  # 重新组合成字符串，如果没有匹配返回 NA
  if (length(amp_samples) == 0) return(NA_character_) else paste(amp_samples, collapse = ";")
})
rownames(top_matched_per_drug_df) <- rownames(CNV_aucell_drugs_merged_sub_1)

library(biomaRt)

# 连接 Ensembl（人类基因组）
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_to_fill <- unique(top_matched_per_drug_df$top_gene[is.na(top_matched_per_drug_df$matched_samples_amp)])

if (length(genes_to_fill) > 0) {
  gene_cytoband <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "band"),
    filters = "hgnc_symbol",
    values = genes_to_fill,
    mart = ensembl
  ) %>%
    mutate(cytoband = paste0(chromosome_name, band),
           matched_samples = paste0(cytoband, "_Amp"),
           matched_samples_amp = matched_samples)
}

# Step 3. 补全 NA
top_matched_per_drug_df <- top_matched_per_drug_df %>%
  left_join(gene_cytoband[, c("hgnc_symbol", "matched_samples", "matched_samples_amp")],
            by = c("top_gene" = "hgnc_symbol"),
            suffix = c("", "_fill")) %>%
  mutate(
    matched_samples = ifelse(is.na(matched_samples), matched_samples_fill, matched_samples),
    matched_samples_amp = ifelse(is.na(matched_samples_amp), matched_samples_amp_fill, matched_samples_amp)
  ) %>%
  dplyr::select(-ends_with("_fill"))

top_matched_per_drug_df <- top_matched_per_drug_df[-9,]
# Step 4. 结果查看
rownames(top_matched_per_drug_df) <- rownames(CNV_aucell_drugs_merged_sub_1)

# 先拆分 chr 和 arm+band
chr_ordered <- sapply(top_matched_per_drug_df$matched_samples_amp, function(x) {
  if (is.na(x)) return(NA_character_)
  # 去掉 _Amp
  x_clean <- sub("_Amp$", "", x)
  # 提取 chr（数字或 X/Y）和位置
  x_clean
})

# 先拆分 chr 和 arm+band
chr_ordered <- sapply(top_matched_per_drug_df$matched_samples_amp, function(x) {
  if (is.na(x)) return(NA_character_)
  # 去掉 _Amp
  x_clean <- sub("_Amp$", "", x)
  # 提取 chr（数字或 X/Y）和位置
  x_clean
})

# 自定义排序函数：先按 chr（数字优先，X=23,Y=24），再按位置
parse_chr_pos <- function(s) {
  # s 示例: "7q21.3"
  m <- regmatches(s, regexec("^([0-9XY]+)(.*)$", s))[[1]]
  chr <- m[2]
  pos <- m[3]
  
  # 数字化 chr
  chr_num <- if (chr %in% c("X","Y")) ifelse(chr=="X", 23, 24) else as.numeric(chr)
  
  # 解析位置部分为可排序数字，去掉 q/p 字母
  pos_num <- as.numeric(gsub("[^0-9.]", "", pos))
  
  c(chr_num, pos_num)
}

# 排序索引
order_idx <- order(sapply(chr_ordered, function(s) parse_chr_pos(s)[1]),
                   sapply(chr_ordered, function(s) parse_chr_pos(s)[2]))

matched_samples_amp_sorted <- unique(top_matched_per_drug_df$matched_samples_amp[order_idx])
top_matched_per_drug_df$matched_samples_amp <- factor(top_matched_per_drug_df$matched_samples_amp,levels = matched_samples_amp_sorted)
top_matched_per_drug_df <- top_matched_per_drug_df[order(top_matched_per_drug_df$matched_samples_amp),]

disease_colors <- c(
  colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(6,7,5)],
  colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[2],
  "#ADD8E6",
  colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[2:3]
)

# Disease 名称需与你的 sample_anno$Disease levels 顺序一致
names(disease_colors) <- c("AF","IMT","IFS","HE","RMS","EWS","MRT")

annotation_colors <- c(
  list(Disease = disease_colors)
)

row_annotation_colors <- list(matched_samples_amp = row_colors)
anno_row <- as.data.frame(top_matched_per_drug_df[,-c(1,2,3,4,5)])
rownames(anno_row) <- rownames(top_matched_per_drug_df)
colnames(anno_row) <- "matched_samples_amp"

#CNV_disease_durgs_sub <- subset(CNV_disease_durgs_sub,logFC != 0)
CNV_disease_durgs_sub <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_CNV_disease_durgs_sub_251203.rds")
df <- CNV_disease_durgs_sub

# 为了美观，可以让疾病顺序按 prop 或 logFC 的平均值排列
df <- df %>%
  mutate(
    Disease = factor(Disease, levels = rev(unique(Disease))),
    Drugs = factor(Drugs, levels = unique(Drugs[order(logFC, decreasing = TRUE)]))
  )

orig_names <- readRDS("/cluster3/yflu/STS/Drug_screen/aucell_druggable_genesets.rds")
orig_names <- rownames(orig_names)
orig_names <- cbind(as.data.frame(orig_names),as.data.frame(orig_names))
colnames(orig_names) <- c("orig","changed")

orig_names$changed <- gsub("[|]",".",  orig_names$changed)
orig_names$changed <- gsub(" ",".",  orig_names$changed)
orig_names$changed <- gsub("-",".",  orig_names$changed)
orig_names$changed <- gsub("[(]",".",  orig_names$changed)
orig_names$changed <- gsub("[)]",".",  orig_names$changed)
orig_names$changed <- gsub(",",".",  orig_names$changed)
orig_names$changed <- gsub("[/]",".",  orig_names$changed)

rownames(orig_names) <- orig_names$changed

df$Drugs <- orig_names[as.character(df$Drugs),]$orig
names_new <- sub("-.*$", "", df$Drugs)

# 首字母大写，其余小写
df$Drugs <- paste0(toupper(substr(names_new, 1, 1)),
                   tolower(substr(names_new, 2, nchar(names_new))))

df$Drugs <- factor(df$Drugs,levels = rownames(top_matched_per_drug_df))
df$Disease <- factor(df$Disease, levels = names(disease_colors))

library(ggplot2)
library(patchwork)

plot_drug_disease_bubble <- function(df,
                                     anno_row,
                                     row_colors,
                                     disease_colors,
                                     size_title = "Scaled_avg_exp",
                                     color_title = "logFC") {
  #-----------------------------
  # 数据准备：保持原始顺序
  #-----------------------------
  #-----------------------------
  # 主图（气泡图）
  #-----------------------------
  p_main <- ggplot(df, aes(x = Drugs, y = Disease)) +
    geom_point(aes(size = av_exp, color = logFC), alpha = 0.9) +
    scale_size_continuous(range = c(3, 10), name = size_title) +
    scale_color_gradient2(
      low = "blue", mid = "gray95", high = "red",
      midpoint = 0, name = color_title
    ) +
    # 倒置 y 轴顺序
    scale_y_discrete(limits = rev(levels(df$Disease))) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(x = "Drugs", y = "Disease") +
    coord_cartesian(clip = "off")
  
  #-----------------------------
  # 行注释（疾病）
  #-----------------------------
  disease_levels <- levels(df$Disease)
  p_rowanno <- ggplot(data.frame(Disease = disease_levels), aes(x = 1, y = Disease, fill = Disease)) +
    geom_tile() +
    scale_fill_manual(values = disease_colors) +
    scale_y_discrete(limits = rev(disease_levels)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(r = 0)
    )
  #-----------------------------
  # 列注释（药物）
  #-----------------------------
  col_anno_df <- data.frame(
    Drugs = rownames(anno_row),
    matched_samples_amp = anno_row$matched_samples_amp
  )
  
  # 按 df$Drugs 顺序排列
  col_anno_df$Drugs <- factor(col_anno_df$Drugs, levels = levels(df$Drugs))
  
  p_colanno <- ggplot(col_anno_df, aes(x = Drugs, y = 1, fill = matched_samples_amp)) +
    geom_tile() +
    scale_fill_manual(values = row_colors, name = "matched_samples_amp") +
    theme_void() +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(b = 0)
    )
  
  #-----------------------------
  # 拼接图形
  #-----------------------------
  final_plot <- p_colanno /
    (p_rowanno + p_main) +
    plot_layout(widths = c(0.06, 1), heights = c(0.12, 1))
  
  return(final_plot)
}

row_colors <- c(
  "1q21.3_Amp"   = "#FC8D62",
  "1p31.1_Amp" = "#FAAA62",
  "2p25.1_Amp"   = "#66C2A5",
  "2q32.1_Amp" = "#A6BA54",
  "2q37.1_Amp"  = "#A6D854",
  "2q37.3_Amp" = "#A6D8AA",
  "5q31.2_Amp" = "#ADD8E6",
  "5q35.1_Amp" = "#FFC0CB",
  "6p22.1_Amp"   = "#D9DFED",  
  "7q21.3_Amp" = "#5E81AC",
  "11q13.2_Amp"  = "#B07CCF",
  "11p15.2_Amp" = "#B0ACCC",
  "14q32.31_Amp" = "#1F78B4",
  "17q21.33_Amp" = "#EAAAC2",
  "17q25.3_Amp"  = "#E377C2",
  "19p13.2_Amp" = "#FFB347",
  "20p13_Amp" = "#1F7874",
  "21q22.13_Amp" = "#8DA0CB",
  "22q11.21_Amp" = "#FFD92F"
)

final_plot <- plot_drug_disease_bubble(
  df = df,
  anno_row = anno_row,
  row_colors = row_colors,
  disease_colors = disease_colors
)
final_plot

col_anno_df <- data.frame(
  Drugs = rownames(anno_row),
  matched_samples_amp = anno_row$matched_samples_amp
)
col_anno_df

top_drug_gene_df <- extract_top_rows(drug_gene_matched)

top_drug_gene_df$drugs_new <- orig_names[top_drug_gene_df$drug,]$orig
names_new <- sub("-.*$", "", top_drug_gene_df$drugs_new)

top_drug_gene_df$drugs_new <- paste0(toupper(substr(names_new, 1, 1)),
                                     tolower(substr(names_new, 2, nchar(names_new))))

top_drug_gene_df$drugs_new <- factor(top_drug_gene_df$drugs_new,levels = levels(df$Drugs))
top_drug_gene_df <- top_drug_gene_df[order(top_drug_gene_df$drugs_new),]

idx <- match(top_drug_gene_df$drugs_new, rownames(anno_row))

# 替换 matched_samples 列
top_drug_gene_df$matched_samples <- anno_row$matched_samples_amp[idx]

top_drug_genes <- unique(top_drug_gene_df$all_genes)

top_drug_genes_exp <- AverageExpression(STS.integrated.pega,features = top_drug_genes,group.by = "Disease")
top_drug_genes_exp <- as.data.frame(top_drug_genes_exp$RNA)

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

top_drug_gene_group <- top_drug_gene_df[!duplicated(top_drug_gene_df$all_genes), ]

top_drug_gene_anno <- top_drug_gene_df[,c(1,3)]
top_drug_gene_anno <- top_drug_gene_anno[!duplicated(top_drug_gene_anno$all_genes), ]
rownames(top_drug_gene_anno) <- top_drug_gene_anno$all_genes
top_drug_gene_anno_1 <- as.data.frame(top_drug_gene_anno[,-1])
rownames(top_drug_gene_anno_1) <- rownames(top_drug_gene_anno)
colnames(top_drug_gene_anno_1) <- colnames(top_drug_gene_anno)[2]

row_colors <- c(
  "1q21.3_Amp"   = "#FC8D62",
  "1p31.1_Amp" = "#FAAA62",
  "2p25.1_Amp"   = "#66C2A5",
  "2q32.1_Amp" = "#A6BA54",
  "2q37.1_Amp"  = "#A6D854",
  "2q37.3_Amp" = "#A6D8AA",
  "5q31.2_Amp" = "#ADD8E6",
  "5q35.1_Amp" = "#FFC0CB",
  "6p22.1_Amp"   = "#D9DFED",  
  "7q21.3_Amp" = "#5E81AC",
  "11q13.2_Amp"  = "#B07CCF",
  "11p15.2_Amp" = "#B0ACCC",
  "14q32.31_Amp" = "#1F78B4",
  "17q21.33_Amp" = "#EAAAC2",
  "17q25.3_Amp"  = "#E377C2",
  "19p13.2_Amp" = "#FFB347",
  "20p13_Amp" = "#1F7874",
  "21q22.13_Amp" = "#8DA0CB",
  "22q11.21_Amp" = "#FFD92F"
)

annotation_colors <- c(
  list(matched_samples = row_colors)
)

group_var <- "matched_samples"

# 按组排序
ord <- order(top_drug_gene_anno_1[[group_var]])
top_drug_genes_exp2 <- top_drug_genes_exp[ord, ]
annotation_row2 <- top_drug_gene_anno_1[ord, , drop = FALSE]

# 生成 gaps_row：每组最后一行的位置
gaps_row <- cumsum(table(annotation_row2[[group_var]]))

# 绘图
pheatmap::pheatmap(
  top_drug_genes_exp2,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = rev(c(colors_combined[1:40], reds, blues, colors_combined[60:100])),
  annotation_row = annotation_row2,
  annotation_colors = annotation_colors,
  clustering_distance_cols = "manhattan",
  gaps_row = gaps_row   # ← 加分隔空隙的关键
)

#TMA violin
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
library(ggbeeswarm)
library(reshape2)

plot_custom_feature_scaled <- function(
    seurat_obj,
    features,
    data_type = c("assay", "meta"),
    assay_name = "RNA",
    slot = "data",
    reduction = "spatial",
    sample_frac = 1,
    color_palette = viridis::inferno(100),
    point_size = 0.5
) {
  
  data_type <- match.arg(data_type)
  plots <- list()
  
  # 获取坐标
  coords <- Embeddings(seurat_obj, reduction)
  coords <- as.data.frame(coords)
  coords$cell <- rownames(coords)
  
  for(feat in features){
    # 获取表达值
    if(data_type == "assay"){
      expr <- GetAssayData(seurat_obj, assay = assay_name, slot = slot)
      if(!feat %in% rownames(expr)) next
      vals <- expr[feat, ]
    } else {
      if(!feat %in% colnames(seurat_obj@meta.data)) next
      vals <- seurat_obj@meta.data[[feat]]
    }
    
    df <- data.frame(
      cell = names(vals),
      value = as.numeric(vals)
    )
    
    # 合并坐标
    df <- df %>% left_join(coords, by = "cell")
    
    # 采样
    if(sample_frac < 1){
      set.seed(123)
      df <- df[sample(nrow(df), ceiling(nrow(df) * sample_frac)), ]
    }
    
    # 线性缩放到 [-2, 2]
    df$scaled_value <- scales::rescale(df$value, to = c(-2, 2))
    
    # ggplot 绘图
    p <- ggplot(df, aes_string(x = colnames(coords)[1], y = colnames(coords)[2], color = "scaled_value")) +
      geom_point(size = point_size) +
      scale_color_gradientn(colors = color_palette, limits = c(-2, 2)) +
      ggtitle(feat) +
      theme_void() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
      )
    
    plots[[feat]] <- p
  }
  
  # 返回 patchwork 合并图或列表
  if(length(plots) == 1) return(plots[[1]])
  return(wrap_plots(plots))
}

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
#plotHeatmap(aurocs_1)

p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')

TMA_merged_drug_score <- readRDS("/cluster3/yflu/STS/TMA/TMA_merged_drug_score.rds")

#drug_mtx <- read.csv("/cluster3/yflu/STS/TMA/drug_score_mtx.csv")
TMA_merged_drug_score@meta.data <- cbind(TMA_merged_drug_score@meta.data,t(TMA_merged_drug_score@assays$AUCell@data))

#TMA_merged_sub <- subset(TMA_merged_drug_score,Sample == "T1620")
TMA_merged_sub <- subset(TMA_merged_drug_score,Disease %in% c("RMS"))
TMA_merged_sub <- subset(TMA_merged_sub,celltype == "Tumor cells")

TMA_merged_tumor <- subset(TMA_merged_drug_score,celltype == "Tumor cells")

plot_density_with_peak_general <- function(
    seurat_obj,
    features,
    color_vec,
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
  
  # 1. 获取数据
  if(data_type == "meta") {
    # metadata
    meta_df <- seurat_obj@meta.data[, features, drop = FALSE]
    meta_df$cell_id <- rownames(meta_df)
    data_long <- melt(meta_df, id.vars = "cell_id", variable.name = "feature", value.name = "value")
  } else if(data_type == "assay") {
    # assay 数据
    assay_data <- GetAssayData(seurat_obj, assay = assay_name, slot = "data")[features, , drop = FALSE]
    data_long <- as.data.frame(t(as.matrix(assay_data)))
    data_long$cell_id <- rownames(data_long)
    data_long <- melt(data_long, id.vars = "cell_id", variable.name = "feature", value.name = "value")
  }
  
  # 2. 随机采样
  set.seed(seed)
  data_long_sample <- data_long %>%
    group_by(feature) %>%
    sample_frac(sample_frac, replace = FALSE) %>%
    filter(!is.na(value)) %>%
    ungroup()
  
  # 3. factor 顺序按 features 输入顺序，保留实际存在
  feature_levels <- features[features %in% unique(data_long_sample$feature)]
  data_long_sample$feature <- factor(data_long_sample$feature, levels = feature_levels)
  
  # 对 color_vec 只保留实际存在的 feature
  color_vec_plot <- color_vec[names(color_vec) %in% feature_levels]
  
  # 4. 临时 ggplot 获取 quasirandom 点位置
  temp_plot <- ggplot(data_long_sample, aes(x = value, y = feature)) +
    geom_quasirandom(groupOnX = FALSE, varwidth = TRUE, orientation = "y")
  gb <- ggplot_build(temp_plot)
  points_df <- gb$data[[1]]
  
  # 5. 计算峰值和 ycenter
  peak_df <- points_df %>%
    group_by(group) %>%
    summarise(
      npoints = n(),
      peak = if(n() >= 2) {
        density(x)$x[which.max(density(x)$y)]
      } else {
        median(x)
      },
      ycenter = median(y),
      .groups = "drop"
    )
  
  # 6. 绘图
  p <- ggplot(data_long_sample, aes(x = value, y = feature, color = feature)) +
    geom_quasirandom(groupOnX = FALSE, varwidth = TRUE, alpha = 0.5, size = 0.5, orientation = "y") +
    geom_segment(data = peak_df,
                 aes(x = peak, xend = peak, y = ycenter - line_height, yend = ycenter + line_height),
                 color = "black", size = 1) +
    scale_color_manual(values = color_vec_plot) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          legend.position = "none")
  
  return(p)
}

library(Seurat)
library(fields)

compute_weighted_density_python_style <- function(
    seurat_obj,
    value_cols = c("Binimetinib"),
    group_col = "celltype",
    tumor_label = "Tumor cells",
    sigma = NULL,       # 高斯平滑参数
    sample_frac = 1     # 下采样比例
) {
  
  coords <- Embeddings(seurat_obj[["spatial"]])
  meta <- seurat_obj@meta.data
  
  mask_tumor <- meta[[group_col]] == tumor_label
  coords_tumor <- coords[mask_tumor, , drop = FALSE]
  
  # 自动 sigma
  if (is.null(sigma)) {
    dists <- dist(coords_tumor)
    sigma <- median(as.numeric(dists)) / 5
  }
  
  gx <- seq(min(coords_tumor[,1]), max(coords_tumor[,1]), length.out = 200)
  gy <- seq(min(coords_tumor[,2]), max(coords_tumor[,2]), length.out = 200)
  
  for (value_col in value_cols) {
    if (!(value_col %in% colnames(meta))) {
      warning(paste("Skipping", value_col, "- not in metadata"))
      next
    }
    
    values_tumor <- meta[[value_col]][mask_tumor]
    values_tumor[is.na(values_tumor)] <- 0
    
    # 下采样
    if (sample_frac < 1 && sample_frac > 0) {
      set.seed(123)
      sel <- sample(seq_len(nrow(coords_tumor)), size = round(nrow(coords_tumor) * sample_frac))
      coords_tumor_sub <- coords_tumor[sel, , drop = FALSE]
      values_tumor_sub <- values_tumor[sel]
    } else {
      coords_tumor_sub <- coords_tumor
      values_tumor_sub <- values_tumor
    }
    
    # 加权 KDE
    weighted_kde <- function(x, y, z, sigma, gridx, gridy) {
      nx <- length(gridx)
      ny <- length(gridy)
      zmat <- matrix(0, nx, ny)
      for (i in seq_len(nx)) {
        for (j in seq_len(ny)) {
          dx <- x - gridx[i]
          dy <- y - gridy[j]
          w <- exp(-(dx^2 + dy^2) / (2*sigma^2))
          zmat[i, j] <- sum(w * z, na.rm = TRUE)
        }
      }
      # ✅ Python 风格归一化：总和 / 样本数
      zmat <- zmat / sum(zmat, na.rm = TRUE) * length(z)
      list(x = gridx, y = gridy, z = zmat)
    }
    
    w_kde <- weighted_kde(coords_tumor_sub[,1], coords_tumor_sub[,2], values_tumor_sub, sigma, gx, gy)
    dens_interp <- fields::interp.surface(w_kde, coords)
    
    new_col <- paste0("density_", value_col)
    seurat_obj[[new_col]] <- dens_interp
    message("✅ Added Python-style weighted density: metadata$", new_col)
  }
  
  return(seurat_obj)
}
TMA_merged_sub <- compute_weighted_density_python_style(
  seurat_obj = TMA_merged_sub,
  value_cols = c("Epalrestat","Niraparib","Pemigatinib","Cetuximab"),
  group_col = "celltype",
  tumor_label = "Tumor cells",
  sigma = 50
)

cols <- c(
  "density_Binimetinib"   = "#FC8D62", 
  "density_Niraparib"   = "#FC8D62",
  "density_Abemaciclib" = "#5E81AC",
  "density_Crizotinib" = "#8DA0CB"
)

cols_1 <- c(
  "density_Epalrestat"   = "#5E81AC", 
  "density_Niraparib"   = "#FC8D62",
  "density_Pemigatinib" = "#D9DFED",
  "density_Cetuximab" = "#B07CCF"
)

plot_density_with_peak_general(
  seurat_obj = TMA_merged_sub,
  features = rev(names(cols_1)),
  color_vec = cols_1,
  data_type = "meta",
  sample_frac = 0.5
)

FeaturePlot(TMA_merged_sub, features = "density_Cetuximab",reduction = "spatial")


#names(disease_colors) <- disease_order
plot_density_by_group(seurat_obj = TMA_merged_sub,
                      features = meta_cols,
                      color_vec = disease_colors,
                      data_type = "meta",
                      sample_frac = 0.2,group_col = "Disease"
)

rownames(TMA_merged_sub@assays$AUCell) <- gsub("[ /-]", ".", rownames(TMA_merged_sub@assays$AUCell))

cols_scores = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:2)],
                colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:2)],
                colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(4)[c(1:2)],
                colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(2)[c(1)],
                colorRampPalette(brewer.pal(12,'Set3')[c(11,12)])(8)[c(1:3)],
                colorRampPalette(brewer.pal(12,'Set3')[c(9,10)])(8)[c(1:3)])
cellnames <- substr(colnames(aurocs_disease),13,nchar(colnames(aurocs_disease)))
cellnames <- as.data.frame(cbind(cellnames,cellnames))
colnames(cellnames) <- c("Var1","Var2")
rownames(cellnames) <- paste("Development|",cellnames$Var1,sep="")
cellnames <- cellnames[colnames(aurocs_disease)[p$tree_col$order],]
names(cols_scores) <- cellnames$Var1

TMA_merged_sub <- subset(TMA_merged_drug_score,Disease %in% c("AS"))
TMA_merged_sub <- subset(TMA_merged_sub,celltype == "Tumor cells")

plot_density_with_peak_general(
  seurat_obj = TMA_merged_sub,
  features = meta_cols,
  color_vec = cols,
  data_type = "meta",
  sample_frac = 0.5
)

meta_cols_score <- rev(names(cols_scores)[c(11,10)])
plot_density_with_peak_general(
  seurat_obj = TMA_merged_sub,
  features = meta_cols_score,
  color_vec = cols_scores,
  data_type = "assay",
  assay_name = "AUCell",
  sample_frac = 0.5
)

corr_2 <- readRDS("/cluster3/yflu/STS/Drug_screen/corr_2.rds")

mat <- corr_2[, -1]
rownames(mat) <- corr_2$Genes

# 取交集基因
common_genes <- intersect(rownames(mat), rownames(TMA_merged))

# 保留交集后的矩阵
mat_sub <- mat[common_genes, ]

origins <- names(cols_scores)[c(7,5,1,9)]
origin_genes <- rownames(mat_sub)[order(mat_sub[,origins[1]],decreasing = T)[1]]
for (i in 2:length(origins)) {
  origin_genes <- c(origin_genes,rownames(mat_sub)[order(mat_sub[,origins[i]],decreasing = T)[1]]) 
}

plot_custom_feature_scaled(TMA_merged_sub,
                           features = origin_genes,
                           data_type = "assay",
                           assay_name = "RNA",
                           slot = "data",
                           reduction = "spatial",
                           color_palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),point_size = 0.2)

Average_markers <- AverageExpression(TMA_merged_tumor,features = common_genes,group.by = "Sample",assays = "RNA")
Average_markers <- as.data.frame(Average_markers$RNA)
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

library(ggplot2)
library(dplyr)
library(tidyr)

plot_gene_expression_boxplot <- function(expr_df, sample_anno, disease, genes,
                                         disease_order, disease_colors,
                                         test = "t", sort_by_p = FALSE) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # -----------------------
  # 1. 样本分组
  # -----------------------
  sample_anno <- sample_anno %>%
    mutate(Group = ifelse(Diseases == disease, disease, paste0("Not_", disease)))
  
  group_samples <- as.character(sample_anno$Samples)
  
  # -----------------------
  # 2. 筛选表达矩阵子集
  # -----------------------
  expr_sub <- expr_df[genes, group_samples, drop=FALSE]
  expr_sub_t <- as.data.frame(t(expr_sub))
  expr_sub_t$Group <- sample_anno$Group[match(rownames(expr_sub_t), sample_anno$Samples)]
  expr_sub_t$Diseases <- sample_anno$Diseases[match(rownames(expr_sub_t), sample_anno$Samples)]
  
  # -----------------------
  # 3. 统计检验
  # -----------------------
  stats_res <- sapply(genes, function(g) {
    group1 <- expr_sub_t %>% filter(Group == disease) %>% pull(g)
    group2 <- expr_sub_t %>% filter(Group == paste0("Not_", disease)) %>% pull(g)
    if (test=="t") {
      t.test(group1, group2)$p.value
    } else if (test=="wilcox") {
      wilcox.test(group1, group2)$p.value
    } else {
      stop("test must be 't' or 'wilcox'")
    }
  })
  
  # -----------------------
  # 4. 准备绘图数据
  # -----------------------
  plot_df <- expr_sub_t %>%
    select(all_of(genes), Group, Diseases) %>%
    pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")
  
  # Gene factor 顺序
  if (sort_by_p) {
    plot_df$Gene <- factor(plot_df$Gene, levels = names(sort(stats_res)))
  } else {
    plot_df$Gene <- factor(plot_df$Gene, levels = genes)
  }
  
  # Group factor 顺序
  plot_df$Group <- factor(plot_df$Group, levels = c(disease, paste0("Not_", disease)))
  
  # -----------------------
  # 5. 显著性标注
  # -----------------------
  get_sig <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  
  sig_df <- data.frame(
    Gene = factor(names(stats_res), levels = levels(plot_df$Gene)),
    y = sapply(names(stats_res), function(g) max(plot_df$Expression[plot_df$Gene==g])*1.15),
    label = sapply(stats_res, get_sig)
  )
  
  # -----------------------
  # 6. 配色
  # -----------------------
  names(disease_colors) <- disease_order
  disease_color_map <- disease_colors
  
  # 箱线图颜色：指定disease vs 其他
  box_fill <- c(disease_color_map[disease], "grey70")
  names(box_fill) <- c(disease, paste0("Not_", disease))
  
  dodge_width <- 0.8
  
  # -----------------------
  # 7. 绘图
  # -----------------------
  p <- ggplot(plot_df, aes(x = Gene, y = Expression)) +
    geom_boxplot(aes(fill = Group),
                 outlier.shape = NA, alpha = 0.6,
                 position = position_dodge(width = dodge_width)) +
    geom_jitter(aes(fill = Diseases, group = Group),
                shape = 21,
                color = "gray",   # 外框黑色
                stroke = 0.3,
                size = 2,
                position = position_jitterdodge(jitter.width = 2, dodge.width = dodge_width),
                alpha = 0.8) +
    geom_segment(data = sig_df,
                 aes(x = as.numeric(Gene) - 0.2, xend = as.numeric(Gene) + 0.2,
                     y = y, yend = y),
                 inherit.aes = FALSE) +
    geom_text(data = sig_df,
              aes(x = Gene, y = y*1.02, label = label),
              inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c(box_fill, disease_color_map)) +  # box + 点颜色都能覆盖
    theme_classic() +
    labs(title = paste("Expression of", paste(genes, collapse=", "), "in", disease, "vs others"))
  
  print(p)
  
  return(stats_res)
}

plot_gene_expression_boxplot(Average_markers, 
                             sample_anno, 
                             disease = "RMS", 
                             genes = origin_genes,
                             disease_order = disease_order,
                             disease_colors = disease_colors,
                             test="t")

plot_gene_expression_boxplot_multi_panel <- function(expr_df, sample_anno, diseases, genes,
                                                     disease_order, disease_colors,
                                                     group_name = NULL, group_color = "#E64B35",
                                                     test = "t", sort_by_p = FALSE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  if(is.null(group_name)){
    group_name <- if(length(diseases)==1) diseases else paste(diseases, collapse="+")
  }
  
  # -----------------------
  # 1. 构建 BoxGroup: 单独疾病 + Other
  # -----------------------
  plot_df <- expr_df[genes, , drop=FALSE] %>% t() %>% as.data.frame()
  plot_df$Sample <- rownames(plot_df)
  plot_df$BoxGroup <- case_when(
    sample_anno$Diseases[match(plot_df$Sample, sample_anno$Samples)] %in% diseases ~ 
      sample_anno$Diseases[match(plot_df$Sample, sample_anno$Samples)],
    TRUE ~ "Other"
  )
  plot_df$Diseases <- sample_anno$Diseases[match(plot_df$Sample, sample_anno$Samples)]
  
  # -----------------------
  # 2. 添加 Merged 组
  # -----------------------
  merged_df <- expr_df[genes, , drop=FALSE] %>% t() %>% as.data.frame()
  merged_df$Sample <- rownames(merged_df)
  merged_df$BoxGroup <- ifelse(sample_anno$Diseases[match(merged_df$Sample, sample_anno$Samples)] %in% diseases,
                               group_name, NA)
  merged_df$Diseases <- sample_anno$Diseases[match(merged_df$Sample, sample_anno$Samples)]
  merged_df <- merged_df[!is.na(merged_df$BoxGroup), ]
  
  plot_df <- bind_rows(plot_df, merged_df)
  
  # -----------------------
  # 3. pivot_longer
  # -----------------------
  plot_df <- plot_df %>%
    pivot_longer(cols=all_of(genes), names_to="Gene", values_to="Expression")
  
  # Gene factor 顺序
  if(sort_by_p){
    stats_res <- sapply(genes, function(g){
      group1 <- plot_df %>% filter(BoxGroup==group_name, Gene==g) %>% pull(Expression)
      group2 <- plot_df %>% filter(BoxGroup=="Other", Gene==g) %>% pull(Expression)
      if(test=="t") t.test(group1, group2)$p.value else wilcox.test(group1, group2)$p.value
    })
    plot_df$Gene <- factor(plot_df$Gene, levels=names(sort(stats_res)))
  } else {
    plot_df$Gene <- factor(plot_df$Gene, levels=genes)
  }
  
  # -----------------------
  # 4. BoxGroup 顺序
  # -----------------------
  plot_df$BoxGroup <- factor(plot_df$BoxGroup, levels=c(diseases, group_name, "Other"))
  
  # -----------------------
  # 5. 统计显著性
  # -----------------------
  sig_list <- list()
  for(g in genes){
    # Merged vs Other
    group1 <- plot_df %>% filter(BoxGroup==group_name, Gene==g) %>% pull(Expression)
    group2 <- plot_df %>% filter(BoxGroup=="Other", Gene==g) %>% pull(Expression)
    p_merged <- if(test=="t") t.test(group1, group2)$p.value else wilcox.test(group1, group2)$p.value
    sig_list[[paste0(g,"_",group_name)]] <- data.frame(Gene=g, Group1=group_name, Group2="Other",
                                                       y=max(plot_df$Expression[plot_df$Gene==g])*1.15,
                                                       label=if(p_merged<0.001) "***" else if(p_merged<0.01) "**" else if(p_merged<0.05) "*" else "ns")
    # 单独疾病 vs Other
    for(d in diseases){
      group1 <- plot_df %>% filter(BoxGroup==d, Gene==g) %>% pull(Expression)
      group2 <- plot_df %>% filter(BoxGroup=="Other", Gene==g) %>% pull(Expression)
      p_val <- if(test=="t") t.test(group1, group2)$p.value else wilcox.test(group1, group2)$p.value
      sig_list[[paste0(g,"_",d)]] <- data.frame(Gene=g, Group1=d, Group2="Other",
                                                y=max(plot_df$Expression[plot_df$Gene==g])*1.05,
                                                label=if(p_val<0.001) "***" else if(p_val<0.01) "**" else if(p_val<0.05) "*" else "ns")
    }
  }
  sig_df <- bind_rows(sig_list)
  
  # -----------------------
  # 6. 配色
  # -----------------------
  names(disease_colors) <- disease_order
  box_fill <- c(setNames(disease_colors[diseases], diseases),
                setNames(group_color, group_name), Other="grey70")
  
  # -----------------------
  # 7. 绘图
  # -----------------------
  # 增大 box dodge 宽度、减小每个 box 自身宽度
  dodge_width <- 1
  box_width <- 0.9
  
  p <- ggplot(plot_df, aes(x=Gene, y=Expression)) +
    geom_boxplot(aes(fill=BoxGroup), outlier.shape=NA, alpha=0.6, color="black",
                 position=position_dodge(width=dodge_width), width=box_width) +
    geom_jitter(aes(fill=Diseases, group=BoxGroup),
                shape=21, color="gray", stroke=0.3, size=2,
                position=position_jitterdodge(jitter.width=0.8, dodge.width=dodge_width), alpha=0.8) +
    geom_segment(data=sig_df,
                 aes(x=Gene, xend=Gene, y=y, yend=y),
                 inherit.aes=FALSE, color="black", size=0.8) +
    geom_text(data=sig_df,
              aes(x=Gene, y=y*1.02, label=label),
              inherit.aes=FALSE, size=4, color="black") +
    scale_fill_manual(values=c(box_fill, disease_colors)) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x="", title=paste("Expression of", paste(genes, collapse=", "),
                           "in selected diseases vs others"))
  
  print(p)
  
  return(sig_df)
}

plot_gene_expression_boxplot_multi_panel(
  expr_df = Average_markers, sample_anno,
  genes  = c("ENG","TGFBR2"),
  diseases = c("LYM","HE", "KHE","AS"),
  group_name = "Endothelial_merged",
  group_color = "#BE8FBE",
  disease_order = disease_order,
  disease_colors = disease_colors,
)

plot_density_with_peak_violin_highlight <- function(
    seurat_obj,
    features,
    color_vec,
    data_type = c("meta", "assay"),
    assay_name = "RNA",
    sample_frac = 0.1,
    highlight_frac = 0.5,
    seed = 123,
    line_height = 0.3,
    sample_col,
    highlight_samples,
    highlight_color_vec,
    highlight_alpha = 0.3,
    highlight_peak_color = "red",
    x_range = NULL       # 新增：限制 x 轴范围，例如 c(-1, 3)
) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  data_type <- match.arg(data_type)
  set.seed(seed)
  
  # 获取数据
  if (data_type == "meta") {
    meta_df <- seurat_obj@meta.data
    df <- meta_df[, c(features, sample_col), drop = FALSE]
    df$cell_id <- rownames(df)
    data_long <- melt(df, id.vars = c("cell_id", sample_col),
                      variable.name = "feature", value.name = "value")
  } else {
    assay_data <- GetAssayData(seurat_obj, assay = assay_name, slot = "data")[features, , drop = FALSE]
    df <- as.data.frame(t(as.matrix(assay_data)))
    df$cell_id <- rownames(df)
    df[[sample_col]] <- seurat_obj@meta.data[[sample_col]]
    data_long <- melt(df, id.vars = c("cell_id", sample_col),
                      variable.name = "feature", value.name = "value")
  }
  
  # 标记高亮
  data_long$highlight_flag <- data_long[[sample_col]] %in% highlight_samples
  
  # 分别抽样
  data_high <- data_long[data_long$highlight_flag, ] %>%
    group_by(feature) %>% sample_frac(highlight_frac) %>% ungroup()
  data_nonhigh <- data_long[!data_long$highlight_flag, ] %>%
    group_by(feature) %>% sample_frac(sample_frac) %>% ungroup()
  data_sampled <- bind_rows(data_high, data_nonhigh)
  
  data_sampled$feature <- factor(data_sampled$feature, levels = features)
  
  # 创建 fill_color
  data_sampled$fill_color <- ifelse(
    data_sampled$highlight_flag,
    highlight_color_vec[as.character(data_sampled[[sample_col]])],
    color_vec[as.character(data_sampled$feature)]
  )
  
  # 计算峰值
  peak_normal <- data_sampled %>%
    filter(!highlight_flag) %>%
    group_by(feature) %>%
    summarise(
      peak = {
        v <- value[!is.na(value)]
        if(length(v) >= 2) density(v)$x[which.max(density(v)$y)]
        else if(length(v) == 1) v else NA_real_
      },
      ycenter = as.numeric(feature),
      .groups = "drop"
    ) %>% filter(!is.na(peak))
  
  peak_highlight <- data_sampled %>%
    filter(highlight_flag) %>%
    group_by(feature) %>%
    summarise(
      peak = {
        v <- value[!is.na(value)]
        if(length(v) >= 2) density(v)$x[which.max(density(v)$y)]
        else if(length(v) == 1) v else NA_real_
      },
      ycenter = as.numeric(feature),
      .groups = "drop"
    ) %>% filter(!is.na(peak))
  
  # 绘图
  p <- ggplot() +
    geom_violin(data = data_sampled[!data_sampled$highlight_flag, ],
                aes(x = value, y = feature, fill = fill_color),
                scale = "width", trim = TRUE, alpha = 1, color = NA) +
    geom_violin(data = data_sampled[data_sampled$highlight_flag, ],
                aes(x = value, y = feature, fill = fill_color),
                scale = "width", trim = TRUE, alpha = highlight_alpha, color = NA) +
    geom_segment(data = peak_normal,
                 aes(x = peak, xend = peak,
                     y = ycenter - line_height,
                     yend = ycenter + line_height),
                 color = "black", size = 0.8) +
    geom_segment(data = peak_highlight,
                 aes(x = peak, xend = peak,
                     y = ycenter - line_height,
                     yend = ycenter + line_height),
                 color = highlight_peak_color, size = 0.8) +
    scale_fill_identity() +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          legend.position = "none") +
    xlab("value")
  
  # 限制 x 轴范围（可选）
  if (!is.null(x_range)) {
    p <- p + coord_cartesian(xlim = x_range)
  }
  
  return(p)
}

cols_highlight <- c("T1620" = "#ADD8E6")

plot_density_with_peak_violin_highlight(
  seurat_obj = TMA_merged_sub,
  features = rev(names(cols_1)),
  color_vec = cols_1,
  data_type = "meta",
  sample_frac = 0.1,
  sample_col = "Sample",
  highlight_samples = c("T1620"),
  highlight_alpha = 0.6,
  highlight_frac = 0.1,       # 高亮组透明度
  highlight_peak_color = "red",
  highlight_color_vec = cols_highlight,
  x_range = c(0, 40)
)
