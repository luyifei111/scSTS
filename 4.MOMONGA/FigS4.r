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
library(dplyr)
library(purrr)
library(tibble)

CNV_disease_durgs_sub <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_CNV_disease_durgs_sub_251203.rds")
CNV_aucell_drugs_top_list <- readRDS("/cluster3/yflu/STS/Drug_screen/CNV_aucell_drugs_top_list_251203.rds")
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

full_names <- c("Hemangioma", "KHE", "Schwannoma", "MPNST", "Undifferentiated sarcoma",
                "RMS", "MRT", "IMT", "Angiosarcoma", "EWS/PNET",
                "NF", "Aggressive fibromatosis", "Liposarcoma", "Spindle cell tumor", "ASPS",
                "Infantile fibrosarcoma", "Synovial sarcoma", "Lipoblastoma", "Pecoma", "Lymphangioma")

# 对应缩写
abbreviations <- c("HE","KHE","SWN","MPNST","US","RMS","MRT","IMT","AS","EWS",
                   "NF","AF","LPS","SCT","ASPS","IFS","SS","LPB","PECOMA","LYM")

# 创建替换映射
name_map <- setNames(abbreviations, full_names)
anno_sample_cluster_extended$Disease <- name_map[anno_sample_cluster_extended$Disease]

get_drug_names <- function(x) {
  
  ## 1. NULL / NA
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(character(0))
  }
  
  ## 2. 空样本：numeric 0（你给的例子）
  if (is.numeric(x) && length(x) == 1 && x == 0) {
    return(character(0))
  }
  
  ## 3. matrix / data.frame：药物名在 rownames
  rn <- rownames(x)
  if (!is.null(rn) && length(rn) > 0) {
    return(as.character(rn))
  }
  
  ## 4. named vector：药物名在 names
  nm <- names(x)
  if (!is.null(nm) && length(nm) > 0) {
    return(as.character(nm))
  }
  
  ## 5. 其他情况：当作无敏感药物
  character(0)
}

library(dplyr)
library(purrr)
library(tibble)

## 1. 过滤有效样本（只保留有疾病药物注释的）
valid_diseases <- unique(CNV_disease_durgs_sub$Disease)

samples_use <- anno_sample_cluster_extended %>%
  rownames_to_column("Sample") %>%
  filter(Disease %in% valid_diseases)

## 2. Disease -> 药物集合（强制 character，解决 factor）
disease2drugs <- CNV_disease_durgs_sub %>%
  mutate(Drugs = as.character(Drugs)) %>%
  group_by(Disease) %>%
  summarise(disease_drugs = list(unique(Drugs)), .groups = "drop")

## 3. 核心计算
res <- samples_use %>%
  left_join(disease2drugs, by = "Disease") %>%
  mutate(
    ## 样本敏感药物（来自 rownames / names，0 -> 空）
    sensitive_drugs = map(
      Sample,
      ~ get_drug_names(CNV_aucell_drugs_top_list[[.x]])
    ),
    
    ## 疾病对应药物
    disease_drugs = map(
      disease_drugs,
      ~ if (is.null(.x) || length(.x) == 0 || all(is.na(.x)))
        character(0)
      else
        as.character(.x)
    ),
    ## overlap / non-overlap
    overlap_drugs    = map2(sensitive_drugs, disease_drugs, intersect),
    nonoverlap_drugs = map2(sensitive_drugs, disease_drugs, setdiff),
    ## 数量
    n_sensitive  = map_int(sensitive_drugs, length),
    n_overlap    = map_int(overlap_drugs, length),
    n_nonoverlap = map_int(nonoverlap_drugs, length),
    ## 比例（无敏感药物 → 0）
    overlap_ratio    = if_else(n_sensitive == 0, 0, n_overlap / n_sensitive),
    nonoverlap_ratio = if_else(n_sensitive == 0, 0, n_nonoverlap / n_sensitive)
)

disease_levels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS",
                    "NF","SWN","LYM","HE","KHE","MPNST","US",
                    "AS","RMS","PECOMA","EWS","MRT")

res$Disease <- factor(res$Disease,levels = disease_levels)

res_plot <- res %>%
  arrange(Disease)   # 按 Disease label 排序
mat_ratio <- res_plot %>%
  dplyr::select(overlap_ratio, nonoverlap_ratio) %>%
  t()

colnames(mat_ratio) <- res_plot$Sample

anno_col <- data.frame(
  Disease = res_plot$Disease
)
rownames(anno_col) <- res_plot$Sample

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

annotation_colors <- list(
  Disease = cols[valid_diseases]
)

pheatmap(
  mat_ratio,
  cluster_rows = FALSE,
  cluster_cols = T,
  annotation_col = anno_col,
  show_colnames = FALSE,
  scale = "none",
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(0, 1, length.out = 101),
  fontsize_row = 12,
  border_color = NA
)

res_plot <- res %>%
  mutate(
    Disease = factor(Disease, levels = names(cols))
  )
res_plot_nonzero <- res_plot %>%
  filter(n_sensitive > 0)

library(ggplot2)

p <- ggplot(
  res_plot_nonzero,
  aes(x = Disease, y = overlap_ratio, fill = Disease)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,
    alpha = 0.6
  ) +
  scale_fill_manual(values = cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Overlap ratio"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p

disease_n <- res_plot_nonzero %>%
  group_by(Disease) %>%
  summarise(
    n = sum(!is.na(overlap_ratio)),
    .groups = "drop"
  )

library(dplyr)
library(tibble)

# 关键：从数据本身取 unique，而不是 levels
diseases <- res_plot_nonzero %>%
  distinct(Disease) %>%
  pull(Disease) %>%
  as.character()

stat_pair <- combn(diseases, 2, simplify = FALSE) %>%
  lapply(function(x) {
    
    d1 <- x[1]
    d2 <- x[2]
    
    v1 <- res_plot_nonzero %>%
      filter(Disease == d1) %>%
      pull(overlap_ratio)
    
    v2 <- res_plot_nonzero %>%
      filter(Disease == d2) %>%
      pull(overlap_ratio)
    
    # 最低安全检查（理论上这里已经不需要，但保留更稳）
    if (length(v1) < 2 || length(v2) < 2) {
      return(tibble(group1 = d1, group2 = d2, p = NA_real_))
    }
    
    wt <- wilcox.test(v1, v2, exact = FALSE)
    
    tibble(
      group1 = d1,
      group2 = d2,
      p = wt$p.value
    )
  }) %>%
  bind_rows() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    p.adj.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

stat_pair

stat_sig <- stat_pair %>%
  filter(!is.na(p) & p < 0.05)

CNV_aucell_drugs_top_list_251203 <- readRDS("/cluster3/yflu/STS/Drug_screen/CNV_aucell_drugs_top_list_251203.rds")

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

# 使用
CNV_aucell_drugs_merged <- combine_aucell_drugs_keep_rownames_full(CNV_aucell_drugs_top_list_251203)
head(CNV_aucell_drugs_merged)

CNV_aucell_drugs_merged$nonzero_count <- apply(
  CNV_aucell_drugs_merged[ , setdiff(colnames(CNV_aucell_drugs_merged), "drug")],
  1,
  function(x) sum(x != 0 & !is.na(x))
)
CNV_aucell_drugs_merged <- as.data.frame(CNV_aucell_drugs_merged)
rownames(CNV_aucell_drugs_merged) <- CNV_aucell_drugs_merged$drug

CNV_aucell_drugs_merged_sub <- subset(CNV_aucell_drugs_merged,nonzero_count > 0)

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

cols_to_check <- setdiff(colnames(CNV_aucell_drugs_merged_sub), c("drug", "nonzero_count"))
# 每行非零列名保存为 list
row_nonzero_cols <- apply(CNV_aucell_drugs_merged_sub[, cols_to_check], 1, function(x) {
  names(x)[x != 0 & !is.na(x)]
})
# 转成 list（apply 返回的是 matrix 时可能自动转成 vector，所以加 as.list 保证 list）
row_nonzero_cols <- as.list(row_nonzero_cols)
names(row_nonzero_cols) <- CNV_aucell_drugs_merged_sub$drug
row_nonzero_cols

drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))

druggenes_list <- drug_hgnc_list_combined[CNV_aucell_drugs_merged_sub$drug]

drug_sample_gene_intersect <- list()

STS_cnv_genes_top_list_251202 <- readRDS("/cluster3/yflu/STS/Drug_screen/STS_cnv_genes_top_list_251202.rds")

for (drug in names(row_nonzero_cols)) {
  
  # 有值的样本
  samples_with_value <- row_nonzero_cols[[drug]]
  
  # drug 的 target genes
  drug_targets <- druggenes_list[[drug]]
  
  # 初始化该药物的结果
  intersect_list <- list()
  
  for (sample in samples_with_value) {
    # 检查 cnv_genes_top_list 中是否有该 sample
    if (!is.null(STS_cnv_genes_top_list_251202[[sample]])) {
      sample_genes <- STS_cnv_genes_top_list_251202[[sample]]$genes
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
    if (sample_id %in% names(STS_cnv_genes_top_list_251202)) {
      cnv_df <- STS_cnv_genes_top_list_251202[[sample_id]]
      
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

top_gene_per_drug <- extract_top_gene_per_drug(drug_gene_counts)

library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)

gene_anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(top_gene_per_drug$top_gene),
  keytype  = "SYMBOL",
  columns  = c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")
)

gene_anno <- gene_anno %>%
  mutate(
    start = abs(CHRLOC),
    end   = abs(CHRLOCEND)
  ) %>%
  dplyr::select(
    top_gene = SYMBOL,
    chr = CHR,
    start,
    end
  )

gene_anno_unique <- gene_anno %>%
  mutate(gene_length = end - start) %>%
  group_by(top_gene) %>%
  slice_max(gene_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-gene_length)

top_gene_per_drug2 <- top_gene_per_drug %>%
  left_join(gene_anno_unique, by = "top_gene")

top_gene_per_drug2$drug_short <- sub("\\..*$", "", top_gene_per_drug2$drug)

drug_count <- list()   # 保存每个样本的计数
drug_global_count <- list()  # 所有样本的总体计数

for (nm in names(CNV_aucell_drugs_top_list_251203)) {
  
  item <- CNV_aucell_drugs_top_list_251203[[nm]]
  
  # 跳过为 0 的元素
  if (is.numeric(item) && item == 0) {
    drug_count[[nm]] <- 0
    next
  }
  
  # 确保是 data.frame 且含 drugnames 列
  if (is.data.frame(item) && "drugnames" %in% colnames(item)) {
    
    # 统计当前样本 drug 数量
    drug_count[[nm]] <- table(item$drugnames)
    
    # 合并到全局计数
    drug_global_count <- c(drug_global_count, as.list(item$drugnames))
  }
}

# 全局 frequency table
drug_global_count <- as.data.frame(table(unlist(drug_global_count)))

library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridisLite)

### 0. 移除 Rle
df <- as.data.frame(drug_global_count)
df[] <- lapply(df, function(x) if (inherits(x, "Rle")) as.vector(x) else x)
df$Var1 <- as.character(df$Var1)
df$Freq <- as.numeric(df$Freq)

### 1. 排序（高→低）
df <- df %>% arrange(desc(Freq))

df$drug_short <- df$Var1

top_gene_per_drug2 <- top_gene_per_drug2 %>%
  left_join(df, by = "drug_short")

top_gene_per_drug2_unique <- top_gene_per_drug2 %>%
  group_by(drug_short) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup()

top_gene_per_drug2_unique <- top_gene_per_drug2_unique %>%
  group_by(top_gene) %>%
  slice_max(Freq, n = 1, with_ties = FALSE) %>%
  ungroup()

top_gene_per_drug2_unique <- top_gene_per_drug2_unique[,-1]

top_gene_per_drug2_unique <- top_gene_per_drug2_unique %>%
  mutate(
    chr = gsub("^chr", "", chr),
    chr = factor(
      chr,
      levels = c(as.character(1:22), "X", "Y")
    )
  )
top10_rows <- df[1:10, ]
top10_label <- top10_rows$Var1

genome_barplot_by_gene_inferno <- function(
    df,
    freq_col   = "Freq",
    chr_col    = "chr",
    start_col  = "start",
    end_col    = "end",
    label_col  = "drug_short",
    top_n      = 10,
    bar_width  = 5e6,      # 固定柱宽（bp）
    n_colors   = 100
) {
  
  library(dplyr)
  library(ggplot2)
  library(viridisLite)
  
  # ---- inferno 调色板（反向）----
  inferno_colors <- inferno(n_colors)
  
  # ---- 染色体长度（hg38）----
  chr_len <- tibble::tibble(
    chr = c(as.character(1:22), "X", "Y"),
    chr_length = c(
      248956422, 242193529, 198295559, 190214555, 181538259,
      170805979, 159345973, 145138636, 138394717, 133797422,
      135086622, 133275309, 114364328, 107043718, 101991189,
      90338345, 83257441, 80373285, 58617616, 64444167,
      46709983, 50818468, 156040895, 57227415
    )
  )
  
  # ---- chr 起点 ----
  chr_len <- chr_len %>%
    mutate(chr_start = cumsum(lag(chr_length, default = 0)))
  
  # ---- 映射到全基因组坐标 ----
  df_genome <- df %>%
    mutate(chr = gsub("^chr", "", .data[[chr_col]])) %>%
    left_join(chr_len, by = "chr") %>%
    mutate(
      genome_center = chr_start + (.data[[start_col]] + .data[[end_col]]) / 2
    )
  
  # ---- Top N 标注 ----
  df_label <- df_genome %>%
    arrange(desc(.data[[freq_col]])) %>%
    slice_head(n = top_n)
  
  # ---- chr 轴 ----
  chr_axis <- chr_len %>%
    mutate(chr_center = chr_start + chr_length / 2)
  
  # ---- 作图 ----
  p <- ggplot(df_genome) +
    geom_col(
      aes(
        x    = genome_center,
        y    = .data[[freq_col]],
        fill = .data[[freq_col]]
      ),
      width = bar_width
    ) +
    
    # Top N drug 标注
    geom_text(
      data = df_label,
      aes(
        x     = genome_center,
        y     = .data[[freq_col]],
        label = .data[[label_col]]
      ),
      angle = 30,
      vjust = -0.4,
      size  = 3
    ) +
    
    # 染色体分界线
    geom_vline(
      data = chr_len,
      aes(xintercept = chr_start),
      color = "grey80",
      linewidth = 0.3
    ) +
    
    # inferno 连续色标
    scale_fill_gradientn(
      colours = inferno_colors,
      name = freq_col
    ) +
    
    scale_x_continuous(
      breaks = chr_axis$chr_center,
      labels = chr_axis$chr
    ) +
    
    labs(
      x = "Genomic position",
      y = freq_col
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right"
    )
  
  return(p)
}
p <- genome_barplot_by_gene_inferno(
  df = top_gene_per_drug2_unique,
  freq_col  = "Freq",
  chr_col   = "chr",
  start_col = "start",
  end_col   = "end",
  label_col = "drug_short",
  top_n     = 10,
  bar_width = 1.2e7,
  n_colors  = 100
)
p

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
#samplenames <- subset(anno_sample_cluster_extended,Malignancy == "Malignant")
#samplenames <- rownames(samplenames)
samplenames <- rownames(anno_sample_cluster_extended)

cnv_genes <- readRDS("/cluster3/yflu/STS/development/target/cnv_genes.rds")

cnvobjpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/run.final.infercnv_obj",sep = "")

gene_order_file= read.delim("/cluster3/yflu/RT/inferCNV/hg38_gencode_v27.txt",header = F)

all_genes <- unlist(drug_hgnc_list_combined, use.names = FALSE)
all_genes <- unique(all_genes)

cnv_genes_intersect <- intersect(cnv_genes$genes,all_genes)
cnv_genes_intersect <- intersect(cnv_genes_intersect,gene_order_file$V1)
rownames(gene_order_file) <- gene_order_file$V1
cnv_genes_counts <- as.data.frame(table(cnv_genes$genes))
rownames(cnv_genes_counts) <- cnv_genes_counts$Var1
cnv_genes_counts <- cnv_genes_counts[cnv_genes_intersect,]

cnv_genes_counts <- cbind(cnv_genes_counts,gene_order_file[rownames(cnv_genes_counts),c(2:4)])

p <- genome_barplot_by_gene_inferno(
  df = cnv_genes_counts,
  freq_col  = "Freq",
  chr_col   = "V2",
  start_col = "V3",
  end_col   = "V4",
  label_col = "Var1",
  top_n     = 10,
  bar_width = 1.2e7,
  n_colors  = 100
)
p

## =============================
## Step 0. 基本一致性检查
## =============================
common_drugs <- intersect(
  names(drug_gene_counts),
  names(drug_hgnc_list_combined)
)

if (length(common_drugs) == 0) {
  stop("❌ drug_gene_counts 和 drug_hgnc_list_combined 没有共同的 drug 名称")
}

## 如果不完全一致，给提示
if (!setequal(names(drug_gene_counts), names(drug_hgnc_list_combined))) {
  warning("⚠️ 两个 list 的 drug 名称不完全一致，仅对交集部分处理")
}

## =============================
## Step 1. 逐 drug 处理
## =============================
drug_gene_counts_labeled <- lapply(
  common_drugs,
  function(drug) {
    
    df <- drug_gene_counts[[drug]]
    target_genes <- drug_hgnc_list_combined[[drug]]
    
    ## 防御式处理
    if (is.null(target_genes)) {
      df$Type <- "TF"
    } else {
      df$Type <- ifelse(
        df$all_genes %in% target_genes,
        "Target",
        "TF"
      )
    }
    
    df
  }
)

## 恢复 drug 名称
names(drug_gene_counts_labeled) <- common_drugs

## =============================
## Step 2. 简单检查
## =============================
# 看一个
head(drug_gene_counts_labeled[[1]])

# 看 Target / TF 分布
lapply(
  drug_gene_counts_labeled,
  function(x) table(x$Type)
)

## =============================
## Step 1. 逐 drug 判断总体类型
## =============================
drug_level_type <- sapply(
  drug_gene_counts_labeled,
  function(df) {
    
    types <- unique(df$Type)
    
    if (length(types) == 1 && types == "Target") {
      "Target"
    } else if (length(types) == 1 && types == "TF") {
      "TF"
    } else {
      "Mixed"
    }
  }
)

## =============================
## Step 2. 保存为 data.frame
## =============================
drug_type_df <- data.frame(
  drug = names(drug_level_type),
  Type = unname(drug_level_type),
  row.names = NULL,
  stringsAsFactors = FALSE
)

## =============================
## Step 3. 快速检查
## =============================
head(drug_type_df)
table(drug_type_df$Type)

library(VennDiagram)
library(grid)

# 数量
target_only <- 65
tf_only     <- 27
both        <- 74

grid.newpage()
draw.pairwise.venn(
  area1 = target_only + both,   # Target 总数
  area2 = tf_only + both,       # TF 总数
  cross.area = both,            # 重叠（Mixed）
  category = c("Target", "TF"),
  fill = c("#FB8072", "#8DD3C7"),
  alpha = 0.6,
  lwd = 1.5,
  cat.cex = 1.2,
  cex = 1.2
)
CNV_aucell_drugs_merged_sub_tf <- subset(CNV_aucell_drugs_merged_sub,drug %in% subset(drug_type_df,Type == "TF")$drug)
CNV_aucell_drugs_merged_sub_tf <- CNV_aucell_drugs_merged_sub_tf[,-c(1,80)]
CNV_aucell_drugs_merged_sub_tf <- CNV_aucell_drugs_merged_sub_tf[ , colSums(CNV_aucell_drugs_merged_sub_tf != 0) > 0 ]

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

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

colors_combined = colorRampPalette(brewer.pal(8,'RdBu'))(500)

pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_tf,
                   cluster_rows = T,
                   color = rev(c(colors_combined[1:200],"#FFFFFF")),
                   annotation_col = anno_sample_cluster_extended[,-3],
                   clustering_distance_cols = "manhattan",
                   annotation_colors = c(annotation_colors))

CNV_aucell_drugs_merged_sub_target <- subset(CNV_aucell_drugs_merged_sub,drug %in% subset(drug_type_df,Type == "Target")$drug)
CNV_aucell_drugs_merged_sub_target <- CNV_aucell_drugs_merged_sub_target[,-c(1,80)]
CNV_aucell_drugs_merged_sub_target <- CNV_aucell_drugs_merged_sub_target[ , colSums(CNV_aucell_drugs_merged_sub_target != 0) > 0 ]

pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_target,
                   cluster_rows = T,
                   color = rev(c(colors_combined[1:200],"#FFFFFF")),
                   annotation_col = anno_sample_cluster_extended[,-3],
                   clustering_distance_cols = "manhattan",
                   annotation_colors = c(annotation_colors))

CNV_aucell_drugs_merged_sub_mix <- subset(CNV_aucell_drugs_merged_sub,drug %in% subset(drug_type_df,Type == "Mixed")$drug)
CNV_aucell_drugs_merged_sub_mix <- CNV_aucell_drugs_merged_sub_mix[,-c(1,80)]
CNV_aucell_drugs_merged_sub_mix <- CNV_aucell_drugs_merged_sub_mix[ , colSums(CNV_aucell_drugs_merged_sub_mix != 0) > 0 ]
pheatmap::pheatmap(CNV_aucell_drugs_merged_sub_mix,
                   cluster_rows = T,
                   color = rev(c(colors_combined[1:200],"#FFFFFF")),
                   annotation_col = anno_sample_cluster_extended[,-3],
                   clustering_distance_cols = "manhattan",
                   annotation_colors = c(annotation_colors))

df_all <- do.call(
  rbind,
  lapply(
    names(drug_gene_counts_labeled),
    function(drug) {
      df <- drug_gene_counts_labeled[[drug]]
      df$drug <- drug
      df
    }
  )
)
gene_type_freq <- aggregate(
  Freq ~ all_genes + Type,
  data = df_all,
  sum
)
gene_type_freq_wide <- reshape(
  gene_type_freq,
  idvar = "all_genes",
  timevar = "Type",
  direction = "wide"
)

# 把 NA 变成 0
gene_type_freq_wide[is.na(gene_type_freq_wide)] <- 0

library(dplyr)
library(ggplot2)
library(tidyr)

## =============================
## Step 1. 整理成长表
## =============================
gene_type_long <- gene_type_freq_wide %>%
  dplyr::select(all_genes, Freq.Target, Freq.TF) %>%
  pivot_longer(
    cols = starts_with("Freq."),
    names_to = "Type",
    values_to = "Freq"
  ) %>%
  mutate(
    Type = gsub("Freq\\.", "", Type)
  )

## =============================
## Step 2. 分别取 Target / TF top10
## =============================
top10_target <- gene_type_long %>%
  filter(Type == "Target") %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 10)

top10_tf <- gene_type_long %>%
  filter(Type == "TF") %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 10)

## =============================
## Step 3. Target top10 柱状图（颜色 = counts）
## =============================
p_target <- ggplot(
  top10_target,
  aes(
    x = reorder(all_genes, Freq),
    y = Freq,
    fill = Freq
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(
    low  = "#FDE0DD",
    high = "#e4007f"
  ) +
  labs(
    title = "top 10 target genes (freq weighted)",
    x = NULL,
    y = "weighted frequency",
    fill = "counts"
  ) +
  theme_classic()

## =============================
## Step 4. TF top10 柱状图（颜色 = counts）
## =============================
p_tf <- ggplot(
  top10_tf,
  aes(
    x = reorder(all_genes, Freq),
    y = Freq,
    fill = Freq
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(
    low  = "#E0F3F1",
    high = "#036eb8"
  ) +
  labs(
    title = "top 10 TF genes (freq weighted)",
    x = NULL,
    y = "weighted frequency",
    fill = "counts"
  ) +
  theme_classic()

## =============================
## Step 5. 画图
## =============================
p_target
p_tf
