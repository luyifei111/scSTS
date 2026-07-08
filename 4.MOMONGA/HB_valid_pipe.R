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

HB_tumor <- readRDS("/cluster3/yflu/STS/Drug_screen/HB/9celltype_HB.rds")
samplenames <- unique(HB_tumor$sample)

cnv_genes <- readRDS("/cluster3/yflu/STS/Drug_screen/HB/cnv_genes.rds")
drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))

cnv_genes_top_list <- list()
CNV_correlated_drugs_list <- list()

degs_tumor_normal_all <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_tumor_markers.rds")
degs_tumor_normal_all <- subset(degs_tumor_normal_all,avg_log2FC > 0&p_val_adj < 0.05&cluster == "Tumor")

regulon_targets <- read.csv("/cluster3/yflu/STS/Drug_screen/HB/scenic/reg.csv")
colnames(regulon_targets) <- c(regulon_targets[2,c(1,2)],regulon_targets[1,-c(1,2)])
regulon_targets <- regulon_targets[-c(1,2),]

#HB_sample_markers <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_sample_markers.rds")

for (i in 1:length(samplenames)) {
  cnv_genes_sub <- subset(cnv_genes,sample == samplenames[i])
  cnv_genes_gain <- subset(cnv_genes_sub,status == "gain")
  cnv_genes_loss <- subset(cnv_genes_sub,status == "loss")
  #markers_up <- subset(NB_sample_markers,cluster == samplenames[i]&avg_log2FC > 1&p_val_adj < 0.05&pct.1 > 0.5)
  #markers_down <- subset(NB_sample_markers,cluster == samplenames[i]&avg_log2FC < -1&p_val_adj < 0.05&pct.2 > 0.5)
  #markers_up_top <- markers_up$gene
  #markers_down_top <- markers_down$gene
  #cnv_genes_gain_top <- subset(cnv_genes_gain,genes %in% markers_up_top)
  #cnv_genes_loss_top <- subset(cnv_genes_loss,genes %in% markers_down_top)
  cnv_genes_gain_top <- cnv_genes_gain
  cnv_genes_loss_top <- cnv_genes_loss
  cnv_genes_top <- rbind(cnv_genes_gain_top,cnv_genes_loss_top)
  cnv_genes_top <- cnv_genes_top[!duplicated(cnv_genes_top$CNV),]
  
  if(length(intersect(cnv_genes_top$genes,regulon_targets$TF)) > 0){
    regulon_targets_sub <- subset(regulon_targets,TF %in% cnv_genes_top$genes)
    TFs <- unique(regulon_targets_sub$TF)
    TF_tragets <- c()
    
    k=TFs[1]
    regulon_targets_sub_1 <- subset(regulon_targets_sub,TF == k)
    regulon_targets_sub_1_long <- regulon_targets_sub_1 %>%
      mutate(row_id = row_number()) %>%  # 保留行号
      group_by(row_id, TF = TF) %>%
      group_modify(~ {
        df <- parse_targetgenes(.x$TargetGenes)
        if(is.null(df)) return(NULL)
        df
      }) %>%
      ungroup()
    regulon_targets_sub_1_targets <- unique(regulon_targets_sub_1_long$Gene)
    #regulon_targets_sub_1_targets <- intersect(markers_up_top,regulon_targets_sub_1_targets)
    if(length(regulon_targets_sub_1_targets) > 0){
      regulon_targets_sub_1_targets <- data.frame(
        CNV = paste(k, "target", sep = "_"),
        sample = samplenames[i],
        status = "TF target",
        genes = regulon_targets_sub_1_targets,
        stringsAsFactors = FALSE
      )
    } else {
      regulon_targets_sub_1_targets <- data.frame(
        CNV = paste(k, "target", sep = "_"),
        sample = samplenames[i],
        status = "TF target",
        genes = "0",
        stringsAsFactors = FALSE
      )
    }
    for (k in TFs[-1]) {
      regulon_targets_sub_1 <- subset(regulon_targets_sub,TF == k)
      regulon_targets_sub_1_long <- regulon_targets_sub_1 %>%
        mutate(row_id = row_number()) %>%  # 保留行号
        group_by(row_id, TF = TF) %>%
        group_modify(~ {
          df <- parse_targetgenes(.x$TargetGenes)
          if(is.null(df)) return(NULL)
          df
        }) %>%
        ungroup()
      regulon_targets_sub_1_targets_1 <- unique(regulon_targets_sub_1_long$Gene)
      #regulon_targets_sub_1_targets_1 <- intersect(markers_up_top,regulon_targets_sub_1_targets_1)
      if(length(regulon_targets_sub_1_targets_1) > 0){
        regulon_targets_sub_1_targets_1 <- data.frame(
          CNV = paste(k, "target", sep = "_"),
          sample = samplenames[i],
          status = "TF target",
          genes = regulon_targets_sub_1_targets_1,
          stringsAsFactors = FALSE
        )
        regulon_targets_sub_1_targets <- rbind(regulon_targets_sub_1_targets,regulon_targets_sub_1_targets_1)
      }
    }
    cnv_genes_top <- rbind(cnv_genes_top,regulon_targets_sub_1_targets)
  }
  cnv_genes_top <- subset(cnv_genes_top,genes %in% degs_tumor_normal_all$gene)
  cnv_genes_top_list[[i]] <- cnv_genes_top
  names(cnv_genes_top_list)[i] <- samplenames[i]
  
  CNV_correlated_drugs <- c()
  for (j in 1:length(drug_hgnc_list_combined)) {
    if(length(intersect(drug_hgnc_list_combined[[j]],subset(cnv_genes_top,status %in% c("gain","TF target"))$genes) > 0)){
      CNV_correlated_drugs <- c(CNV_correlated_drugs,names(drug_hgnc_list_combined)[j])
    }
  }
  if(length(CNV_correlated_drugs)==0){
    CNV_correlated_drugs <- c(0)
  }
  
  CNV_correlated_drugs_list[[i]] <- CNV_correlated_drugs
  names(CNV_correlated_drugs_list)[i] <- samplenames[i]
  print(paste(i,"/33",sep = ""))
}
saveRDS(cnv_genes_top_list,"/cluster3/yflu/STS/Drug_screen/HB_cnv_genes_top_list_251204.rds")
saveRDS(CNV_correlated_drugs_list,"/cluster3/yflu/STS/Drug_screen/HB_CNV_correlated_drugs_list_251204.rds")

aucell_druggable_genesets <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_aucell_druggable_genesets.rds")
aucell_druggable_genesets <- t(as.matrix(aucell_druggable_genesets))
aucell_druggable_genesets_scaled <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_aucell_druggable_genesets_scaled.rds")

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

samplenames_new <- gsub("_", "-", samplenames)

CNV_aucell_drugs_top_list <- list()
for (i in 1:length(samplenames)) {
  if(CNV_correlated_drugs_list[[i]][1] == c(0)){
    CNV_aucell_drugs_top <- c(0)
  } else {
    CNV_correlated_drugs_sub <- CNV_correlated_drugs_list[[i]]
    CNV_aucell_drugs_top <- as.data.frame(aucell_druggable_genesets[samplenames_new[i],CNV_correlated_drugs_sub])
    colnames(CNV_aucell_drugs_top) <- "aucell_druggable_genesets"
    rownames(CNV_aucell_drugs_top) <- CNV_correlated_drugs_sub
    CNV_aucell_drugs_top$drugnames <- rownames(CNV_aucell_drugs_top)
    for (j in 1:nrow(CNV_aucell_drugs_top)) {
      CNV_aucell_drugs_top$drugnames[j] <-  paste(strsplit(rownames(CNV_aucell_drugs_top)[j],"[.]")[[1]][1],collapse = " ")
    }
    CNV_aucell_drugs_top <- CNV_aucell_drugs_top[order(CNV_aucell_drugs_top$aucell_druggable_genesets,decreasing = T),]
    if(nrow(CNV_aucell_drugs_top) > 20){
      CNV_aucell_drugs_top <- subset(CNV_aucell_drugs_top,aucell_druggable_genesets >= CNV_aucell_drugs_top$aucell_druggable_genesets[20])
    }
    drugs_kept <- c()
    for (j in 1:nrow(CNV_aucell_drugs_top)) {
      CNV_D2C_scale <- as.data.frame(aucell_druggable_genesets_scaled[,rownames(CNV_aucell_drugs_top)[j]])
      colnames(CNV_D2C_scale) <- "drug"
      CNV_D2C_scale$sample <- rownames(CNV_D2C_scale)
      CNV_D2C_scale_ordered <- as.data.frame(CNV_D2C_scale[order(CNV_D2C_scale$drug,decreasing = T),])
      if(samplenames_new[i] %in% subset(CNV_D2C_scale_ordered,drug >= CNV_D2C_scale_ordered$drug[as.numeric(round(quantile(1:length(samplenames))[2]))])$sample){
        drugs_kept <- c(drugs_kept,CNV_aucell_drugs_top$drugnames[j])
      }
    }
    if(length(drugs_kept) == 0){
      CNV_aucell_drugs_top <- c(0)
    } else {
      CNV_aucell_drugs_top <- subset(CNV_aucell_drugs_top,drugnames %in% drugs_kept)
    } 
  }
  CNV_aucell_drugs_top_list[[i]] <- CNV_aucell_drugs_top
  names(CNV_aucell_drugs_top_list)[i] <- samplenames[i]
  print(paste(i,"/33",sep = ""))
}

saveRDS(CNV_aucell_drugs_top_list,"HB_CNV_aucell_drugs_top_list_251204.rds")

HB_CNV_aucell_drugs_top_list_251204 <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_CNV_aucell_drugs_top_list_251204.rds")
HB_cnv_genes_top_list_251204 <- readRDS("/cluster3/yflu/STS/Drug_screen/HB_cnv_genes_top_list_251204.rds")

drug_count <- list()   # 保存每个样本的计数
drug_global_count <- list()  # 所有样本的总体计数

for (nm in names(HB_CNV_aucell_drugs_top_list_251204)) {
  
  item <- HB_CNV_aucell_drugs_top_list_251204[[nm]]
  
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

### 2. top5
top10_rows <- df[1:10, ]

### 3. inferno 配色（反转，使高频颜色变亮）
n_colors <- nrow(df)
inferno_colors <- rev(inferno(n_colors))

### 4. 绘图
ggplot(df, aes(x = Freq, y = factor(Var1, levels = rev(Var1)))) +
  geom_col(fill = inferno_colors, width = 1) +  
  scale_y_discrete(expand = c(0, 0)) +           
  
  geom_text_repel(
    data = top10_rows,
    aes(label = Var1),
    xlim = c(max(df$Freq) * 1.1, NA),
    direction = "y",
    hjust = 0,
    size = 4,
    box.padding = 0.4,
    min.segment.length = 0,
    segment.color = "grey40"
  ) +
  
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 60, 5, 5),
    panel.grid = element_blank()
  ) +
  labs(x = "Freq", y = "")

drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))

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
CNV_aucell_drugs_merged <- combine_aucell_drugs_keep_rownames_full(HB_CNV_aucell_drugs_top_list_251204)
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
CNV_aucell_drugs_merged_sub_1 <- CNV_aucell_drugs_merged_sub[,-c(1,35)]

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
    if (!is.null(HB_cnv_genes_top_list_251204[[sample]])) {
      sample_genes <- HB_cnv_genes_top_list_251204[[sample]]$genes
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
    if (sample_id %in% names(HB_cnv_genes_top_list_251204)) {
      cnv_df <- HB_cnv_genes_top_list_251204[[sample_id]]
      
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

top10_label <- top_gene_per_drug2_unique %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 10)
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
