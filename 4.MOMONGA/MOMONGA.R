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

cnv_genes_top_list <- list()
CNV_correlated_drugs_list <- list()

degs_tumor_normal_all <- readRDS("/cluster3/yflu/STS/Drug_screen/degs_tumor_normal_all.rds")
degs_tumor_normal_all <- subset(degs_tumor_normal_all,avg_log2FC > 0&p_val_adj < 0.05)
degs_tumor_normal_all$gene <- rownames(degs_tumor_normal_all)

degs_tumor_normal_sample_list <- readRDS("/cluster3/yflu/STS/Drug_screen/degs_tumor_normal_sample_list.rds")

regulon_targets <- read.csv("/cluster3/yflu/STS/scenic/all_reg_cell_40.csv")
colnames(regulon_targets) <- c(regulon_targets[2,c(1,2)],regulon_targets[1,-c(1,2)])
regulon_targets <- regulon_targets[-c(1,2),]

for (i in 1:length(samplenames)) {
  cnv_genes_sub <- subset(cnv_genes,sample == samplenames[i])
  cnv_genes_gain <- subset(cnv_genes_sub,status == "gain")
  cnv_genes_loss <- subset(cnv_genes_sub,status == "loss")
  markers_up <- openxlsx::read.xlsx("/cluster3/yflu/STS/louvain_labels_channel_de20241230.xlsx",paste(samplenames[i],"|up",sep = ""))
  markers_down <- openxlsx::read.xlsx("/cluster3/yflu/STS/louvain_labels_channel_de20241230.xlsx",paste(samplenames[i],"|dn",sep = ""))
  markers_up_top <- subset(markers_up,log2FC >1)$featurekey
  markers_down_top <- subset(markers_down,log2FC < -1)$featurekey
  cnv_genes_gain_top <- subset(cnv_genes_gain,genes %in% markers_up_top)
  cnv_genes_loss_top <- subset(cnv_genes_loss,genes %in% markers_down_top)
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
    regulon_targets_sub_1_targets <- intersect(markers_up_top,regulon_targets_sub_1_targets)
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
      regulon_targets_sub_1_targets_1 <- intersect(markers_up_top,regulon_targets_sub_1_targets_1)
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
  }
  
  cnv_genes_top <- rbind(cnv_genes_top,regulon_targets_sub_1_targets)
  
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
  print(paste(i,"/78",sep = ""))
}

saveRDS(cnv_genes_top_list,"STS_cnv_genes_top_list_251202.rds")
saveRDS(CNV_correlated_drugs_list,"STS_CNV_correlated_drugs_list_251202.rds")

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

CNV_aucell_drugs_top_list <- list()
for (i in 1:length(samplenames)) {
  if(CNV_correlated_drugs_list[[i]][1] == c(0)){
    CNV_aucell_drugs_top <- c(0)
  } else {
    CNV_correlated_drugs_sub <- CNV_correlated_drugs_list[[i]]
    CNV_aucell_drugs_top <- as.data.frame(aucell_druggable_genesets[samplenames[i],CNV_correlated_drugs_sub])
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
      if(samplenames[i] %in% subset(CNV_D2C_scale_ordered,drug >= CNV_D2C_scale_ordered$drug[as.numeric(round(quantile(1:length(samplenames))[2]))])$sample){
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
  print(paste(i,"/78",sep = ""))
}

saveRDS(CNV_aucell_drugs_top_list,"CNV_aucell_drugs_top_list_251203.rds")
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
diseasenames <- unique(anno_sample_cluster_extended$Disease)

i = 1
disease_samplenames <- rownames(subset(anno_sample_cluster_extended,Disease == diseasenames[i]))
CNV_disease_durgs <- c()
disease_length <- length(disease_samplenames)

for (j in 1:disease_length) {
  if (ncol(as.data.frame(CNV_aucell_drugs_top_list[disease_samplenames[j]][[1]])) > 1){
    CNV_disease_durgs <- c(CNV_disease_durgs,rownames(CNV_aucell_drugs_top_list[disease_samplenames[j]][[1]]))
  }
}
CNV_disease_durgs <- as.data.frame(table(CNV_disease_durgs))
CNV_disease_durgs$prop <- CNV_disease_durgs$Freq/disease_length
CNV_disease_durgs$Disease <- rep(diseasenames[i],nrow(CNV_disease_durgs))
colnames(CNV_disease_durgs)[1] <- "Drugs"

for (i in 2:length(diseasenames)) {
  disease_samplenames <- rownames(subset(anno_sample_cluster_extended,Disease == diseasenames[i]))
  CNV_disease_durgs_1 <- c()
  disease_length <- length(disease_samplenames)
  
  for (j in 1:disease_length) {
    if (ncol(as.data.frame(CNV_aucell_drugs_top_list[disease_samplenames[j]][[1]])) > 1){
      CNV_disease_durgs_1 <- c(CNV_disease_durgs_1,rownames(CNV_aucell_drugs_top_list[disease_samplenames[j]][[1]]))
    }
  }
  CNV_disease_durgs_1 <- as.data.frame(table(CNV_disease_durgs_1))
  CNV_disease_durgs_1$prop <- CNV_disease_durgs_1$Freq/disease_length
  CNV_disease_durgs_1$Disease <- rep(diseasenames[i],nrow(CNV_disease_durgs_1))
  colnames(CNV_disease_durgs_1)[1] <- "Drugs"
  CNV_disease_durgs <- rbind(CNV_disease_durgs,CNV_disease_durgs_1)
}

CNV_disease_durgs_sub <- subset(CNV_disease_durgs,prop > 0.25)
CNV_disease_durgs_sub <- subset(CNV_disease_durgs_sub,Freq > 2)

CNV_disease_durgs_wide <- dcast(CNV_disease_durgs_sub,Drugs~Disease,value.var = 'Freq')
CNV_disease_durgs_wide[is.na(CNV_disease_durgs_wide)] <- 0
rownames(CNV_disease_durgs_wide) <- CNV_disease_durgs_wide$Drugs
CNV_disease_durgs_wide <- CNV_disease_durgs_wide[,-1]

CNV_disease_durgs.clust<-hclust(dist(CNV_disease_durgs_wide))
p2<-ggtree(CNV_disease_durgs.clust)
p2+
  geom_tiplab()+
  xlim(NA,7)

p1 <- ggplot(data = CNV_disease_durgs_sub,aes(x=Disease,y=Drugs))+geom_point(aes(size=prop,color=Freq))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightblue",high="blue")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+
  scale_y_discrete(position = "right")
p1
drugs_screened <- unique(CNV_disease_durgs_sub$Drugs)
write.csv(drugs_screened,"/cluster3/yflu/STS/Drug_screen/CNV_based_drugnames_anno_251203.csv")

# STS.integrated.pega_aucell_druggable_genesets <- readRDS("/cluster3/yflu/STS/Drug_screen/STS.integrated.pega_aucell_druggable_genesets.rds")
# druggable_logFC <- FindAllMarkers(STS.integrated.pega_aucell_druggable_genesets, assay = "AUCell",group.by = "Disease")
# druggable_average <- AverageExpression(STS.integrated.pega_aucell_druggable_genesets,assay = "AUCell",group.by = "Disease")
# druggable_average <- as.data.frame(druggable_average$AUCell)

#saveRDS(druggable_logFC,"druggable_logFC.rds")
#saveRDS(druggable_average,"druggable_average.rds")

druggable_logFC <- readRDS("/cluster3/yflu/STS/Drug_screen/druggable_logFC.rds")
druggable_average <- readRDS("/cluster3/yflu/STS/Drug_screen/druggable_average.rds")

druggable_average_scaled <- t(scale(t(as.matrix(druggable_average))))
druggable_average_scaled <- as.data.frame(druggable_average_scaled)
druggable_average_scaled$geneset <- rownames(druggable_average_scaled)

druggable_logFC$gene <- gsub("[\\| (),/-]", ".", druggable_logFC$gene)
druggable_average_scaled$geneset <- gsub("[\\| (),/-]", ".", druggable_average_scaled$geneset)

CNV_disease_durgs_sub$logFC <- CNV_disease_durgs_sub$Disease
CNV_disease_durgs_sub$av_exp <- CNV_disease_durgs_sub$Disease

for (i in 1:nrow(CNV_disease_durgs_sub)) {
  if(nrow(subset(druggable_logFC,cluster == CNV_disease_durgs_sub$Disease[i]&gene == CNV_disease_durgs_sub$Drugs[i]))>0){
    CNV_disease_durgs_sub$logFC[i] <- subset(druggable_logFC,cluster == CNV_disease_durgs_sub$Disease[i]&gene == CNV_disease_durgs_sub$Drugs[i])$avg_log2FC
  } else {
    CNV_disease_durgs_sub$logFC[i] <- "0"
  }
  CNV_disease_durgs_sub$av_exp[i] <- subset(druggable_average_scaled,geneset == CNV_disease_durgs_sub$Drugs[i])[,CNV_disease_durgs_sub$Disease[i]]
}

CNV_disease_durgs_sub

CNV_disease_durgs_sub$logFC <- as.numeric(CNV_disease_durgs_sub$logFC)
CNV_disease_durgs_sub$av_exp <- as.numeric(CNV_disease_durgs_sub$av_exp)

library(ggplot2)

full_names <- c("Hemangioma", "KHE", "Schwannoma", "MPNST", "Undifferentiated sarcoma",
                "RMS", "MRT", "IMT", "Angiosarcoma", "EWS/PNET",
                "NF", "Aggressive fibromatosis", "Liposarcoma", "Spindle cell tumor", "ASPS",
                "Infantile fibrosarcoma", "Synovial sarcoma", "Lipoblastoma", "Pecoma", "Lymphangioma")

# 对应缩写
abbreviations <- c("HE","KHE","SWN","MPNST","US","RMS","MRT","IMT","AS","EWS",
                   "NF","AF","LPS","SCT","ASPS","IFS","SS","LPB","PECOMA","LYM")

# 创建替换映射
name_map <- setNames(abbreviations, full_names)

# 替换 CNV_disease_durgs_sub$Disease
CNV_disease_durgs_sub$Disease <- name_map[CNV_disease_durgs_sub$Disease]

library(ggplot2)
library(RColorBrewer)

# Disease 缩写顺序
disease_levels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS","NF","SWN","LYM",
                    "HE","KHE","MPNST","US","AS","RMS","PECOMA","EWS","MRT")

# 对应自定义颜色
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

# 将 Disease 列转为因子，确保顺序正确
CNV_disease_durgs_sub$Disease <- factor(CNV_disease_durgs_sub$Disease, levels = disease_levels)

# 绘图示例（保留之前高亮逻辑）
CNV_disease_durgs_sub$plot_group <- with(CNV_disease_durgs_sub, 
                                         ifelse(logFC > 1 & av_exp > 2, "both",
                                                ifelse(logFC > 1 | av_exp > 2, "one", "none")))
# fill_factor 保留 Disease 色，none 为灰色
CNV_disease_durgs_sub$fill_factor <- ifelse(CNV_disease_durgs_sub$plot_group=="none", "Other", as.character(CNV_disease_durgs_sub$Disease))
CNV_disease_durgs_sub$fill_factor <- factor(CNV_disease_durgs_sub$fill_factor, levels = c(disease_levels, "Other"))

# 对应颜色向量，最后加灰色给 Other
fill_colors <- c(setNames(cols, disease_levels), Other="grey")

# 绘图
library(ggplot2)
library(RColorBrewer)

# Disease 缩写顺序
disease_levels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS","NF","SWN","LYM",
                    "HE","KHE","MPNST","US","AS","RMS","PECOMA","EWS","MRT")

# 对应自定义颜色
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

CNV_disease_durgs_sub$Disease <- factor(CNV_disease_durgs_sub$Disease, levels=disease_levels)

# 分组
CNV_disease_durgs_sub$plot_group <- with(CNV_disease_durgs_sub, 
                                         ifelse(logFC > 1 & av_exp > 1, "both",
                                                ifelse(logFC > 1 | av_exp > 1, "one", "none")))

# fill_factor: 高亮和单条件用 Disease 色，none 为灰色
CNV_disease_durgs_sub$fill_factor <- ifelse(CNV_disease_durgs_sub$plot_group=="none", "Other",
                                            as.character(CNV_disease_durgs_sub$Disease))
CNV_disease_durgs_sub$fill_factor <- factor(CNV_disease_durgs_sub$fill_factor, levels=c(disease_levels,"Other"))

# color_edge: 高亮点边框红色，其余按 fill_color
CNV_disease_durgs_sub$color_edge <- ifelse(CNV_disease_durgs_sub$plot_group=="both", "red",
                                           as.character(CNV_disease_durgs_sub$fill_factor))
CNV_disease_durgs_sub$color_edge <- factor(CNV_disease_durgs_sub$color_edge, levels=c(disease_levels,"Other","red"))

# 对应颜色向量
fill_colors <- c(setNames(cols, disease_levels), Other="grey")
edge_colors <- c(setNames(cols, disease_levels), Other="grey", red="red")

# 绘图
ggplot(CNV_disease_durgs_sub,
       aes(y=logFC, x=av_exp, size=prop, alpha=Freq,
           fill=fill_factor, color=color_edge)) +
  geom_point(shape=21, stroke=1.2) +   # stroke 控制边框粗细
  # 0 刻度线
  geom_hline(yintercept=0, color="black", size=0.7) +
  geom_vline(xintercept=0, color="black", size=0.7) +
  # 阈值红色虚线
  geom_vline(xintercept=1, color="red", linetype="dashed", size=0.8) +
  geom_hline(yintercept=1, color="red", linetype="dashed", size=0.8) +
  scale_size_continuous(range=c(2,10)) +
  scale_alpha_continuous(range=c(0.6,1)) +
  scale_fill_manual(values=fill_colors) +
  scale_color_manual(values=edge_colors) +
  theme_minimal() +
  theme(panel.grid=element_line(color="grey90")) +
  labs(y="logFC", x="Average expression",
       size="Prop", alpha="Freq (transparency)",
       fill="Disease", color="Edge color",
       title="Drug–Disease CNV visualization with red highlight edge")

# anno_1 <- drugs_anno[,c(2,4)] %>% as.data.frame()
# anno_1$cat <- rep("Drugs_anno",nrow(anno_1))
# anno_1 <- ggplot(data = anno_1, aes(x=cat,y=Drugs,fill=Anno_1))+
#   geom_tile() + theme_bw() + theme(panel.grid=element_blank())+
#   theme(axis.text.y = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks=element_blank()) +
#   ylab(NULL)
# anno_1
# 
# anno_2 <- drugs_anno[,c(2,6)] %>% as.data.frame()
# anno_2$cat <- rep("Drugs_anno",nrow(anno_2))
# anno_2 <- ggplot(data = anno_2, aes(x=cat,y=Drugs,fill=Top_cross_ref))+
#   geom_tile() + theme_bw() + theme(panel.grid=element_blank())+
#   theme(axis.text.y = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks=element_blank()) +
#   ylab(NULL)
# anno_2
saveRDS(CNV_disease_durgs_sub,"/cluster3/yflu/STS/Drug_screen/STS_CNV_disease_durgs_sub_251203.rds")

p1%>%
  insert_left(anno,width = 0.05)%>% 
  insert_left(p2,width = 0.2) 

i = 1
disease_samplenames <- rownames(subset(anno_sample_cluster_extended,Disease == diseasenames[i]))
CNV_disease_genes <- c()
disease_length <- length(disease_samplenames)
for (j in 1:disease_length) {
  CNV_disease_genes <- c(CNV_disease_genes,cnv_genes_top_list[disease_samplenames[j]][[1]]$CNV)
}
CNV_disease_genes <- as.data.frame(table(CNV_disease_genes))
CNV_disease_genes$prop <- CNV_disease_genes$Freq/disease_length
CNV_disease_genes$Disease <- rep(diseasenames[i],nrow(CNV_disease_genes))
colnames(CNV_disease_genes)[1] <- "CNVs"

for (i in 2:length(diseasenames)) {
  disease_samplenames <- rownames(subset(anno_sample_cluster_extended,Disease == diseasenames[i]))
  CNV_disease_genes_1 <- c()
  disease_length <- length(disease_samplenames)
  for (j in 1:disease_length) {
    CNV_disease_genes_1 <- c(CNV_disease_genes_1,cnv_genes_top_list[disease_samplenames[j]][[1]]$CNV)
  }
  CNV_disease_genes_1 <- as.data.frame(table(CNV_disease_genes_1))
  CNV_disease_genes_1$prop <- CNV_disease_genes_1$Freq/disease_length
  CNV_disease_genes_1$Disease <- rep(diseasenames[i],nrow(CNV_disease_genes_1))
  colnames(CNV_disease_genes_1)[1] <- "CNVs"
  CNV_disease_genes <- rbind(CNV_disease_genes,CNV_disease_genes_1)
}

CNV_disease_genes$CNVs <- as.character(CNV_disease_genes$CNVs)
CNV_disease_genes$status <- CNV_disease_genes$CNVs
CNV_disease_genes$genes <- CNV_disease_genes$CNVs

for (i in 1:nrow(CNV_disease_genes)) {
  CNV_disease_genes$status[i] <- strsplit(CNV_disease_genes$CNVs[i],"_")[[1]][2]
  CNV_disease_genes$genes[i] <- strsplit(CNV_disease_genes$CNVs[i],"_")[[1]][1]
  print(i)
}

saveRDS(CNV_disease_genes,"/cluster3/yflu/STS/Drug_screen/CNV_disease_genes_250915.rds")
saveRDS(CNV_disease_durgs,"/cluster3/yflu/STS/Drug_screen/CNV_disease_durgs_250915.rds")

CNV_disease_durgs_sub$geneset <- rownames(CNV_disease_durgs_sub$geneset)
for (i in 1:nrow(CNV_disease_durgs_sub)) {
  CNV_disease_durgs_sub$geneset[i] <- paste(drug_hgnc_list_combined[[as.character(CNV_disease_durgs_sub$Drugs[i])]], collapse = ", ")
}
write.xlsx(CNV_disease_durgs_sub,"CNV_based_drugs_250915.xlsx")


