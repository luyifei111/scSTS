library(Seurat)
library(scKidney)
library(tximport)  # 用于获取基因长度（推荐）
library(GenomicFeatures)
library(org.Hs.eg.db)
library(limma)
library(sva)
library(ggplot2)
library(ggrepel)
annot_file <- "/cluster3/yflu/STS/normal/Homo_sapiens.GRCh38.114.gtf"
# 1. 创建TxDb对象（从GTF/GFF文件）
txdb <- makeTxDbFromGFF(annot_file, format = "gtf")

# 2. 获取转录本长度信息
transcript_lengths <- transcriptLengths(txdb, 
                                        with.cds_len = FALSE,
                                        with.utr5_len = FALSE,
                                        with.utr3_len = FALSE)
counts_to_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)  # 每千碱基读数
  scaling_factor <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  return(tpm)
}

objlist <- read.csv("STS_sample_select_24.5.5.csv")
objlist <- subset(objlist,Disease1 %in% c("Benign","Malignant"))
seuobjnames <- objlist$Sample
seuobjpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",seuobjnames,"/",seuobjnames,".rds",sep = "")
outpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",seuobjnames,"/",seuobjnames,".h5ad",sep = "")

for (i in 1:length(seuobjpath)) {
  seuobj <- readRDS(seuobjpath[i])
  seuobj@meta.data <- seuobj@meta.data[colnames(seuobj@assays$RNA@counts),]
  kidneyH5(seuobj, outpath[i])
  print(paste(seuobjnames[i],"_",i,"/85",sep = ""))
}

celltypepath <- paste("/cluster3/yflu/STS/normal/celltype/",seuobjnames,"_","celltype.integrated.xlsx",sep = "")
celltypistpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",seuobjnames,"/",seuobjnames,"_celltypist.csv",sep = "")
celltype_merged_out <- paste("/cluster3/yflu/STS/normal/celltype_merged/",seuobjnames,"_","celltype.integrated_merged.xlsx",sep = "")
for (i in 80:length(seuobjnames)) {
  celltype <- openxlsx::read.xlsx(celltypepath[i],"Sheet1")
  celltypist <- read.csv(celltypistpath[i])
  celltype_merged <- cbind(celltype,celltypist)
  print(seuobjnames[i])
  openxlsx::write.xlsx(celltype_merged,celltype_merged_out[i])
}

#distance
IH_ref <- read.delim("/cluster3/yflu/STS/normal/PUBLIC/GSE225591_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rownames(IH_ref) <- IH_ref$GeneID

anno <- read.delim("/cluster3/yflu/STS/normal/PUBLIC/Human.GRCh38.p13.annot.tsv")
anno_1 <- subset(anno,GeneID %in% rownames(IH_ref))
rownames(anno_1) <- anno_1$GeneID
anno_1 <- anno_1[rownames(IH_ref),]
IH_ref$Symbol <- anno_1$Symbol
IH_ref <- IH_ref[!duplicated(IH_ref$Symbol),]
rownames(IH_ref) <- IH_ref$Symbol

IH_ref <- IH_ref[,-c(1,11:14)]
IH_ref_average <- as.data.frame(rowMeans(IH_ref))

normal.pega_20250624 <- readRDS("/cluster3/yflu/STS/normal/normal.pega_20250624.rds")
pdfpath <- paste("/cluster3/yflu/STS/normal/distance/",seuobjnames,"_distance.pdf")
normal.pega_endo <- subset(normal.pega_20250624,celltype %in% c("Endothelial cells"))
match("CH3",seuobjnames)

i=15
seuobj <- readRDS(seuobjpath[i])
seuobj@meta.data <- seuobj@meta.data[colnames(seuobj@assays$RNA@counts),]

celltype_merged <- openxlsx::read.xlsx(celltype_merged_out[i],"Sheet 1")
View(celltype_merged)
suspected <- c(6,10,13)
seuobj$seurat_clusters <- Idents(seuobj)
seuobj_1 <- subset(seuobj,seurat_clusters %in% suspected)
#seuobj_2 <- subset(seuobj,seurat_clusters %in% names(table(seuobj$seurat_clusters))[-(suspected+1)])
#Idents(seuobj_2) <- "all"
DimPlot(seuobj_1)

# seuobj_merged <- merge(seuobj_1,normal.pega_endo)
# seuobj_merged <- NormalizeData(seuobj_merged)
# seuobj_merged <- FindVariableFeatures(seuobj_merged, selection.method = "vst", nfeatures = 8000)
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# seuobj_merged <- CellCycleScoring(seuobj_merged, s.features = s.genes, g2m.features = g2m.genes)
# seuobj_merged <- ScaleData(seuobj_merged, vars.to.regress = c("S.Score", "G2M.Score"))
# seuobj_merged <- RunPCA(pc.genes = seuobj_merged@var.genes, npcs = 20 , object = seuobj_merged)
# ElbowPlot(seuobj_merged)
# seuobj_merged <- FindNeighbors(seuobj_merged, dims = 1:15)
# seuobj_merged <- FindClusters(seuobj_merged, resolution = 0.8)
# head(Idents(seuobj_merged), 5)
# seuobj_merged <- RunUMAP(seuobj_merged, dims = 1:15)
# DimPlot(seuobj_merged,group.by = "orig.ident")
# markers <- FindMarkers(seuobj_merged,ident.1 = "CH1",group.by = "orig.ident")

#markers <- FindMarkers(seuobj,ident.1 = 9,group.by = "seurat_clusters")

Idents(seuobj_1) <- seuobj_1$seurat_clusters

# 提取原始 counts
counts_mat <- GetAssayData(seuobj_1, slot = "counts")

# 提取细胞的 cluster 信息
cluster_ids <- Idents(seuobj_1)

# pseudobulk 计算：按 cluster 对 counts 矩阵列求和
pseudobulk_counts <- sapply(levels(cluster_ids), function(cluster) {
  cells_in_cluster <- names(cluster_ids[cluster_ids == cluster])
  Matrix::rowSums(counts_mat[, cells_in_cluster, drop = FALSE])
})

# 转为 data.frame（行名是基因，列名是 cluster）
pseudobulk_counts <- as.data.frame(pseudobulk_counts)

#seuobj_2_pseudobulk_counts <- as.data.frame(rowSums(GetAssayData(seuobj_2, slot = "counts")))
normal.pega_endo_pseudobulk_counts <- as.data.frame(rowSums(GetAssayData(normal.pega_endo, slot = "counts")))
genes <- intersect(rownames(pseudobulk_counts),rownames(normal.pega_endo_pseudobulk_counts))

# 3. 按基因ID聚合（取最长转录本的长度）
gene_lengths <- aggregate(tx_len ~ gene_id, 
                          data = transcript_lengths, 
                          FUN = max)
id_result <- select(
  org.Hs.eg.db,
  keys = gene_lengths$gene_id,
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID"),
  keytype = "ENSEMBL"
)

id_result <- id_result[!duplicated(id_result$ENSEMBL),]

gene_lengths$symbol <- id_result$SYMBOL
common_genes <- intersect(genes, gene_lengths$symbol)
gene_lengths <- gene_lengths[match(common_genes, gene_lengths$symbol), ]

pseudobulk_counts_1 <- as.data.frame(pseudobulk_counts[common_genes, ])
#seuobj_2_pseudobulk_counts_1 <- as.data.frame(seuobj_2_pseudobulk_counts[common_genes, ])
normal.pega_endo_pseudobulk_counts_1 <- as.data.frame(normal.pega_endo_pseudobulk_counts[common_genes, ])
rownames(pseudobulk_counts_1) <- common_genes
#rownames(seuobj_2_pseudobulk_counts_1) <- common_genes
rownames(normal.pega_endo_pseudobulk_counts_1) <- common_genes

pseudobulk_counts_1_tpm <- counts_to_tpm(pseudobulk_counts_1, gene_lengths$tx_len)
#seuobj_2_tpm <- counts_to_tpm(seuobj_2_pseudobulk_counts_1, gene_lengths$tx_len)
normal.pega_endo_tpm <- counts_to_tpm(normal.pega_endo_pseudobulk_counts_1, gene_lengths$tx_len)

genes_intersect <- intersect(rownames(IH_ref_average),rownames(pseudobulk_counts_1_tpm))

average_merged <- cbind(IH_ref_average[genes_intersect,],normal.pega_endo_tpm[genes_intersect,])
#average_merged <- cbind(average_merged,seuobj_2_tpm[genes_intersect,])
average_merged <- cbind(average_merged,pseudobulk_counts_1_tpm[genes_intersect,])

# 批次信息 (例如: batch <- c(1,1,1,2,2,2))
batch <- factor(c(rep("Batch1", 1), rep("Batch2", 1+length(suspected)))) # 替换为你的批次信息

# 对TPM数据进行log2转换 (加1伪计数避免log(0))
log_tpm <- log2(average_merged + 1)

# 使用ComBat进行批次校正
combat_tpm <- ComBat(
  dat = as.matrix(log_tpm),
  batch = batch,  # 保留生物学差异
  par.prior = TRUE,  # 使用参数先验
  prior.plots = FALSE
)

# 转换回线性空间 (可选)
corrected_tpm <- 2^combat_tpm - 1
rownames(corrected_tpm) <- genes_intersect

corrected_tpm_DF <- as.data.frame(corrected_tpm)
corrected_tpm_DF$gene <- rownames(corrected_tpm_DF)
# #markers_1 <- rownames(markers)[1:500]
# genes_intersect_1 <- intersect(genes_intersect,markers_1)
# 
# corrected_tpm_1 <- corrected_tpm[genes_intersect_1,]

distance_matrix <- dist(t(combat_tpm),method = "euclidean")
distance_matrix <- as.matrix(distance_matrix)
distance_matrix
saveRDS(distance_matrix,paste("/cluster3/yflu/STS/normal/distance/",seuobjnames[i],"_distance_matrix.rds",sep = ""))

distance_ref <- distance_matrix[1,3]
distance_normal <- distance_matrix[2,3]
distance_merged <- as.data.frame(cbind(distance_ref,distance_normal))

if(length(suspected) > 1){
  for (j in 2:length(suspected)) {
    distance_ref <- distance_matrix[1,j+2]
    distance_normal <- distance_matrix[2,j+2]
    distance_merged_1 <- as.data.frame(cbind(distance_ref,distance_normal))
    distance_merged <- rbind(distance_merged,distance_merged_1)
  }
}
rownames(distance_merged) <- paste("c",suspected,sep = "")

# Calculate T position on a normalized scale (0 = T_ref, 1 = T_normal)
T_positions <- distance_merged$distance_ref / (distance_merged$distance_ref + distance_merged$distance_normal)

# Create data frame
points <- data.frame(
  Label = c("ref", "normal",rownames(distance_merged)),
  Position = c(0, 1,T_positions)
)
midpoint <- 0.5

labels <- points$Label
colors <- setNames(
  c("red", "blue", rep("purple", length(labels) - 2)),
  labels
)
points$ColorLabel <- ifelse(points$Position >= 0.45 & points$Position <= 0.55,
                            "gray", points$Label)
points$TextColor <- ifelse(points$Position >= 0.45 & points$Position <= 0.55,
                           "gray40", "purple")
color_values <- colors
color_values["gray"] <- "gray40"
# Horizontal plot

p <- ggplot(points, aes(x = Position, y = 0))+
  
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "gray60", size = 1) +
  geom_point(aes(color = ColorLabel), size = 4) +
  
  geom_text_repel(
    data = points[!points$Label %in% c("ref", "normal"), ],
    aes(label = Label, color = NULL),  # 取消与 aes(color = ...) 冲突
    size = 5,
    nudge_y = 0.01,
    max.overlaps = Inf,
    color = points$TextColor[!points$Label %in% c("ref", "normal")]
  ) +
  
  # 数值
  geom_text_repel(
    data = points[!points$Label %in% c("ref", "normal"), ],
    aes(label = sprintf("%.2f", Position), color = NULL),
    size = 4,
    nudge_y = -0.01,
    max.overlaps = Inf,
    color = points$TextColor[!points$Label %in% c("ref", "normal")]
  )+
  
  annotate("text", x = 0 - 0.03, y = 0.03, label = "ref", size = 5, color = "red", hjust = 1) +
  annotate("text", x = 1 + 0.03, y = 0.03, label = "normal", size = 5, color = "blue", hjust = 0) +
  
  geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.01), color = "red", size = 0.5) +
  geom_segment(aes(x = 1, xend = 1, y = -0.01, yend = 0.01), color = "blue", size = 0.5) +
  geom_segment(aes(x = midpoint, xend = midpoint, y = -0.005, yend = 0.005), color = "darkgray", linetype = "dashed") +
  annotate("text", x = midpoint, y = -0.025, label = "Midpoint", size = 4, color = "darkgray") +
  
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Position of suspected clusters Between ref and normal") +
  scale_color_manual(values = colors) + 
  geom_rect(aes(xmin = 0.45, xmax = 0.55, ymin = -0.01, ymax = 0.01),
                                                fill = "gray80", alpha = 0.03, inherit.aes = FALSE)
p
ggsave(pdfpath[i], plot = p, width = 8, height = 3, units = "in")
