library(data.table)
library(ggplot2)
library(RColorBrewer)
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
library(ggcorrplot)
library(tidyr)
library(maftools)
library(pbapply)
library(S4Vectors)
library(pheatmap)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
samplenames <- rownames(anno_sample_cluster_extended)
cnvscorepath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.xlsx",sep = "")

i=1
cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
cnv_score$group2 <- factor(cnv_score$group2,levels = c("Tumor","Normal"))
cnv_stat <- cnv_score %>%
  summarise(
    mean_group1 = mean(cnvscore[group2 == "Tumor"], na.rm = TRUE),
    mean_group2 = mean(cnvscore[group2 == "Normal"], na.rm = TRUE),
    p_value = t.test(cnvscore ~ group2)$p.value
  )
cnv_stat$mean <- cnv_stat[1,1] - cnv_stat[1,2]
for (i in 2:length(samplenames)) {
  cnv_score <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
  cnv_stat_1 <- cnv_score %>%
    summarise(
      mean_group1 = mean(cnvscore[group2 == "Tumor"], na.rm = TRUE),
      mean_group2 = mean(cnvscore[group2 == "Normal"], na.rm = TRUE),
      p_value = t.test(cnvscore ~ group2)$p.value
    )
  cnv_stat_1$mean <- cnv_stat_1[1,1] - cnv_stat_1[1,2]
  cnv_stat <- rbind(cnv_stat,cnv_stat_1)
  print(i)
}
cnv_stat <- as.data.frame(cnv_stat)
rownames(cnv_stat) <- samplenames

cnv_score_sample <- cbind(cnv_stat,anno_sample_cluster_extended)
cnv_score_sample <- cnv_score_sample %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))
cnv_score_sample$Disease <- factor(cnv_score_sample$Disease,levels = labels)

cnv_score_sample <- cnv_score_sample %>%
  arrange(Disease, desc(mean))

mat_mean <- matrix(
  cnv_score_sample$mean,
  ncol = 1
)
rownames(mat_mean) <- rownames(cnv_score_sample)
colnames(mat_mean) <- "CNV mean"
mat_sig <- matrix(
  cnv_score_sample$sig,
  ncol = 1
)

anno_row <- cnv_score_sample %>%
  dplyr::select(Disease)
rownames(anno_row) <- rownames(cnv_score_sample)

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

disease_levels <-labels
# 如果颜色数 ≥ Disease 数（推荐）
names(disease_palette) <- disease_levels

annotation_colors <- list(
  Disease = disease_palette
)

min_val <- min(mat_mean, na.rm = TRUE)
max_val <- max(mat_mean, na.rm = TRUE)
n_col <- 100
n_half <- n_col / 2

breaks <- c(
  seq(min_val, 0, length.out = n_half + 1),
  seq(0, max_val, length.out = n_half + 1)[-1]
)
cols <- colorRampPalette(
  c("#2166AC", "#FFFFBF", "#B2182B")
)(n_col)

pheatmap(
  mat_mean,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cols,
  breaks = breaks,
  annotation_row = anno_row,
  annotation_colors = annotation_colors,
  display_numbers = mat_sig,
  number_color = "black",
  fontsize_number = 10
)

markers_sub <- openxlsx::read.xlsx("/cluster3/yflu/STS/pegasus/markers_membrane_new.xlsx","Sheet 1")
aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

diseasenames <- labels
diseasenames <- gsub("[ ]",".",  diseasenames)
diseasenames[19] <- "EWS_PNET"
i = 1
index <- match(diseasenames[i],colnames(markers_sub))
markers_sub_index <- markers_sub[,c(index,index+1,index+2)]
colnames(markers_sub_index)[2] <- "cdq"
markers_sub_index[,3] <- factor(markers_sub_index[,3],levels = c(1:100))
markers_sub_index <- markers_sub_index[order(markers_sub_index[,3]),]
#markers_sub_index <- subset(markers_sub_index,cdq == "Membrane")
markers_sub_index <- markers_sub_index[c(1:3),1]
for (i in 2:length(diseasenames)) {
  index <- match(diseasenames[i],colnames(markers_sub))
  markers_sub_index_1 <- markers_sub[,c(index,index+1,index+2)]
  colnames(markers_sub_index_1)[2] <- "cdq"
  markers_sub_index_1[,3] <- factor(markers_sub_index_1[,3],levels = c(1:100))
  markers_sub_index_1 <- markers_sub_index_1[order(markers_sub_index_1[,3]),]
 #markers_sub_index_1 <- subset(markers_sub_index_1,cdq == "Membrane")
  max <- nrow(markers_sub_index_1)
  markers_sub_index_2 <- markers_sub_index_1[c(1:3),1]
  markers_sub_index <- c(markers_sub_index,markers_sub_index_2)
  markers_sub_index <- unique(markers_sub_index)
  for (j in 1:(max-3)) {
    if(length(markers_sub_index) < i*3){
      markers_sub_index <- c(markers_sub_index,markers_sub_index_1[c(3+j),1])
      markers_sub_index <- unique(markers_sub_index)
    }
  }
}
markers_sub_index[46] <- "IGF2"

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5ad")
STS.pega.tumor <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.tumor@meta.data <- metadata$obs

STS.pega.tumor$barcode <- rownames(STS.pega.tumor@meta.data)
meta_sampled <- STS.pega.tumor@meta.data %>% group_by(louvain_labels) %>% slice_sample(prop = .05)
STS.pega.tumor_sampled <- subset(STS.pega.tumor,barcode %in% meta_sampled$barcode)

STS.pega.tumor_sampled$Disease_ordered <- factor(STS.pega.tumor_sampled$Disease,levels = labels)
DotPlot(STS.pega.tumor_sampled,features = markers_sub_index,group.by = "Disease_ordered")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))

P = DotPlot(STS.pega.tumor_sampled,features = markers_sub_index,group.by = "Disease_ordered",cols = c("#4575B4", "#D73027"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
P
data <- P$data
order_id <- rev(levels(data$id))
data$id <- as.character(data$id)
data$id <- factor(data$id,levels = order_id)
levels(data$features.plot) <- rev(levels(data$features.plot))

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])
diseasenames <- unique(as.character(STS.pega.tumor$Disease))
diseasenames <- as.data.frame(cbind(diseasenames,diseasenames))
colnames(diseasenames) <- c("Var1","Var2")
rownames(diseasenames) <- paste("STS|",diseasenames$Var1,sep="")
diseasenames <- diseasenames[rownames(aurocs_disease)[p$tree_row$order],]
names(cols) <- diseasenames$Var1
data$cols <- as.character(data$id)
for (i in 1:1200) {
  data$cols[i] <- cols[data$id[i]]
}
ggplot(data,aes(x=features.plot,y=id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_tile(data=data, aes(x=-0.5, y=id, fill=id), width=1, height=1)+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))+
  labs(x=NULL,y=NULL)

benign_order <- subset(cnv_score_sample,Malignancy == "Benign")
distance_file <- paste("/cluster3/yflu/STS/normal/distance/",rownames(benign_order),"_distance_matrix.rds",sep = "")

i=1

distance_matrix <- readRDS(distance_file[i])
distance_ref <- distance_matrix[1,3]
distance_normal <- distance_matrix[2,3]
distance_merged <- as.data.frame(cbind(distance_ref,distance_normal))
if(nrow(distance_matrix) > 3){
  for (j in 2:(nrow(distance_matrix)-2)) {
    distance_ref <- distance_matrix[1,j+2]
    distance_normal <- distance_matrix[2,j+2]
    distance_merged_1 <- as.data.frame(cbind(distance_ref,distance_normal))
    distance_merged <- rbind(distance_merged,distance_merged_1)
  }
}
rownames(distance_merged) <- paste("c",rownames(as.data.frame(distance_matrix))[3:nrow(distance_matrix)],sep = "")

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
  
  geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.01), color = "red", size = 0.5) +
  geom_segment(aes(x = 1, xend = 1, y = -0.01, yend = 0.01), color = "blue", size = 0.5) +
  geom_segment(aes(x = midpoint, xend = midpoint, y = -0.005, yend = 0.005), color = "darkgray", linetype = "dashed") +

  theme_minimal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = rownames(benign_order)[i]) +
  scale_color_manual(values = colors) + 
  geom_rect(aes(xmin = 0.45, xmax = 0.55, ymin = -0.01, ymax = 0.01),
            fill = "gray80", alpha = 0.03, inherit.aes = FALSE)
p
plot_list <- list()
plot_list[[1]] <- p

for (i in 2:length(rownames(benign_order))) {
  distance_matrix <- readRDS(distance_file[i])
  distance_ref <- distance_matrix[1,3]
  distance_normal <- distance_matrix[2,3]
  distance_merged <- as.data.frame(cbind(distance_ref,distance_normal))
  if(nrow(distance_matrix) > 3){
    for (j in 2:(nrow(distance_matrix)-2)) {
      distance_ref <- distance_matrix[1,j+2]
      distance_normal <- distance_matrix[2,j+2]
      distance_merged_1 <- as.data.frame(cbind(distance_ref,distance_normal))
      distance_merged <- rbind(distance_merged,distance_merged_1)
    }
  }
  rownames(distance_merged) <- paste("c",rownames(as.data.frame(distance_matrix))[3:nrow(distance_matrix)],sep = "")
  
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
  
  p1 <- ggplot(points, aes(x = Position, y = 0))+
    
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "gray60", size = 1) +
    geom_point(aes(color = ColorLabel), size = 4) +
    
    geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.01), color = "red", size = 0.5) +
    geom_segment(aes(x = 1, xend = 1, y = -0.01, yend = 0.01), color = "blue", size = 0.5) +
    geom_segment(aes(x = midpoint, xend = midpoint, y = -0.005, yend = 0.005), color = "darkgray", linetype = "dashed") +

    theme_minimal() +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = rownames(benign_order)[i]) +
    scale_color_manual(values = colors) + 
    geom_rect(aes(xmin = 0.45, xmax = 0.55, ymin = -0.01, ymax = 0.01),
              fill = "gray80", alpha = 0.03, inherit.aes = FALSE)
  plot_list[[i]] <- p1
}

final_plot <- wrap_plots(plot_list, ncol = 1)

ggsave(
  filename = "distance_position_vertical.pdf",
  plot = final_plot,
  limitsize = FALSE,
  width = 12,
  height = 2 * length(plot_list)  # 每张图 3 英寸
)

metadata <- read_h5ad("adata_TMA_merged_tumor.h5ad")
# Convert("adata_TMA_merged_tumor.h5ad", dest = "h5seurat", overwrite = T)
# f <- H5File$new("adata_TMA_merged_tumor.h5seurat", "r+")
# groups <- f$ls(recursive = TRUE)
# 
# for (name in groups$name[grepl("categories", groups$name)]) {
#   names <- strsplit(name, "/")[[1]]
#   names <- c(names[1:length(names) - 1], "levels")
#   new_name <- paste(names, collapse = "/")
#   f[[new_name]] <- f[[name]]
# }
# 
# for (name in groups$name[grepl("codes", groups$name)]) {
#   names <- strsplit(name, "/")[[1]]
#   names <- c(names[1:length(names) - 1], "values")
#   new_name <- paste(names, collapse = "/")
#   f[[new_name]] <- f[[name]]
#   grp <- f[[new_name]]
#   grp$write(args = list(1:grp$dims), value = grp$read() + 1)
# }
# 
# f$close_all()

metadata <- read_h5ad("/cluster3/yflu/STS/TMA/adata_TMA_merged_tumor.h5ad")

TMA_merged_tumor <- LoadH5Seurat("/cluster3/yflu/STS/TMA/adata_TMA_merged_tumor.h5seurat",meta.data = FALSE, misc = FALSE)
TMA_merged_tumor@meta.data <- metadata$obs

panel_5k <- read.csv("/cluster3/yflu/STS/TMA/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")

markers_sub <- openxlsx::read.xlsx("/cluster3/yflu/STS/pegasus/markers_membrane_new.xlsx","Sheet 1")

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

diseasenames <- labels
diseasenames <- gsub("[ ]",".",  diseasenames)
diseasenames[19] <- "EWS_PNET"
i = 1
index <- match(diseasenames[i],colnames(markers_sub))
markers_sub_index <- markers_sub[,c(index,index+1,index+2)]
colnames(markers_sub_index)[2] <- "cdq"
#markers_sub_index <- subset(markers_sub_index,cdq == "Membrane")
markers_sub_index$gene <- markers_sub_index[,1]
markers_sub_index <- subset(markers_sub_index,gene %in% panel_5k$gene_name)
markers_sub_index[,3] <- factor(markers_sub_index[,3],levels = c(1:100))
markers_sub_index <- markers_sub_index[order(markers_sub_index[,3]),]
markers_sub_index <- markers_sub_index[c(1:2),1]

for (i in 2:length(diseasenames)) {
  index <- match(diseasenames[i],colnames(markers_sub))
  markers_sub_index_1 <- markers_sub[,c(index,index+1,index+2)]
  colnames(markers_sub_index_1)[2] <- "cdq"
  #markers_sub_index_1 <- subset(markers_sub_index_1,cdq == "Membrane")
  markers_sub_index_1$gene <- markers_sub_index_1[,1]
  markers_sub_index_1 <- subset(markers_sub_index_1,gene %in% panel_5k$gene_name)
  max <- nrow(markers_sub_index_1)
  markers_sub_index_1[,3] <- factor(markers_sub_index_1[,3],levels = c(1:100))
  markers_sub_index_1 <- markers_sub_index_1[order(markers_sub_index_1[,3]),]
  markers_sub_index_2 <- markers_sub_index_1[c(1:2),1]
  markers_sub_index <- c(markers_sub_index,markers_sub_index_2)
  markers_sub_index <- unique(markers_sub_index)
  for (j in 1:(max-2)) {
    if(length(markers_sub_index) < i*2){
      markers_sub_index <- c(markers_sub_index,markers_sub_index_1[c(2+j),1])
      markers_sub_index <- unique(markers_sub_index)
    }
  }
}
labels <- c("LPB","SS","LPS","SCT","IFS","AF","IMT","ASPS","NF","SWN","LYM","HE","KHE",
            "MPNST","US","AS","RMS","PECOMA","EWS","MRT")
TMA_merged_tumor$Disease_ordered <- factor(TMA_merged_tumor$Disease,levels = labels)
DotPlot(TMA_merged_tumor,features = markers_sub_index,group.by = "Disease_ordered")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
TMA_merged_tumor$Disease_ordered <- factor(TMA_merged_tumor$Disease,levels = rev(labels))
#markerdata <- ScaleData(TMA_merged_tumor, features = markers_sub_index, assay = "RNA")
P = DotPlot(TMA_merged_tumor,features = rev(markers_sub_index),group.by = "Disease_ordered",cols = c("#4575B4", "#D73027"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
data <- P$data
order_id <- levels(data$id)
data$id <- as.character(data$id)
data$id <- factor(data$id,levels = order_id)
data$features.plot <- factor(data$features.plot,levels = rev(levels(data$features.plot)))
P
cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])
diseasenames <- unique(as.character(TMA_merged_tumor$Disease))
diseasenames <- as.data.frame(cbind(diseasenames,diseasenames))
colnames(diseasenames) <- c("Var1","Var2")
#rownames(diseasenames) <- paste("STS|",diseasenames$Var1,sep="")
rownames(diseasenames) <- diseasenames[,1]
diseasenames <- diseasenames[labels,]
names(cols) <- diseasenames$Var1
data$cols <- as.character(data$id)
for (i in 1:400) {
  data$cols[i] <- cols[data$id[i]]
}
ggplot(data,aes(x=features.plot,y=id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_tile(data=data, aes(x=-0.5, y=id, fill=id), width=1, height=1)+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))+
  labs(x=NULL,y=NULL)

library(RColorBrewer)
library(EnhancedVolcano)
library(scales)
library(DoubletFinder)
library(caret)
library(hdf5r)
library(ggpubr)
library(irGSEA)

metadata <- read_h5ad("/cluster3/yflu/STS/TMA/TMA_normal_pegasus_250823.h5ad")
Convert("/cluster3/yflu/STS/TMA/TMA_normal_pegasus_250823.h5ad", dest = "h5seurat", overwrite = T)
f <- H5File$new("/cluster3/yflu/STS/TMA/TMA_normal_pegasus_250823.h5seurat", "r+")
groups <- f$ls(recursive = TRUE)

for (name in groups$name[grepl("categories", groups$name)]) {
  names <- strsplit(name, "/")[[1]]
  names <- c(names[1:length(names) - 1], "levels")
  new_name <- paste(names, collapse = "/")
  f[[new_name]] <- f[[name]]
}

for (name in groups$name[grepl("codes", groups$name)]) {
  names <- strsplit(name, "/")[[1]]
  names <- c(names[1:length(names) - 1], "values")
  new_name <- paste(names, collapse = "/")
  f[[new_name]] <- f[[name]]
  grp <- f[[new_name]]
  grp$write(args = list(1:grp$dims), value = grp$read() + 1)
}
f$close_all()

TMA_merged_normal <- LoadH5Seurat("/cluster3/yflu/STS/TMA/TMA_normal_pegasus_250823.h5seurat",meta.data = FALSE, misc = FALSE)
TMA_merged_normal@meta.data <- metadata$obs

TMA_merged_normal$celltype <- factor(TMA_merged_normal$celltype,levels = unique(TMA_merged_normal$celltype)[c(1,9,11,3,14,6,7,2,4,8,10,5,12,13)])

DimPlot(TMA_merged_normal,reduction = "tsne",raster = F,group.by = "celltype",
        cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(6)[c(1:2)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(10)[c(1:3)],
                 colorRampPalette(brewer.pal(12,'Set3')[c(9,10)])(4)[c(1:3)]))

markers <- c("COL5A1","MMP14","SERPINH1",
             "MYH2","ACTN2","PDGFA",
             "RGS5","NOTCH3","COL5A2",
             "PECAM1","PLVAP","CD34",
             "CXCL12","SLC40A1","NR1H4",
             "CD79A","CD19","MS4A1",
             "XBP1","PIM2","DERL3",
             "CD3E","CD8A","GZMA",
             "CD68","CD14","CD163",
             "CLEC9A","LILRA4","CLEC4C",
             "KIT","CTSG","IL1RL1",
             "SOX1","SIM2","GRP",
             "ADIPOQ","PLIN1","LPL",
             "CFB","SERPINA3","PLG")

markers_average <- AverageExpression(TMA_merged_normal,features = markers,group.by = "celltype")
markers_average <- markers_average$RNA
markers_average <- t(markers_average)
markers_average <- scale(markers_average)
pheatmap::pheatmap(markers_average,cluster_rows = F,cluster_cols = F)
i=1
anno <- as.data.frame(rep(rownames(markers_average)[i],3))
colnames(anno) <- "Markers"
rownames(anno) <- colnames(markers_average)[(i*3-2):(3*i)]
for (i in 2:nrow(markers_average)) {
  anno_1 <- as.data.frame(rep(rownames(markers_average)[i],3))
  colnames(anno_1) <- "Markers"
  rownames(anno_1) <- colnames(markers_average)[(i*3-2):(3*i)]
  anno <- rbind(anno,anno_1)
  print(i)
}
ann_colors = c(colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(10)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(6)[c(1:2)],
               colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(10)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(9,10)])(4)[c(1:3)])
names(ann_colors) <- levels(TMA_merged_normal$celltype)
ann_colors <- list(Markers = ann_colors)

pheatmap::pheatmap(markers_average,cluster_rows = F,cluster_cols = F,annotation_col = anno,annotation_colors = ann_colors)