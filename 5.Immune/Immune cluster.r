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

metadata <- read_h5ad("/cluster3/yflu/STS/microenvironmnet/data_normal_B.h5ad")
STS_pega_B <- LoadH5Seurat("/cluster3/yflu/STS/microenvironmnet/data_normal_B.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega_B@meta.data <- metadata$obs
Idents(STS_pega_B) <- STS_pega_B@meta.data$louvain_labels_2
celltype <- read.xlsx("/cluster3/yflu/STS/microenvironmnet/B/celltype.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS_pega_B)
STS_pega_B <- RenameIdents(STS_pega_B, new.cluster.ids)
STS_pega_B@meta.data$Celltype_new <- Idents(STS_pega_B)
celltypes <- unique(levels(STS_pega_B@meta.data$Celltype_new))[-c(5,10)]
STS_pega_B <- subset(STS_pega_B,Celltype_new %in% celltypes)
STS_pega_B@meta.data$Celltype_new <- droplevels(STS_pega_B@meta.data$Celltype_new)
DimPlot(STS_pega_B, reduction = "umap",label = TRUE,raster = F)

metadata <- read_h5ad("/cluster3/yflu/STS/microenvironmnet/data_normal_TNK.h5ad")
STS_pega_TNK <- LoadH5Seurat("/cluster3/yflu/STS/microenvironmnet/data_normal_TNK.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega_TNK@meta.data <- metadata$obs
Idents(STS_pega_TNK) <- STS_pega_TNK@meta.data$louvain_labels_4
celltype <- read.xlsx("/cluster3/yflu/STS/microenvironmnet/T/celltypes res4.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS_pega_TNK)
STS_pega_TNK <- RenameIdents(STS_pega_TNK, new.cluster.ids)
STS_pega_TNK@meta.data$Celltype_new <- Idents(STS_pega_TNK)
celltypes <- unique(levels(STS_pega_TNK@meta.data$Celltype_new))[-15]
STS_pega_TNK <- subset(STS_pega_TNK,Celltype_new %in% celltypes)
STS_pega_TNK@meta.data$Celltype_new <- droplevels(STS_pega_TNK@meta.data$Celltype_new)
DimPlot(STS_pega_TNK, reduction = "umap",label = TRUE,raster = F)
FeaturePlot(STS_pega_TNK,feature = c("KIT"))
metadata <- read_h5ad("/cluster3/yflu/STS/microenvironmnet/data_normal_M.h5ad")
STS_pega_M <- LoadH5Seurat("/cluster3/yflu/STS/microenvironmnet/data_normal_M.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega_M@meta.data <- metadata$obs
Idents(STS_pega_M) <- STS_pega_M@meta.data$louvain_labels_3
celltype <- read.xlsx("/cluster3/yflu/STS/microenvironmnet/M/celltypes res3.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS_pega_M)
STS_pega_M <- RenameIdents(STS_pega_M, new.cluster.ids)
STS_pega_M@meta.data$Celltype_new <- Idents(STS_pega_M)
celltypes <- unique(levels(STS_pega_M@meta.data$Celltype_new))[-9]
STS_pega_M <- subset(STS_pega_M,Celltype_new %in% celltypes)
STS_pega_M@meta.data$Celltype_new <- droplevels(STS_pega_M@meta.data$Celltype_new)
DimPlot(STS_pega_M, reduction = "umap",label = TRUE,raster = F)

STS_pega_immune <- merge(STS_pega_B,c(STS_pega_TNK,STS_pega_M))
average_score_sample <- readRDS("/cluster3/yflu/STS/nmf/average_score_sample.rds")
average_score_sample <- as.data.frame(average_score_sample)

cellnumber <- read.csv("/cluster/home/yflu/STS/pegasus/cellnumber_240722.csv")
rownames(cellnumber) <- cellnumber$X
tumorcells <- read.csv("/cluster/home/yflu/STS/pegasus/cellnumber_tumor_240505.csv")
tumorcells <- subset(tumorcells,tumorcells$Channel > 0)
tumornames <- tumorcells$X
cellnumber_tumor <- cellnumber[tumornames,]
tumorcell_proportion <- tumorcells
tumorcell_proportion$proportion <- tumorcell_proportion$Channel
for (i in 1:79) {
  tumorcell_proportion$proportion[i] <- tumorcell_proportion$Channel[i]/cellnumber_tumor$Channel[i]
}

#programgenes <- read.csv("/cluster/home/yflu/STS/pegasus/nmf/allprogramgene_top50.csv")

rownames(tumorcells) <- tumorcells$X
anno <- read.csv("/cluster/home/yflu/STS/pegasus/STS_sample_select_24.5.5.csv")

#table_T <- as.data.frame.array(table(STS_pega_TNK@meta.data$Celltype_new,STS_pega_TNK@meta.data$Channel))
#table_T <- table_T[,tumornames]
#table_M <- as.data.frame.array(table(STS_pega_M@meta.data$Celltype_new,STS_pega_M@meta.data$Channel))
#table_M <- table_M[,tumornames]
#table_B <- as.data.frame.array(table(STS_pega_B@meta.data$Celltype_new,STS_pega_B@meta.data$Channel))
#table_B <- table_B[,tumornames]

#table_merged <- rbind(table_T,table_M)
#table_merged <- rbind(table_merged,table_B)
table_merged <- as.data.frame.array(table(STS_pega_immune@meta.data$Celltype_new,STS_pega_immune@meta.data$Channel))
table_merged <- t(table_merged)
table_merged <- table_merged[tumornames,]
#table_merged <- scale(table_merged)
table_merged <- na.omit(t(table_merged))
table_merged <- as.matrix(table_merged)

pheatmap(table_merged)
table_merged_proportion <- table_merged
for (i in 1:length(colnames(table_merged))) {
  table_merged_proportion[,i] <- table_merged_proportion[,i]/cellnumber[tumornames[i],2]
}
table_merged_proportion <- t(table_merged_proportion)
tumorcell_proportion <- tumorcell_proportion[colnames(average_score_sample),]
table_merged_proportion <- table_merged_proportion[colnames(average_score_sample),]
table_merged_proportion <- cbind(table_merged_proportion,tumorcell_proportion$proportion)
colnames(table_merged_proportion)[40] <- "Tumor"
rownames(average_score_sample) <- substr(rownames(average_score_sample),1,8)
table_merged_proportion <- cbind(table_merged_proportion,t(average_score_sample))
saveRDS(table_merged_proportion,"table_merged_proportion_raw.rds")

#restart here
table_merged_proportion <- readRDS("/cluster3/yflu/STS/cpdb/table_merged_proportion_raw.rds")
table_merged_proportion <- scale(table_merged_proportion)
table_merged_proportion <- na.omit(t(table_merged_proportion))
table_merged_proportion <- as.matrix(table_merged_proportion)

colors = colorRampPalette(brewer.pal(8,'RdYlBu'))(90)
colors1 = colorRampPalette(brewer.pal(8,'RdYlBu'))(155)
color = c(colors1[1:85],colors[75:90])
pheatmap(table_merged_proportion,color = rev(color))

rownames(anno) <- anno$Sample
anno <- anno[tumornames,]
anno_sample <- cbind(anno$Disease,anno$Disease1)
rownames(anno_sample) <- rownames(anno)
colnames(anno_sample) <- c("Disease","Malignancy")
anno_sample <- as.data.frame(anno_sample)

anno_immune <- as.character(unique(STS_pega_B@meta.data$Celltype_new))
anno_immune <- anno_immune[order(anno_immune)]
anno_immune <- as.data.frame(cbind(anno_immune,rep("B cells",length(anno_immune))))
colnames(anno_immune) <- c("celltype","Group")
rownames(anno_immune) <- anno_immune$celltype

anno_immune_1 <- as.character(unique(STS_pega_TNK@meta.data$Celltype_new))
anno_immune_1 <- anno_immune_1[order(anno_immune_1)]
anno_immune_1 <- as.data.frame(cbind(anno_immune_1,rep("T/NK cells",length(anno_immune_1))))
colnames(anno_immune_1) <- c("celltype","Group")
rownames(anno_immune_1) <- anno_immune_1$celltype
anno_immune <- rbind(anno_immune,anno_immune_1)

anno_immune_1 <- as.character(unique(STS_pega_M@meta.data$Celltype_new))
anno_immune_1 <- anno_immune_1[order(anno_immune_1)]
anno_immune_1 <- as.data.frame(cbind(anno_immune_1,rep("Myeloid cells",length(anno_immune_1))))
colnames(anno_immune_1) <- c("celltype","Group")
rownames(anno_immune_1) <- anno_immune_1$celltype
anno_immune <- rbind(anno_immune,anno_immune_1)

anno_immune <- subset(anno_immune,celltype %in% rownames(table_merged_proportion))

#jump here
anno_immune <- readRDS("/cluster3/yflu/STS/microenvironmnet/anno_immune.rds")

anno_immune_1 <- anno_immune[,-1]
anno_immune_1 <- as.data.frame(anno_immune_1)
rownames(anno_immune_1) <- rownames(anno_immune)
colnames(anno_immune_1) <- "Group"

anno_immune_nmf <- as.data.frame(rep("NMF program",8))
colnames(anno_immune_nmf) <- "Group"
rownames(anno_immune_nmf) <- rownames(table_merged_proportion)[c(40:47)]
anno_immune_1 <- rbind(anno_immune_1,anno_immune_nmf)
anno_immune_1 <- rbind(anno_immune_1,"Tumor cells")
rownames(anno_immune_1)[47] <- "Tumor"
anno_immune_1$Group_1 <- c("Bgc","Bm","Bpro","Bn","Bn","Plasma","Bn","Plasma","Bn",
                         "Tm","NK","Tem","Tem","Tex","Tfh","T ISG","NK","ILC3","Tc17","Temra","Tn","Treg","Tem","gdT",
                         "Mono/Macro","cDC","Mono/Macro","cDC","Mono/Macro","Mono/Macro","Mono/Macro","Mono/Macro","Mono/Macro","cDC","Mast",
                         "pDC","Mono/Macro","Mono/Macro",rownames(anno_immune_1)[39:46],"Tumor")

colors = colorRampPalette(brewer.pal(8,'RdYlBu'))(90)
colors1 = colorRampPalette(brewer.pal(8,'RdYlBu'))(155)
color = c(colors1[1:85],colors[75:90])
table_merged_proportion <- table_merged_proportion[rownames(anno_immune_1),]
P <- pheatmap(table_merged_proportion,color = rev(color),annotation_col = anno_sample,
              annotation_row = anno_immune_1,cluster_rows = F,clustering_distance_cols = "canberra")
P
table_merged_proportion_pca <- prcomp(t(table_merged_proportion))
table_merged_proportion_pca <- table_merged_proportion_pca$x
table_merged_proportion_pca <- table_merged_proportion_pca[,1:10]

plot(table_merged_proportion_pca[,1:2], pch=20,cex=0.7)

custom.config = umap.defaults
#custom.config$random_state = 123 ## 设定随机数
custom.config$n_neighbors = 5 ## 设定邻接数目
table_merged_proportion_umap <- umap(table_merged_proportion_pca,config = custom.config)
table_merged_proportion_kmeans <- kmeans(table_merged_proportion_pca, centers=5)
plot.data = as.data.frame(table_merged_proportion_umap$layout)
plot(plot.data, pch=20,cex=0.7)

dist <-as.matrix(dist(table_merged_proportion_pca))
edges <- mat.or.vec(0,2)

for (i in 1:nrow(dist)){
  # find closes neighbours(matches即表示最近细胞的编号)
  matches <- setdiff(order(dist[i,],decreasing = F)[1:6],i) #去除细胞自己与自己的距离
  # add edges
  edges <- rbind(edges,cbind(i,matches))  
}
edges <- na.omit(edges)

graph <- graph_from_edgelist(edges,directed=F)
graph

table_merged_proportion_cluster <- cluster_louvain(graph,resolution = 0.5)
table_merged_proportion_cluster <- as.data.frame(table_merged_proportion_cluster$membership)
rownames(table_merged_proportion_cluster) <- rownames(t(table_merged_proportion))
colnames(table_merged_proportion_cluster) <- "louvain_cluster"

plot.data = as.data.frame(table_merged_proportion_umap$layout)
colnames(plot.data)=c("umap1","umap2")
head(plot.data)
plot(plot.data,col = as.factor(table_merged_proportion_cluster$louvain_cluster), pch=20,cex=0.7)
anno_cluster <- as.character(table_merged_proportion_cluster$louvain_cluster)
ggplot(plot.data,aes(x=umap1,y=umap2,col = anno_cluster)) + geom_point(shape=19)

table_merged_proportion_kmeans_1 <- table_merged_proportion_kmeans
table_merged_proportion_kmeans_1 <- table_merged_proportion_kmeans_1$cluster
anno_cluster_kmeans <- as.character(table_merged_proportion_kmeans_1)
ggplot(plot.data,aes(x=umap1,y=umap2,col = anno_cluster_kmeans)) + geom_point(shape=19)

table_merged_proportion_cluster_ordered <- cbind(rownames(table_merged_proportion_cluster),table_merged_proportion_cluster)
table_merged_proportion_cluster_ordered <- table_merged_proportion_cluster_ordered[order(table_merged_proportion_cluster_ordered$louvain_cluster),]
anno_sample_cluster <- cbind(anno_sample[rownames(table_merged_proportion_cluster_ordered),],table_merged_proportion_cluster_ordered$louvain_cluster)
colnames(anno_sample_cluster)[3] <- "Louvain"

#table_merged_proportion_cluster_ordered_2_3 <- subset(table_merged_proportion_cluster_ordered,louvain_cluster %in% c(2,3))
#table_merged_proportion_umap_2_3 <- plot.data[rownames(table_merged_proportion_cluster_ordered_2_3),]
#table_merged_proportion_umap_2_3 <- table_merged_proportion_umap_2_3[order(table_merged_proportion_umap_2_3$umap1),]

anno_sample_cluster <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

sample_order <- anno_sample_cluster
sample_order$Louvain <- factor(sample_order$Louvain,levels = c(4,5,1,2,3))
sample_order <- sample_order[order(sample_order$Louvain),]
sample_order <- rownames(sample_order)
#anno_sample_cluster$Louvain <- as.character(anno_sample_cluster$Louvain)

pheatmap(table_merged_proportion[,sample_order],color = rev(color),
         annotation_col = anno_sample_cluster,
         annotation_row = anno_immune_1,cluster_rows = F,cluster_cols = F)

aurocs_disease <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
order <- p$tree_row$order
labels <- p$tree_row$labels
labels <- labels[order]

labels <- substr(labels,5,nchar(labels))

ann_colors = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
               colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
               colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
               colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
               colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])
names(ann_colors) <- labels
ann_colors <- list(Disease = ann_colors)
malignancy <- c("lightblue","pink")
names(malignancy) <- c("Benign","Malignant")
ann_colors$Malignancy <- malignancy
pheatmap(table_merged_proportion[,sample_order],color = rev(color),
         annotation_col = anno_sample_cluster,
         annotation_row = anno_immune_1,cluster_rows = F,cluster_cols = F,annotation_colors = ann_colors)

#plot(table_merged_proportion_pca[,1:2],col = as.factor(table_merged_proportion_cluster_ordered$louvain_cluster), pch=20,cex=0.7)
saveRDS(anno_sample_cluster,"anno_sample_cluster_extended.rds")

anno_sample_cluster <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

table_merged_proportion_merged <- as.data.frame(table_merged_proportion)
table_merged_proportion_merged <- cbind(table_merged_proportion_merged,anno_immune_1$Group)
colnames(table_merged_proportion_merged)[79] <- "cell_group"
table_merged_proportion_merged <- as.data.frame(table_merged_proportion_merged)
table_merged_proportion_merged <- table_merged_proportion_merged[c(1:38),]
result <- table_merged_proportion_merged %>%
  group_by(cell_group = .[[79]]) %>%  # 用第79列作为分组
  summarise(across(.cols = 1:78, .fns = mean, na.rm = TRUE))  # 对前78列取均值

result <- as.data.frame(result)
rownames(result) <- result$cell_group
result <- result[c(1,3,2),]
pheatmap(result[,sample_order],
         annotation_col = anno_sample_cluster,
         cluster_rows = F,cluster_cols = F,annotation_colors = ann_colors)

colors_new <- colorRampPalette(brewer.pal(8, "YlGnBu"))(90)
colors1_new = colorRampPalette(brewer.pal(8,'YlGnBu'))(155)
color_new = c(colors1_new[1:65],colors_new[82:90])

pheatmap(result[, sample_order],
         annotation_col = anno_sample_cluster,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors,
         color = rev(color_new))

table_merged_proportion_corr <- t(table_merged_proportion)
table_merged_proportion_corr <- as.data.frame(table_merged_proportion_corr)

i = 39
j = 1
corr_test <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
corr_test_corr <- corr_test$r[2,1]
corr_test_p <- corr_test$P[2,1]
corr_test <- as.data.frame(cbind(corr_test_corr,corr_test_p))
colnames(corr_test) <- c("corr","p")
rownames(corr_test) <- paste(colnames(table_merged_proportion_corr)[i],colnames(table_merged_proportion_corr)[j],sep = "_")
for (j in c(2:38,47)) {
  corr_test_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_corr <- corr_test_1$r[2,1]
  corr_test_p <- corr_test_1$P[2,1]
  corr_test_1 <- as.data.frame(cbind(corr_test_corr,corr_test_p))
  colnames(corr_test_1) <- c("corr","p")
  rownames(corr_test_1) <- paste(colnames(table_merged_proportion_corr)[i],colnames(table_merged_proportion_corr)[j],sep = "_")
  corr_test <- rbind(corr_test,corr_test_1)
}
for (i in 40:46) {
  j = 1
  corr_test_2 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_corr <- corr_test_2$r[2,1]
  corr_test_p <- corr_test_2$P[2,1]
  corr_test_2 <- as.data.frame(cbind(corr_test_corr,corr_test_p))
  colnames(corr_test_2) <- c("corr","p")
  rownames(corr_test_2) <- paste(colnames(table_merged_proportion_corr)[i],colnames(table_merged_proportion_corr)[j],sep = "_")
  for (j in c(2:38,47)) {
    corr_test_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
    corr_test_corr <- corr_test_1$r[2,1]
    corr_test_p <- corr_test_1$P[2,1]
    corr_test_1 <- as.data.frame(cbind(corr_test_corr,corr_test_p))
    colnames(corr_test_1) <- c("corr","p")
    rownames(corr_test_1) <- paste(colnames(table_merged_proportion_corr)[i],colnames(table_merged_proportion_corr)[j],sep = "_")
    corr_test_2 <- rbind(corr_test_2,corr_test_1)
  }
  corr_test <- rbind(corr_test,corr_test_2)
}

corr_test_sig <- subset(corr_test,p<0.05)

i = 39
j = 1
corr_test_mtx <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
corr_test_mtx_corr <- corr_test_mtx$r[2,1]
corr_test_mtx <- as.data.frame(corr_test_mtx_corr)
colnames(corr_test_mtx) <- colnames(table_merged_proportion_corr)[i]
rownames(corr_test_mtx) <- colnames(table_merged_proportion_corr)[j]

for (j in c(2:38,47)) {
  corr_test_mtx_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_mtx_corr <- corr_test_mtx_1$r[2,1]
  corr_test_mtx_1 <- as.data.frame(corr_test_mtx_corr)
  colnames(corr_test_mtx_1) <- colnames(table_merged_proportion_corr)[i]
  rownames(corr_test_mtx_1) <- colnames(table_merged_proportion_corr)[j]
  corr_test_mtx <- rbind(corr_test_mtx,corr_test_mtx_1)
}
for (i in 40:46) {
  j=1
  corr_test_mtx_2 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_mtx_corr <- corr_test_mtx_2$r[2,1]
  corr_test_mtx_2 <- as.data.frame(corr_test_mtx_corr)
  colnames(corr_test_mtx_2) <- colnames(table_merged_proportion_corr)[i]
  rownames(corr_test_mtx_2) <- colnames(table_merged_proportion_corr)[j]
  for (j in c(2:38,47)) {
    corr_test_mtx_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
    corr_test_mtx_corr <- corr_test_mtx_1$r[2,1]
    corr_test_mtx_1 <- as.data.frame(corr_test_mtx_corr)
    colnames(corr_test_mtx_1) <- colnames(table_merged_proportion_corr)[i]
    rownames(corr_test_mtx_1) <- colnames(table_merged_proportion_corr)[j]
    corr_test_mtx_2 <- rbind(corr_test_mtx_2,corr_test_mtx_1)
  }
  corr_test_mtx <- cbind(corr_test_mtx,corr_test_mtx_2)
}

i = 39
j = 1
corr_test_mtx_p <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
corr_test_mtx_p_corr <- corr_test_mtx_p$P[2,1]
corr_test_mtx_p <- as.data.frame(corr_test_mtx_p_corr)
colnames(corr_test_mtx_p) <- colnames(table_merged_proportion_corr)[i]
rownames(corr_test_mtx_p) <- colnames(table_merged_proportion_corr)[j]

for (j in c(2:38,47)) {
  corr_test_mtx_p_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_mtx_p_corr <- corr_test_mtx_p_1$P[2,1]
  corr_test_mtx_p_1 <- as.data.frame(corr_test_mtx_p_corr)
  colnames(corr_test_mtx_p_1) <- colnames(table_merged_proportion_corr)[i]
  rownames(corr_test_mtx_p_1) <- colnames(table_merged_proportion_corr)[j]
  corr_test_mtx_p <- rbind(corr_test_mtx_p,corr_test_mtx_p_1)
}
for (i in 40:46) {
  j=1
  corr_test_mtx_p_2 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
  corr_test_mtx_p_corr <- corr_test_mtx_p_2$P[2,1]
  corr_test_mtx_p_2 <- as.data.frame(corr_test_mtx_p_corr)
  colnames(corr_test_mtx_p_2) <- colnames(table_merged_proportion_corr)[i]
  rownames(corr_test_mtx_p_2) <- colnames(table_merged_proportion_corr)[j]
  for (j in c(2:38,47)) {
    corr_test_mtx_p_1 <- Hmisc::rcorr(as.matrix(table_merged_proportion_corr[,i]),as.matrix(table_merged_proportion_corr[,j]),type = "spearman")
    corr_test_mtx_p_corr <- corr_test_mtx_p_1$P[2,1]
    corr_test_mtx_p_1 <- as.data.frame(corr_test_mtx_p_corr)
    colnames(corr_test_mtx_p_1) <- colnames(table_merged_proportion_corr)[i]
    rownames(corr_test_mtx_p_1) <- colnames(table_merged_proportion_corr)[j]
    corr_test_mtx_p_2 <- rbind(corr_test_mtx_p_2,corr_test_mtx_p_1)
  }
  corr_test_mtx_p <- cbind(corr_test_mtx_p,corr_test_mtx_p_2)
}

pheatmap(t(corr_test_mtx),show_rownames=T,show_colnames=T,
         color=rev(colorRampPalette(brewer.pal(11,'RdBu'))(50)),
         breaks = seq(-0.5,0.5,0.02),
         display_numbers = t(matrix(ifelse(corr_test_mtx_p < 0.05, "*", ""), nrow = nrow(corr_test_mtx_p))), number_color = "white",fontsize_number = 10 
)

average_score_sample <- readRDS("/cluster3/yflu/STS/nmf/average_score_sample.rds")
average_score_sample <- as.data.frame(average_score_sample)
average_score_sample <- t(average_score_sample)
anno_sample_cluster

average_score_sample_merged <- cbind(anno_sample_cluster,average_score_sample[rownames(anno_sample_cluster),])


i=4
average_score_sample_merged_1 <- average_score_sample_merged[,c(i,1,2,3)]
colnames(average_score_sample_merged_1)[1] <- "Program_score"
average_score_sample_merged_1$Louvain <- as.character(average_score_sample_merged_1$Louvain)

P1 <- ggplot(average_score_sample_merged_1, mapping=aes(x=Louvain,y=Program_score,fill=Louvain))+ ##设置图形的纵坐标横坐标和分组
  stat_boxplot(mapping=aes(x=Louvain,y=Program_score),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
  geom_boxplot(aes(fill=Louvain),                             ##分组比较的变量
               position=position_dodge(0.8),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  ggtitle(colnames(average_score_sample_merged)[i])+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none",    
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
P1
louviannames <- unique(average_score_sample_merged_1$Louvain)
louviannames <- louviannames[order(louviannames)]
my_comparisons <- combn(louviannames,2,simplify = F)
my_comparisons_sig <- list()
for (i in 1:length(my_comparisons)) {
  por1 <- subset(average_score_sample_merged_1,Louvain == my_comparisons[[i]][1])$Program_score
  por2 <- subset(average_score_sample_merged_1,Louvain == my_comparisons[[i]][2])$Program_score
  
  if(length(por1) > 1 & length(por2) > 1){
    test <- t.test(por1, 
                   por2)
    if(is.na(test$p.value)) {
      test$p.value <- 1
    }
    if(test$p.value < 0.05){
      my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[i]]))
    }
  }
}

P2 <- P1 + stat_compare_means(comparisons=my_comparisons_sig,
                              label.y = seq(from=max(average_score_sample_merged_1$Program_score)+2, to=80, by=2),
                              method="t.test",
                              label="p.signif",hide.ns = T)
P2

for (i in 5:11) {
  average_score_sample_merged_1 <- average_score_sample_merged[,c(i,1,2,3)]
  colnames(average_score_sample_merged_1)[1] <- "Program_score"
  average_score_sample_merged_1$Louvain <- as.character(average_score_sample_merged_1$Louvain)
  
  P1 <- ggplot(average_score_sample_merged_1, mapping=aes(x=Louvain,y=Program_score,fill=Louvain))+ ##设置图形的纵坐标横坐标和分组
    stat_boxplot(mapping=aes(x=Louvain,y=Program_score),
                 geom ="errorbar",                             ##添加箱子的bar为最大、小值
                 width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
    geom_boxplot(aes(fill=Louvain),                             ##分组比较的变量
                 position=position_dodge(0.8),                 ##因为分组比较，需设组间距
                 width=0.6,                                    ##箱子的宽度
                 outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
    geom_jitter(aes(fill=Louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
    ggtitle(colnames(average_score_sample_merged)[i])+ #设置总的标题
    theme_bw()+ #背景变为白色
    theme(legend.position="none",    
          panel.grid.major = element_blank(), #不显示网格线
          panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  P1
  louviannames <- unique(average_score_sample_merged_1$Louvain)
  louviannames <- louviannames[order(louviannames)]
  my_comparisons <- combn(louviannames,2,simplify = F)
  my_comparisons_sig <- list()
  for (i in 1:length(my_comparisons)) {
    por1 <- subset(average_score_sample_merged_1,Louvain == my_comparisons[[i]][1])$Program_score
    por2 <- subset(average_score_sample_merged_1,Louvain == my_comparisons[[i]][2])$Program_score
    
    if(length(por1) > 1 & length(por2) > 1){
      test <- t.test(por1, 
                     por2)
      if(is.na(test$p.value)) {
        test$p.value <- 1
      }
      if(test$p.value < 0.05){
        my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[i]]))
      }
    }
  }
  
  P2_1 <- P1 + stat_compare_means(comparisons=my_comparisons_sig,
                                label.y = seq(from= max(average_score_sample_merged_1$Program_score)+2, to=max(average_score_sample_merged_1$Program_score)+40, by=2),
                                method="t.test",
                                label="p.signif",hide.ns = T)
  P2 <- P2 + P2_1
}
P2

table_merged_proportion_raw <- readRDS("/cluster3/yflu/STS/cpdb/table_merged_proportion_raw.rds")
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")

table_merged_proportion_raw <- table_merged_proportion_raw[rownames(anno_sample_cluster_extended),]
table_merged_proportion_raw <- as.data.frame(table_merged_proportion_raw)
table_merged_proportion_raw$louvain <- anno_sample_cluster_extended$Louvain
table_merged_proportion_raw$louvain <- as.character(table_merged_proportion_raw$louvain)

i=1
table_merged_proportion_raw_1 <- table_merged_proportion_raw[,c(i,49)]
colnames(table_merged_proportion_raw_1)[1] <- "proportion"
P1 <- ggplot(table_merged_proportion_raw_1, mapping=aes(x=louvain,y=proportion,fill=louvain))+ ##设置图形的纵坐标横坐标和分组
  stat_boxplot(mapping=aes(x=louvain,y=proportion),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
  geom_boxplot(aes(fill=louvain),                             ##分组比较的变量
               position=position_dodge(0.8),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  ggtitle(colnames(table_merged_proportion_raw)[i])+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none",    
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
P1
diseasenames <- unique(table_merged_proportion_raw_1$louvain)
diseasenames <- diseasenames[order(diseasenames)]
my_comparisons <- combn(diseasenames,2,simplify = F)
my_comparisons_sig <- list()
for (i in 1:length(my_comparisons)) {
  por1 <- subset(table_merged_proportion_raw_1,louvain == my_comparisons[[i]][1])$proportion
  por2 <- subset(table_merged_proportion_raw_1,louvain == my_comparisons[[i]][2])$proportion
  
  if(length(por1) > 1 & length(por2) > 1){
    test <- t.test(por1, 
                   por2)
    if(is.na(test$p.value)) {
      test$p.value <- 1
    }
    if(test$p.value < 0.05){
      my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[i]]))
    }
  }
}

P2 <- P1 + stat_compare_means(comparisons=my_comparisons_sig,
                              label.y = seq(from=max(table_merged_proportion_raw_1$proportion)+0.002, to=max(table_merged_proportion_raw_1$proportion)+30, by=max(table_merged_proportion_raw_1$proportion)/10),
                              method="t.test",
                              label="p.signif",hide.ns = T)
P2

for (i in 2:48) {
  table_merged_proportion_raw_1 <- table_merged_proportion_raw[,c(i,49)]
  colnames(table_merged_proportion_raw_1)[1] <- "proportion"
  
  P1_1 <- ggplot(table_merged_proportion_raw_1,mapping=aes(x=louvain,y=proportion,fill=louvain))+ ##设置图形的纵坐标横坐标和分组
    stat_boxplot(mapping=aes(x=louvain,y=proportion),
                 geom ="errorbar",                             ##添加箱子的bar为最大、小值
                 width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
    geom_boxplot(aes(fill=louvain),                             ##分组比较的变量
                 position=position_dodge(0.8),                 ##因为分组比较，需设组间距
                 width=0.6,                                    ##箱子的宽度
                 outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
    geom_jitter(aes(fill=louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
    ggtitle(colnames(table_merged_proportion_raw)[i])+ #设置总的标题
    theme_bw()+ #背景变为白色
    theme(legend.position="none",  
          panel.grid.major = element_blank(), #不显示网格线
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  P1_1
  
  my_comparisons_sig <- list()
  for (j in 1:length(my_comparisons)) {
    por1 <- subset(table_merged_proportion_raw_1,louvain == my_comparisons[[j]][1])$proportion
    por2 <- subset(table_merged_proportion_raw_1,louvain == my_comparisons[[j]][2])$proportion
    
    if(length(por1) > 1 & length(por2) > 1){
      test <- t.test(por1, 
                     por2)
      if(is.na(test$p.value)) {
        test$p.value <- 1
      }
      if(test$p.value < 0.05){
        my_comparisons_sig <- append(my_comparisons_sig,list(my_comparisons[[j]]))
      }
    }
  }
  P2_1 <- P1_1 + stat_compare_means(comparisons=my_comparisons_sig,
                                    label.y = seq(from=max(table_merged_proportion_raw_1$proportion)+0.002, to=max(table_merged_proportion_raw_1$proportion)+30, by=(max(table_merged_proportion_raw_1$proportion)+0.01)/10),
                                    method="t.test",
                                    label="p.signif",hide.ns = F)
  P2 <- P2+P2_1
}
P2