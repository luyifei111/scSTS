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

metadata <- read_h5ad("/cluster3/yflu/STS/development/data_development_umap.h5ad")
Convert("/cluster3/yflu/STS/development/data_development_umap.h5ad", dest = "h5seurat", overwrite = F)
f <- H5File$new("/cluster3/yflu/STS/development/data_development_umap.h5seurat", "r+")
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
STS.pega.development <- LoadH5Seurat("/cluster3/yflu/STS/development/data_development_umap.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.development@meta.data <- metadata$obs

DimPlot(STS.pega.development,group.by = "celltype",raster = T,label = T)
DimPlot(STS.pega.development,group.by = "louvain_labels",raster = T,label = T,reduction = "tsne")
DimPlot(STS.pega.development,group.by = "louvain_labels",raster = T,label = T,reduction = "umap")

celltype <- as.data.frame.array(table(STS.pega.development$celltype,STS.pega.development$louvain_labels))
celltype <- t(celltype)
celltype <- as.data.frame(celltype)
celltype$sum <- rowSums(celltype)
for (i in 1:nrow(celltype)) {
  for (j in 1:ncol(celltype)) {
    celltype[i,j] <- celltype[i,j]/celltype$sum[i]
  }
}
celltype <- t(celltype)
celltype <- celltype[-41,]
celltype <- as.data.frame(celltype)
Idents(STS.pega.development) <- STS.pega.development$louvain_labels
STS.pega.development.markers <- FindAllMarkers(STS.pega.development,only.pos = T,logfc.threshold = 0.5)

celltype <- read.xlsx("/cluster3/yflu/STS/development/integrated/celltype.xlsx","Sheet1")
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(STS.pega.development)
STS.pega.development <- RenameIdents(STS.pega.development, new.cluster.ids)
DimPlot(STS.pega.development, reduction = "umap",label = F,raster = F)
DimPlot(STS.pega.development, reduction = "tsne",label = F,raster = F)
STS.pega.development$celltype_new <- Idents(STS.pega.development)
STS.pega.development$dataset <- as.character(STS.pega.development$celltype)
STS.pega.development$dataset <- substr(STS.pega.development$dataset,1,3)
STS.pega.development$dataset_celltype <- paste(STS.pega.development$dataset, STS.pega.development$celltype_new,sep = "_")

dataset_celltype_table <- as.data.frame(table(STS.pega.development$dataset_celltype))

STS.pega.development_1 <- subset(STS.pega.development,dataset_celltype %in% dataset_celltype_table$Var1)

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5ad")
STS.pega.tumor <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.tumor@meta.data <- metadata$obs

STS.pega.tumor$celltype_new <- STS.pega.tumor@meta.data$Channel

STS.pega.tumor.sce <- as.SingleCellExperiment(STS.pega.tumor)
STS.pega.development.sce <- as.SingleCellExperiment(STS.pega.development)

samplelist <- list(STS.pega.tumor.sce,STS.pega.development.sce)
names(samplelist) <- c("STS","Development")

fused_data = mergeSCE(samplelist)

global_hvgs = variableGenes(dat = fused_data, exp_labels = fused_data$study_id)

aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = fused_data,
                        study_id = fused_data$study_id,
                        cell_type = fused_data$celltype_new,
                        fast_version = TRUE)
saveRDS(aurocs,"aurocs_disease_250210.rds")
saveRDS(aurocs,"aurocs_channel_250210.rds")

STS.pega.development_1.sce <- as.SingleCellExperiment(STS.pega.development_1)
global_hvgs = variableGenes(dat = STS.pega.development_1.sce, exp_labels = STS.pega.development_1.sce$dataset)

aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = STS.pega.development_1.sce,
                        study_id = STS.pega.development_1.sce$dataset,
                        cell_type = STS.pega.development_1.sce$celltype_new,
                        fast_version = TRUE)
saveRDS(aurocs,"aurocs_development_250227.rds")

STS.pega.development_1$dataset_celltype <- paste(STS.pega.development_1$dataset,STS.pega.development_1$celltype_new,sep = "|")
celltype_order <- as.data.frame(table(STS.pega.development_1$dataset_celltype,STS.pega.development_1$celltype_new))
celltype_order <- subset(celltype_order,Freq > 0)
celltype_order$Var2 <- factor(celltype_order$Var2,levels = cellnames$Var1)
celltype_order <- celltype_order[order(celltype_order$Var2),]

aurocs <- as.data.frame(aurocs)
aurocs <- aurocs[as.character(celltype_order$Var1),as.character(celltype_order$Var1)]

aurocs_long <- gather(as.data.frame(aurocs))
aurocs_long$x <- rep(colnames(aurocs),52)
aurocs_long$x_index <- rep(c(1:52),52)
aurocs_long$y_index <- rep(c(1:52),each = 52)
aurocs_long <- aurocs_long[aurocs_long$x_index>=aurocs_long$y_index,]

aurocs_long$x <- factor(aurocs_long$x,levels = rev(celltype_order$Var1))
aurocs_long$key <- factor(aurocs_long$key,levels = rev(celltype_order$Var1))

p= ggplot(data=aurocs_long,aes(x=x,y=key))+
  geom_tile(aes(fill=value))+
  coord_equal(clip = "off")+
  theme_minimal() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
  theme(panel.grid=element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank())
p

aurocs_disease <- readRDS("aurocs_disease_250210.rds")
aurocs_disease <- as.data.frame(aurocs_disease)
aurocs_disease <- aurocs_disease[c(1:20),c(21:38)]
#plotHeatmap(aurocs_1)
anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
samplenames <- read.csv("/cluster/home/yflu/STS/pegasus/STS_sample_select_24.5.5.csv")
sample_tumors <- subset(sample_tumors,sample_tumors$`table(STS.pega.tumor@meta.data$Channel)` > 0)
samplenames_tumor <- subset(samplenames,Sample %in% rownames(sample_tumors))
samplenames_tumor_disease <- samplenames_tumor$Disease
samplenames_tumor_disease <- as.data.frame(samplenames_tumor_disease)
rownames(samplenames_tumor_disease) <- samplenames_tumor$Sample
colnames(samplenames_tumor_disease) <- "disease"
rownames(samplenames_tumor_disease) <- paste("STS|",rownames(samplenames_tumor_disease),sep = "")

diseasecolor <- colorRampPalette(brewer.pal(n = 12, name = "Set3"))(20)
names(diseasecolor) <- names(table(samplenames_tumor_disease$disease))
ann_colors <- list(disease=diseasecolor)

p = pheatmap(aurocs_disease,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')

aurocs_channel <- readRDS("/cluster3/yflu/STS/development/aurocs_channel_250210.rds")
aurocs_channel <- as.data.frame(aurocs_channel)
aurocs_channel <- aurocs_channel[c(1:78),c(79:96)]
#plotHeatmap(aurocs_1)

pheatmap(aurocs_channel,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean',
         annotation_row = samplenames_tumor_disease,annotation_colors = ann_colors)
samplenames_tumor_disease$disease <- paste("STS|",samplenames_tumor_disease$disease,sep = "")
samplenames_tumor_disease$disease <- factor(samplenames_tumor_disease$disease,levels = rownames(aurocs_disease)[p$tree_row$order])

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:8)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(6,7)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(4,5)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(10,11)])(8)[c(1:3)])

names(cols) <- levels(samplenames_tumor_disease$disease)
ann_colors <- list(disease=cols)

pheatmap(aurocs_channel,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean',
         annotation_row = samplenames_tumor_disease,annotation_colors = ann_colors)

SaveH5Seurat(CreateSeuratObject(counts = STS.pega.development@assays$RNA@counts, project = "STS.pega.development", meta.data = STS.pega.development@meta.data),
             filename = "/cluster3/yflu/STS/development/loom/STS.pega.development.h5Seurat")
Convert("/cluster3/yflu/STS/development/loom/STS.pega.development.h5Seurat", dest = "h5ad")

colnames(aurocs_disease)[p$tree_col$order]
cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:2)],
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
cellnames <- unique(as.character(STS.pega.development$celltype_new))
cellnames <- as.data.frame(cbind(cellnames,cellnames))
colnames(cellnames) <- c("Var1","Var2")
rownames(cellnames) <- paste("Development|",cellnames$Var1,sep="")
cellnames <- cellnames[colnames(aurocs_disease)[p$tree_col$order],]
names(cols) <- cellnames$Var1

STS.pega.development$celltype_new <- factor(STS.pega.development$celltype_new,levels = cellnames$Var1)
DimPlot(STS.pega.development,group.by = "celltype_new",raster = F,reduction = "umap",label = F,cols = cols)

celltype <- levels(unique(STS.pega.development$celltype_new))
celltype[18] <- "Muscle_adipose MSC"
rb.genes <- rownames(STS.pega.development)[grep("^RP[SL]",rownames(STS.pega.development))]
mt.genes <- rownames(STS.pega.development)[grep("^MT-",rownames(STS.pega.development))]

intersect_genes <- readRDS("/cluster3/yflu/STS/development/intersect_genes.rds")

i=1
markers <- openxlsx::read.xlsx("/cluster3/yflu/STS/development/development_markers_250213.xlsx",paste(celltype[i],"|up",sep = ""))
markers <- markers[order(markers$log2FC,decreasing = T),]
markers <- subset(markers,!(featurekey %in% rb.genes))
markers <- subset(markers,!(featurekey %in% mt.genes))
markers <- subset(markers,featurekey %in% intersect_genes)
markers_select <- markers$featurekey[1:3]

for (i in 2:length(celltype)) {
  markers <- openxlsx::read.xlsx("/cluster3/yflu/STS/development/development_markers_250213.xlsx",paste(celltype[i],"|up",sep = ""))
  markers <- markers[order(markers$log2FC,decreasing = T),]
  markers <- subset(markers,!(featurekey %in% markers_select))
  markers <- subset(markers,!(featurekey %in% rb.genes))
  markers <- subset(markers,!(featurekey %in% mt.genes))
  markers <- subset(markers,featurekey %in% intersect_genes)
  markers_select_1 <- markers$featurekey[1:3]
  markers_select <- c(markers_select,markers_select_1)
  print(i)
}

markers_average <- AverageExpression(STS.pega.development,features = markers_select,group.by = "celltype_new")
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
}
ann_colors = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:2)],
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
names(ann_colors) <- levels(unique(STS.pega.development$celltype_new))
ann_colors <- list(Markers = ann_colors)

pheatmap::pheatmap(markers_average,cluster_rows = F,cluster_cols = F,annotation_col = anno,annotation_colors = ann_colors)

aurocs_disease_250210 <- readRDS("/cluster3/yflu/STS/development/aurocs_disease_250210.rds")
aurocs_disease_250210 <- aurocs_disease_250210[c(21:38),c(21:38)]

aurocs_disease_250210 <- as.data.frame(aurocs_disease_250210)
aurocs_disease_250210 <- aurocs_disease_250210[rownames(cellnames),rownames(cellnames)]

aurocs_disease_250210_long <- gather(as.data.frame(aurocs_disease_250210))
aurocs_disease_250210_long$x <- rep(colnames(aurocs_disease_250210),18)
aurocs_disease_250210_long$x_index <- rep(c(1:18),18)
aurocs_disease_250210_long$y_index <- rep(c(1:18),each = 18)
aurocs_disease_250210_long <- aurocs_disease_250210_long[aurocs_disease_250210_long$x_index>=aurocs_disease_250210_long$y_index,]

aurocs_disease_250210_long$x <- factor(aurocs_disease_250210_long$x,levels = rev(rownames(cellnames)))
aurocs_disease_250210_long$key <- factor(aurocs_disease_250210_long$key,levels = rev(rownames(cellnames)))

p= ggplot(data=aurocs_disease_250210_long,aes(x=x,y=key))+
  geom_tile(aes(fill=value))+
  coord_equal(clip = "off")+
  theme_minimal() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
  theme(panel.grid=element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank())
p

RMS_atlas_final_20240130 <- readRDS("/cluster3/yflu/STS/public_data/RMS_atlas_final_20240130.rds")
RMS_atlas_final_20240130$celltype_new <- RMS_atlas_final_20240130$name

RMS_atlas_final_20240130.sce <- as.SingleCellExperiment(RMS_atlas_final_20240130)
STS.pega.development.sce <- as.SingleCellExperiment(STS.pega.development)

samplelist <- list(RMS_atlas_final_20240130.sce,STS.pega.development.sce)
names(samplelist) <- c("RMS","Development")

fused_data = mergeSCE(samplelist)

global_hvgs = variableGenes(dat = fused_data, exp_labels = fused_data$study_id)

aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = fused_data,
                        study_id = fused_data$study_id,
                        cell_type = fused_data$celltype_new,
                        fast_version = TRUE)
saveRDS(aurocs,"aurocs_external_RMS_250304.rds")

RMS_names <- as.character(unique(RMS_atlas_final_20240130$name))
infer_path <- paste("/cluster3/yflu/STS/development/RMS_external/",RMS_names,sep = "")

normal.pega <- readRDS("/cluster/home/yflu/STS/normal/normal.pega_1.rds")
normal.pega$group <- paste("ref-",normal.pega$celltype,"_",normal.pega$seurat_clusters, sep = "")
normal.pega$group

normal.pega_2 <- subset(normal.pega, celltype %in% c("Muscle cells"))
library(infercnv)
for (i in 32:72) {
  seuobj <- subset(RMS_atlas_final_20240130,name == RMS_names[i])
  seuobj$group <- seuobj$seurat_clusters
  seuobj_merged <- merge(seuobj,normal.pega_2)
  table_cellnumber <- as.data.frame(table(seuobj_merged$group))
  seuobj_merged <- subset(seuobj_merged,group %in% as.character(subset(table_cellnumber, Freq > 1)$Var1))
  counts <- seuobj_merged@assays$RNA@counts
  anno <- data.frame(seuobj_merged@meta.data$group)
  rownames(anno) <- rownames(seuobj_merged@meta.data)
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                      annotations_file=anno,
                                      delim="\t",
                                      gene_order_file="/cluster3/yflu/RT/inferCNV/hg38_gencode_v27.txt",
                                      ref_group_names=c(names(table(normal.pega_2$group)))) 
  infercnv_obj_2 = infercnv::run(infercnv_obj,
                                 cutoff = 0.1,
                                 out_dir = infer_path[i], 
                                 cluster_by_groups = T,
                                 cluster_references = F,
                                 HMM = T,
                                 analysis_mode = c('samples'),
                                 denoise = TRUE,
                                 num_threads = 40,
                                 inspect_subclusters = F)
  gc()
  print(paste(i,"72",sep = "/"))
  write.table(i,"/cluster3/yflu/STS/development/RMS_external/order.txt")
}

cnvobjpath <- paste("/cluster3/yflu/STS/development/RMS_external/",RMS_names,"/run.final.infercnv_obj",sep = "")
cnvgenespath <- paste("/cluster3/yflu/STS/development/RMS_external/",RMS_names,"/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_genes.dat",sep = "")
cnvscorepath <- paste("/cluster3/yflu/STS/development/RMS_external/CNVSCORE/sum_new/",RMS_names,"_cnvscore_sum_noscale.xlsx",sep = "")
pdfpath <- paste("/cluster3/yflu/STS/development/RMS_external/CNVSCORE/sum_new/",RMS_names,"_cnvscore_sum_noscale.pdf",sep = "")
pngpath <- paste("/cluster3/yflu/STS/development/RMS_external/CNVSCORE/sum_new/",RMS_names,"_cnvscore_sum_noscale.png",sep = "")
proportionpath <- paste("/cluster3/yflu/STS/development/RMS_external/CNVSCORE/sum_new/",RMS_names,"_proportion.xlsx",sep = "")

for (i in c(1:26)) {
  CNV_genes <- read.delim(cnvgenespath[i])
  #CNV_genes_1 <- subset(CNV_genes,cell_group_name %in% paste("all_observations.all_observations",c("1.1.1.1","1.1.1.2","1.1.2.1","1.1.2.2","1.2.1.1","1.2.1.2",
  #                                                                                                 "1.2.2.1","1.2.2.2"),sep = "."))
  CNV_genes_1 <- subset(CNV_genes, !(state == 3))
  CNV_genes_1 <- unique(CNV_genes_1$gene)
  
  cnvobject <- readRDS(cnvobjpath[i])
  #cnvscore <- rescale(cnvobject@expr.data, to = c(-1,1))
  cnvscore <- apply(cnvobject@expr.data[CNV_genes_1,],2,function(x){sum((x-1)^2)})
  cnvscore <- as.data.frame(cnvscore)
  cnvscore$barcode <- rownames(cnvscore)
  cnvscore$group <- substr(cnvscore$barcode,1,nchar(cnvscore$barcode)-17)
  cnvscore$group <- ifelse(cnvscore$group %in% c("T1745","T1753N","T943N","T969N"),"Normal","Tumor")
  cnvscore$group2 <- cnvscore$group
  openxlsx::write.xlsx(cnvscore,cnvscorepath[i])
  
  pdf(pdfpath[i])
  plot <- ggplot(cnvscore,aes(x=group2,y=cnvscore,fill=group2)) + geom_boxplot(outlier.shape = NA) + 
    geom_signif(comparisons = list(c("Tumor","Normal")), 
                map_signif_level = TRUE, test = t.test, 
                tip_length = c(0.05,0.05)) +
    stat_summary(fun.data = function(x) data.frame(y=max(x)*0.9,label = paste("Mean=", round(mean(x), 4))), geom="text", size = 5)
  print(plot)
  dev.off()
  
  png(pngpath[i],width = 618,height = 506,units = "px")
  plot <- ggplot(cnvscore,aes(x=group2,y=cnvscore,fill=group2)) + geom_boxplot(outlier.shape = NA) + 
    geom_signif(comparisons = list(c("Tumor","Normal")), 
                map_signif_level = TRUE, test = t.test, 
                tip_length = c(0.05,0.05)) +
    stat_summary(fun.data = function(x) data.frame(y=max(x)*0.9,label = paste("Mean=", round(mean(x), 4))), geom="text", size = 5)
  print(plot)
  dev.off()
  print(i)
}
