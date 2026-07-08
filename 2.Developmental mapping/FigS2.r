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
STS.pega.development <- LoadH5Seurat("/cluster3/yflu/STS/development/data_development_umap.h5seurat",meta.data = FALSE, misc = FALSE)
STS.pega.development@meta.data <- metadata$obs

DimPlot(STS.pega.development,group.by = "celltype",raster = T,label = T)
DimPlot(STS.pega.development,group.by = "louvain_labels",raster = T,label = T,reduction = "tsne")
DimPlot(STS.pega.development,group.by = "celltype",raster = T,label = T,reduction = "umap")

Idents(STS.pega.development) <- STS.pega.development$louvain_labels

STS.pega.development$prefix <-
  sub("_.*", "", STS.pega.development$celltype)

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

order <- names(table(STS.pega.development$celltype_new))[c(1,7,5,18,8,12,4,10,6,3,13,2,16,9,17,11,14,15)]
order <- order[c(2:4,14,15,12,13,6,11,7:10,1,5,16:18)]

STS.pega.development$celltype_new <-
  factor(STS.pega.development$celltype_new, levels = order)

# 新列：大类
STS.pega.development$celltype_group <- NA

STS.pega.development$celltype_group[
  STS.pega.development$celltype_new %in% order[c(1:3, 8, 14, 15)]
] <- "Mesenchymal"

STS.pega.development$celltype_group[
  STS.pega.development$celltype_new %in% order[c(4, 5)]
] <- "Supporting"

STS.pega.development$celltype_group[
  STS.pega.development$celltype_new %in% order[c(6, 7)]
] <- "Endothelial"

STS.pega.development$celltype_group[
  STS.pega.development$celltype_new %in% order[c(9:13)]
] <- "Myogenetic"

STS.pega.development$celltype_group[
  STS.pega.development$celltype_new %in% order[c(16:18)]
] <- "Neural"

cols = c(colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(1:3)],
         colorRampPalette(brewer.pal(12,'Set3')[c(2,3)])(8)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(8,9)])(4)[c(1:2)],
         colorRampPalette(brewer.pal(12,'Set3')[c(7,8)])(8)[c(1)],
         colorRampPalette(brewer.pal(12,'Set3')[c(5,6)])(20)[c(1:5)],
         colorRampPalette(brewer.pal(12,'Set3')[c(1,2)])(12)[c(9,10)],
         colorRampPalette(brewer.pal(12,'Set3')[c(3,4)])(12)[c(1:3)])

STS.pega.development$celltype_new <- factor(STS.pega.development$celltype_new,levels = order)
DimPlot(STS.pega.development, reduction = "tsne",label = F,raster = T,cols = cols,group.by = "celltype_new")

STS.pega.development_1 <- subset(STS.pega.development,dataset_celltype %in% dataset_celltype_table$Var1)

STS.pega.development_1$dataset <- "Development"

STS.pega.development_1.sce <- as.SingleCellExperiment(STS.pega.development_1)
global_hvgs = variableGenes(dat = STS.pega.development_1.sce, exp_labels = STS.pega.development_1.sce$dataset)

aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = STS.pega.development_1.sce,
                        study_id = STS.pega.development_1.sce$dataset,
                        cell_type = STS.pega.development_1.sce$celltype_new,
                        fast_version = TRUE)
saveRDS(aurocs,"aurocs_development_260114.rds")

aurocs <- as.data.frame(aurocs)
rownames(aurocs) <- substr(rownames(aurocs),13,nchar(rownames(aurocs)))
colnames(aurocs) <- substr(colnames(aurocs),13,nchar(colnames(aurocs)))
aurocs <- aurocs[order,order]

aurocs_long <- gather(as.data.frame(aurocs))
aurocs_long$x <- rep(colnames(aurocs),18)
aurocs_long$x_index <- rep(c(1:18),18)
aurocs_long$y_index <- rep(c(1:18),each = 18)
aurocs_long <- aurocs_long[aurocs_long$x_index>=aurocs_long$y_index,]

aurocs_long$x <- factor(aurocs_long$x,levels = rev(order))
aurocs_long$key <- factor(aurocs_long$key,levels = rev(order))

p= ggplot(data=aurocs_long,aes(x=x,y=key))+
  geom_tile(aes(fill=value))+
  coord_equal(clip = "off")+
  theme_minimal() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
  theme(panel.grid=element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank())
p

intersect_genes <- readRDS("/cluster3/yflu/STS/development/intersect_genes.rds")
markers <- c("PDGFRA","THY1",
             "APOD","GPX3",
             "CXCL12","LEPR",
             "RGS5","NOTCH3",
             "MCAM","CSPG4",
             "PLVAP","FLT1",
             "CCL21","TFF3",
             "TWIST1","TWIST2",
             "TNMD","THBS4",
             "EEF1G","NACA2",
             "MYOD1","MYOG",
             "MYF5","FGFR4",
             "SOX9","RPS14P3",
             "CNTNAP2","MEIS1",
             "COL3A1","SELENOH",
             "SOX2","HES5",
             "S100B","SOX10",
             "STMN2","SNCG")
order_1 <- order
order_1[2] <- "Muscle_adipose MSC"

STS.pega.development$celltype_new <- factor(STS.pega.development$celltype_new,levels = order)
P = DotPlot(STS.pega.development,features = markers,group.by = "celltype_new",cols = c("#4575B4", "#D73027"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
P

data <- P$data
order_id <- rev(levels(data$id))
data$id <- as.character(data$id)
data$id <- factor(data$id,levels = order_id)
levels(data$features.plot) <- levels(data$features.plot)

diseasenames <- unique(as.character(STS.pega.development$celltype_new))
diseasenames <- as.data.frame(cbind(diseasenames,diseasenames))
colnames(diseasenames) <- c("Var1","Var2")
rownames(diseasenames) <- paste("Development|",diseasenames$Var1,sep="")
diseasenames <- diseasenames[paste("Development|",order,sep = ""),]
names(cols) <- diseasenames$Var1
data$cols <- as.character(data$id)

for (i in 1:648) {
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
