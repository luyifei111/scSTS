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
library(scCustomize)
library(irGSEA)
library(infercnv)
library(ggpubr)
library(ggsignif)
library(corrplot)

samples <- read.xlsx("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/cnvscore.xlsx","Sheet3")
samplenames <- samples$Sample
samplenames <- samplenames[order(samplenames)]
cnvobjpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/run.final.infercnv_obj",sep = "")
cnvgenespath <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_genes.dat",sep = "")
cnvscorepath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.xlsx",sep = "")
scobjpath <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/",samplenames,".rds",sep = "")
celltype3path <- paste("/cluster3/yflu/STS/separated_orig/separated/",samplenames,"/inferCNV_ref/celltype.integrated.xlsx",sep = "")
pdfpath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.pdf",sep = "")
pngpath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_cnvscore_sum_noscale.png",sep = "")
proportionpath <- paste("/cluster3/yflu/STS/separated_orig/separated/CNVSCORE/sum_new/",samplenames,"_proportion.xlsx",sep = "")

normal.pega <- readRDS("/cluster/home/yflu/STS/normal/normal.pega_1.rds")
normal.pega$group <- paste("ref-",normal.pega$celltype,"_",normal.pega$seurat_clusters, sep = "")
normal.pega$barcode <- rownames(normal.pega@meta.data)
LE_ref <- readRDS("~/STS/normal/LE_ref.rds")
SC_ref <- readRDS("~/STS/normal/SC_ref.rds")
normal.pega <- merge(normal.pega,c(SC_ref,LE_ref))
normal.pega$group2 <- rep("Normal",length(colnames(normal.pega)))

for (i in 1:length(samplenames)) {
  CNV_genes <- read.delim(cnvgenespath[i])
  #CNV_genes_1 <- subset(CNV_genes,cell_group_name %in% paste("all_observations.all_observations",c("1.1.1.1","1.1.1.2","1.1.2.1","1.1.2.2","1.2.1.1","1.2.1.2",
  #                                                                                                 "1.2.2.1","1.2.2.2"),sep = "."))
  CNV_genes_1 <- subset(CNV_genes, !(state == 3))
  CNV_genes_1 <- unique(CNV_genes_1$gene)
  
  cnvobject <- readRDS(cnvobjpath[i])
  #cnvscore <- rescale(cnvobject@expr.data, to = c(-1,1))
  cnvscore <- apply(cnvobject@expr.data[CNV_genes_1,],2,function(x){sum((x-1)^2)})
  cnvscore <- as.data.frame(cnvscore)
  
  scobj <- readRDS(scobjpath[i])
  scobj@meta.data <- scobj@meta.data[match(colnames(scobj), scobj@meta.data$barcode), ] 
  scobj$group <- Idents(scobj)
  
  celltype3 <- openxlsx::read.xlsx(celltype3path[i],"Sheet1")
  new.cluster.ids <- celltype3$celltype3
  names(new.cluster.ids) <- levels(scobj)
  
  scobj <- RenameIdents(scobj, new.cluster.ids)
  
  scobj$group2 <- Idents(scobj)
  proprotion <- table(scobj$group2)
  proprotion <- proprotion[2]/(proprotion[1]+proprotion[2])
  names(proprotion) <- samplenames[i]
  openxlsx::write.xlsx(proprotion,proportionpath[i])
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
  print(samplenames[i])
}

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
print(samplenames[i])

normal_samples <- c("T1314N","T1746","T943N","T810N","T1254N",
                    "T924T1","T924T2","T976","T1091pi","T1745",
                    "T969N","T614L","T888","T1753N","T1100N")
normal_grouppath <- paste("/cluster/home/yflu/STS/separated/CNVSCORE/",normal_samples,"_group.xlsx",sep = "")
normal_grouppath

proportion_all <- read.xlsx(proportionpath[1],"Sheet 1",header = F)
colnames(proportion_all) <- "tumor cell proportion"
rownames(proportion_all) <- samplenames[1]
proportion_all$sample <- rownames(proportion_all)
for (i in 2:length(samplenames)) {
  proportion_all_1 <- read.xlsx(proportionpath[i],"Sheet 1",header = F)
  colnames(proportion_all_1) <- "tumor cell proportion"
  rownames(proportion_all_1) <- samplenames[i]
  proportion_all_1$sample <- rownames(proportion_all_1)
  proportion_all <- rbind(proportion_all,proportion_all_1)
}

meanpath <- paste("/cluster/home/yflu/STS/separated/",samplenames,"/inferCNV_ref/",samplenames,"_cnvscore.xlsx",sep = "")
mean <- openxlsx::read.xlsx(meanpath[1],"Sheet1")
mean <- subset(mean,group2 == "Tumor")
mean <- mean(mean$cnvscore)
mean <- as.data.frame(mean)
rownames(mean) <- samplenames[1]
mean$sample <- rownames(mean)

for (i in 2:length(samplenames)) {
  mean_1 <- openxlsx::read.xlsx(meanpath[i],"Sheet1")
  mean_1 <- subset(mean_1,group2 == "Tumor")
  mean_1 <- mean(mean_1$cnvscore)
  mean_1 <- as.data.frame(mean_1)
  colnames(mean_1) <- "mean"
  rownames(mean_1) <- samplenames[i]
  mean_1$sample <- rownames(mean_1)
  mean <- rbind(mean,mean_1)
}

mean_sum <- openxlsx::read.xlsx(cnvscorepath[1],"Sheet 1")
mean_sum <- subset(mean_sum,group2 == "Tumor")
mean_sum <- mean(mean_sum$cnvscore)
mean_sum <- as.data.frame(mean_sum)
rownames(mean_sum) <- samplenames[1]
mean_sum$sample <- rownames(mean_sum)

for (i in 2:length(samplenames)) {
  mean_sum_1 <- openxlsx::read.xlsx(cnvscorepath[i],"Sheet 1")
  mean_sum_1 <- subset(mean_sum_1,group2 == "Tumor")
  mean_sum_1 <- mean(mean_sum_1$cnvscore)
  mean_sum_1 <- as.data.frame(mean_sum_1)
  colnames(mean_sum_1) <- "mean_sum"
  rownames(mean_sum_1) <- samplenames[i]
  mean_sum_1$sample <- rownames(mean_sum_1)
  mean_sum <- rbind(mean_sum,mean_sum_1)
}

proportion_score <- merge(proportion_all,c(mean,mean_sum),by = "sample")

rownames(proportion_score) <- proportion_score$sample
proportion_score <- proportion_score[,c(2,3,4)]

res1 <- cor.mtest(proportion_score, conf.level=0.95,use="complete.obs")
p <- res1$p
res2 <- cor(proportion_score,use="complete.obs" ,method = "pearson")
corrplot(res2,
         p.mat=p, 
         method = "color",
         col = rev(COL2('RdBu', 10)),
         insig ="label_sig", 
         sig.level=c(0.001, 0.01, 0.05), 
         pch.col = "black",
         tl.col = "black") 
