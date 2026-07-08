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
base_dir <- "/cluster3/yflu/STS/Drug_screen/CRC"

# 列出所有一级子目录（假设每个样本一个文件夹）
# subdirs <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
# subdirs <- subdirs[-c(1:3)]
# samples <- substr(subdirs,36,nchar(subdirs))
# samples <- substr(samples,1,nchar(samples)-23)
# 
# saveRDS(samples,"/cluster3/yflu/STS/Drug_screen/CRC/samples_crc.rds")

samples <- readRDS("/cluster3/yflu/STS/Drug_screen/CRC/samples_crc.rds")

counts <- read.delim("/cluster3/yflu/STS/Drug_screen/pan-cancer/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt")
rownames(counts) <- counts$Index
counts <- counts[,-1]
crc <- CreateSeuratObject(counts)

meta <- read.delim("/cluster3/yflu/STS/Drug_screen/pan-cancer/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")
rownames(meta) <- meta$Index
rownames(meta) <- gsub("-", ".", rownames(meta))

crc@meta.data <- meta[colnames(crc),]
crc <- subset(crc,Class == "Tumor")

crc_tumor <- subset(crc,Cell_type %in% c("Epithelial cells")) 
crc_tumor <- subset(crc_tumor,Sample %in% samples)

drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))
crc_tumor <- NormalizeData(crc_tumor)
crc_tumor <- irGSEA.score(object = crc_tumor, assay = "RNA", slot = "counts", 
                         seeds = 123, ncores = 20,msigdb=F, 
                         custom = T, geneset = drug_hgnc_list_combined, method = c("AUCell", "UCell", "singscore",
                                                                                   "ssgsea")[c(1)], 
                         kcdf = 'Gaussian')
saveRDS(crc_tumor,"/cluster3/yflu/STS/Drug_screen/CRC_tumor_aucell_druggable_genesets.rds")

aucell_drugs <- AverageExpression(crc_tumor,assays = "AUCell",group.by = "Sample")
aucell_drugs <- aucell_drugs$AUCell
saveRDS(aucell_drugs,"/cluster3/yflu/STS/Drug_screen/CRC_aucell_druggable_genesets.rds")
aucell_drugs_scaled <- t(aucell_drugs)
aucell_drugs_scaled <- scale(aucell_drugs_scaled)
saveRDS(aucell_drugs_scaled,"/cluster3/yflu/STS/Drug_screen/CRC_aucell_druggable_genesets_scaled.rds")

#HB_sample_markers <- FindAllMarkers(HB_tumor,group.by = "sample")

crc$group <- ifelse(crc$Cell_type %in% c("Epithelial cells"), "Tumor", "Normal")
crc <- NormalizeData(crc)
CRC_tumor_markers <- FindAllMarkers(crc,group.by = "group")

#saveRDS(HB_sample_markers,"/cluster3/yflu/STS/Drug_screen/HB_sample_markers.rds")
saveRDS(CRC_tumor_markers,"/cluster3/yflu/STS/Drug_screen/CRC_tumor_markers.rds")

CNV_region_path <- paste("/cluster3/yflu/STS/Drug_screen/CRC/hmm/",samples,"_HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat",sep = "")
infercnv_gene_order_path <- paste("/cluster3/yflu/STS/Drug_screen/CRC/",samples,"_run.final.infercnv_obj",sep = "")

i=1
cnvobject <- readRDS(infercnv_gene_order_path[i])
gene_order <- cnvobject@gene_order
#cnvexpr <- cnvobject@expr.data
CNV_region <- read.table(CNV_region_path[i],header=T)
CNV_region <- subset(CNV_region,state != 3)
CNV_region$cnv_status <- CNV_region$state

CNV_region$cluster <- CNV_region$cell_group_name
for (j in 1:nrow(CNV_region)) {
  CNV_region$cluster[j] <- strsplit(CNV_region$cluster[j],"[.]")[[1]][1]
}
for (j in 1:nrow(CNV_region)) {
  if(CNV_region$cnv_status[j] > 3){
    CNV_region$cnv_status[j] <- "gain"
  } else {
    CNV_region$cnv_status[j] <- "loss"
  }
}
cnv_genes <- c()
for (j in unique(CNV_region$chr)) {
  cnv_chr <- subset(CNV_region,chr == j)
  gene_chr <- subset(gene_order,chr == j)
  cnv_genes_chr <- c()
  for (k in 1:nrow(cnv_chr)) {
    status <- cnv_chr$cnv_status[k]
    cnv_start <- cnv_chr$start[k]
    cnv_end <- cnv_chr$end[k]
    gene_chr_sub <- subset(gene_chr,start >= cnv_start)
    gene_chr_sub <- subset(gene_chr_sub,stop <= cnv_end)
    gene_chr_sub <- paste(rownames(gene_chr_sub),status,sep = "_")
    cnv_genes_chr <- c(cnv_genes_chr,gene_chr_sub)
  }
  cnv_genes <- c(cnv_genes,cnv_genes_chr)
}
cnv_genes <- as.data.frame(cnv_genes)
colnames(cnv_genes) <- "CNV"
cnv_genes$sample <- rep(samples[i],nrow(cnv_genes))
cnv_genes <- cnv_genes[!duplicated(cnv_genes$CNV),]

for (i in 2:length(samples)) {
  cnvobject <- readRDS(infercnv_gene_order_path[i])
  gene_order <- cnvobject@gene_order
  #cnvexpr <- cnvobject@expr.data
  CNV_region <- read.table(CNV_region_path[i],header=T)
  CNV_region <- subset(CNV_region,state != 3)
  CNV_region$cnv_status <- CNV_region$state
  CNV_region$cluster <- CNV_region$cell_group_name
  for (j in 1:nrow(CNV_region)) {
    CNV_region$cluster[j] <- strsplit(CNV_region$cluster[j],"[.]")[[1]][1]
  }
  for (j in 1:nrow(CNV_region)) {
    if(CNV_region$cnv_status[j] > 3){
      CNV_region$cnv_status[j] <- "gain"
    } else {
      CNV_region$cnv_status[j] <- "loss"
    }
  }
  cnv_genes_1 <- c()
  for (j in unique(CNV_region$chr)) {
    cnv_chr <- subset(CNV_region,chr == j)
    gene_chr <- subset(gene_order,chr == j)
    cnv_genes_chr <- c()
    for (k in 1:nrow(cnv_chr)) {
      status <- cnv_chr$cnv_status[k]
      cnv_start <- cnv_chr$start[k]
      cnv_end <- cnv_chr$end[k]
      gene_chr_sub <- subset(gene_chr,start >= cnv_start)
      gene_chr_sub <- subset(gene_chr_sub,stop <= cnv_end)
      gene_chr_sub <- paste(rownames(gene_chr_sub),status,sep = "_")
      cnv_genes_chr <- c(cnv_genes_chr,gene_chr_sub)
    }
    cnv_genes_1 <- c(cnv_genes_1,cnv_genes_chr)
  }
  cnv_genes_1 <- as.data.frame(cnv_genes_1)
  colnames(cnv_genes_1) <- "CNV"
  cnv_genes_1$sample <- rep(samples[i],nrow(cnv_genes_1))
  cnv_genes_1 <- cnv_genes_1[!duplicated(cnv_genes_1$CNV),]
  cnv_genes <- rbind(cnv_genes,cnv_genes_1)
  print(paste(i,"/19",sep = ""))
}
cnv_genes$status <- cnv_genes$CNV
cnv_genes$genes <- cnv_genes$CNV

for (i in 1:nrow(cnv_genes)) {
  cnv_genes$status[i] <- strsplit(cnv_genes$CNV[i],"_")[[1]][2]
  cnv_genes$genes[i] <- strsplit(cnv_genes$CNV[i],"_")[[1]][1]
  print(i)
}
saveRDS(cnv_genes,"/cluster3/yflu/STS/Drug_screen/CRC/cnv_genes.rds")
