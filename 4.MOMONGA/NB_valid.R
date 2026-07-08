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
base_dir <- "/cluster3/yflu/STS/Drug_screen/NB"

subdirs <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[-c(1,25,42,43)]
samples <- substr(subdirs,35,nchar(subdirs))
samples <- substr(samples,1,nchar(samples)-23)

saveRDS(samples,"/cluster3/yflu/STS/Drug_screen/NB/samples_nb.rds")

seuratObj_NBAtlas_share_v20240130 <- readRDS("/cluster3/yflu/STS/Drug_screen/NB/seuratObj_NBAtlas_share_v20240130.rds")
NB_tumor <- subset(seuratObj_NBAtlas_share_v20240130,Cell_type == "Neuroendocrine") 
NB_tumor <- subset(NB_tumor,Sample %in% samples)

drug_hgnc_list_combined <- readRDS("/cluster3/yflu/STS/Drug_screen/drug_hgnc_list_combined.rds")

names(drug_hgnc_list_combined) <- gsub("[|]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(" ",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("-",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[(]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[)]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub(",",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("[/]",".",  names(drug_hgnc_list_combined))
names(drug_hgnc_list_combined) <- gsub("_",".",  names(drug_hgnc_list_combined))

NB_tumor <- irGSEA.score(object = NB_tumor, assay = "RNA", slot = "data", 
                                    seeds = 123, ncores = 20,msigdb=F, 
                                    custom = T, geneset = drug_hgnc_list_combined, method = c("AUCell", "UCell", "singscore",
                                                                                              "ssgsea")[c(1)], 
                                    kcdf = 'Gaussian')
saveRDS(NB_tumor,"/cluster3/yflu/STS/Drug_screen/NB_tumor_aucell_druggable_genesets.rds")

aucell_drugs <- AverageExpression(NB_tumor,assays = "AUCell",group.by = "Sample")
aucell_drugs <- aucell_drugs$AUCell
saveRDS(aucell_drugs,"/cluster3/yflu/STS/Drug_screen/NB_aucell_druggable_genesets.rds")
aucell_drugs_scaled <- t(aucell_drugs)
aucell_drugs_scaled <- scale(aucell_drugs_scaled)
saveRDS(aucell_drugs_scaled,"/cluster3/yflu/STS/Drug_screen/NB_aucell_druggable_genesets_scaled.rds")

NB_sample_markers <- FindAllMarkers(NB_tumor,group.by = "Sample")

seuratObj_NBAtlas_share_v20240130$group <- ifelse(seuratObj_NBAtlas_share_v20240130$Cell_type == "Neuroendocrine", "Tumor", "Normal")
NB_tumor_markers <- FindAllMarkers(seuratObj_NBAtlas_share_v20240130,group.by = "group")
saveRDS(NB_sample_markers,"/cluster3/yflu/STS/Drug_screen/NB_sample_markers.rds")
saveRDS(NB_tumor_markers,"/cluster3/yflu/STS/Drug_screen/NB_tumor_markers.rds")

samples <- readRDS("/cluster3/yflu/STS/Drug_screen/NB/samples_nb.rds")
CNV_region_path <- paste("/cluster3/yflu/STS/Drug_screen/NB/hmm/",samples,"_HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat",sep = "")
infercnv_gene_order_path <- paste("/cluster3/yflu/STS/Drug_screen/NB/",samples,"_run.final.infercnv_obj",sep = "")

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
  print(paste(i,"/58",sep = ""))
}
cnv_genes$status <- cnv_genes$CNV
cnv_genes$genes <- cnv_genes$CNV

for (i in 1:nrow(cnv_genes)) {
  cnv_genes$status[i] <- strsplit(cnv_genes$CNV[i],"_")[[1]][2]
  cnv_genes$genes[i] <- strsplit(cnv_genes$CNV[i],"_")[[1]][1]
  print(i)
}
saveRDS(cnv_genes,"/cluster3/yflu/STS/Drug_screen/NB/cnv_genes.rds")

