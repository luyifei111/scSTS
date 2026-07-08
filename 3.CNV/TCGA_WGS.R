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
library(GenomicRanges)
library(plyranges)
library(R.matlab)

files <- list.files(
  path = "/cluster3/yflu/STS/WES_CNV/TCGA-SARC/WGS",
  pattern = "\\.seg\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)

files
i=1
seg <- read.delim(files[i])
colnames(seg)[1] <- "Sample"
seg <- subset(seg,Chromosome %in% names(table(seg$Chromosome))[1:22])
seg <- subset(seg,Num_Probes > 10)

for (i in 2:length(files)) {
  seg_1 <- read.delim(files[i])
  colnames(seg_1)[1] <- "Sample"
  seg_1 <- subset(seg_1,Chromosome %in% names(table(seg_1$Chromosome))[1:22])
  seg_1 <- subset(seg_1,Num_Probes > 10)
  seg <- rbind(seg,seg_1)
  print(i)
}

seg <- seg[
  seg$Segment_Mean > -10 & seg$Segment_Mean < 10 &
    (seg$End - seg$Start) > 100,
]

write.table(seg[,],"/cluster3/yflu/STS/WES_CNV/TCGA-SARC/segments_merged_new.seg",sep = "\t",row.names = F)

mat <- readMat("/cluster3/yflu/STS/GISTIC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat")

infer.gistic <- readGistic(gisticAllLesionsFile="/cluster3/yflu/STS/WES_CNV/gistic_TCGA/all_lesions.conf_90.txt", 
                           gisticAmpGenesFile="/cluster3/yflu/STS/WES_CNV/gistic_TCGA/amp_genes.conf_90.txt", 
                           gisticDelGenesFile="/cluster3/yflu/STS/WES_CNV/gistic_TCGA/del_genes.conf_90.txt", 
                           gisticScoresFile="/cluster3/yflu/STS/WES_CNV/gistic_TCGA/scores.gistic", isTCGA=F)

gisticChromPlot(gistic=infer.gistic, markBands="all",ref.build = "hg38",fdrCutOff = 0.05)

infer.gistic <- readGistic(gisticAllLesionsFile="/cluster3/yflu/STS/WES_CNV/infer_gistic/all_lesions.conf_90.txt", 
                           gisticAmpGenesFile="/cluster3/yflu/STS/WES_CNV/infer_gistic/amp_genes.conf_90.txt", 
                           gisticDelGenesFile="/cluster3/yflu/STS/WES_CNV/infer_gistic/del_genes.conf_90.txt", 
                           gisticScoresFile="/cluster3/yflu/STS/WES_CNV/infer_gistic/scores.gistic", isTCGA=F)