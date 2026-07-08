library("tigeR")
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
library(GEOquery)
library(ggbeeswarm)
library(shapviz)
library(GSVA)

metadata <- read_h5ad("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5ad")
STS_pega <- LoadH5Seurat("/cluster3/yflu/STS/pegasus/STS_95samples_nomiro_harmony_filter20250114.h5seurat",meta.data = FALSE, misc = FALSE,assays = "RNA")
STS_pega@meta.data <- metadata$obs

bulk <- AverageExpression(STS_pega,group.by = "Channel",assays = "RNA")
bulk <- bulk$RNA
bulk <- as.data.frame(bulk)

FPKM_data <- read.table("/cluster3/yflu/STS/TIGER/GSE213065_Sarcoma_ecotyper_tpm_bulksamples.txt", header=T, sep="\t")
rownames(FPKM_data) <- FPKM_data$gene
FPKM_data <- FPKM_data[,-c(1)]

GSE213065 <- readRDS("/cluster3/yflu/STS/TIGER/GSE213065.rds")
FPKM_meta <- GSE213065$`GSE213065-GPL20301_series_matrix.txt.gz`$status
meta_1 <- read.xlsx("/cluster3/yflu/STS/TIGER/NIHMS1972087-supplement-Supplementary_Tables_1-22.xlsx","Table S20")
meta_2 <- read.xlsx("/cluster3/yflu/STS/TIGER/NIHMS1972087-supplement-Supplementary_Tables_1-22.xlsx","Table S19")
colnames(meta_1) <- meta_1[2,]
meta_1 <- meta_1[-c(1:2),]
meta_1 <- subset(meta_1,`Used for Immunotherapy Analysis*` == "Yes")
meta_1$patient <- substr(meta_1$Timepoint,1,nchar(meta_1$Timepoint)-3)

colnames(meta_2) <- meta_2[1,]
meta_2 <- meta_2[-c(1),]

meta_2$`Best Response` <- as.character(meta_2$`Best Response`)

# Replace values using logical conditions
meta_2$`Best Response` <- ifelse(
  meta_2$`Best Response` %in% c("SD", "PD", "NE"), 
  "N",   # Replace SD/PD/NE with "N"
  "R"    # Replace all others with "R"
)
colnames(meta_2)[9] <- "response_NR"
meta_2$Treatment <- rep("PRE",nrow(meta_2))
meta_1 <- meta_1 %>%
  mutate(Timepoint = str_replace_all(Timepoint, fixed("-"), "."))
samplenames <- meta_1$Timepoint
samplenames <- c(paste(samplenames,"a",sep = ""),paste(samplenames,"b",sep = ""))
FPKM_data_1 <- as.data.frame(t(FPKM_data))
FPKM_data_1$sample <- rownames(FPKM_data_1)
FPKM_data_1 <- subset(FPKM_data_1,sample %in% samplenames)
FPKM_data_1 <- t(FPKM_data_1)
FPKM_meta <- as.data.frame(colnames(FPKM_data_1))
rownames(FPKM_meta) <- FPKM_meta$`colnames(FPKM_data_1)`
colnames(FPKM_meta) <- "sample_id"
FPKM_meta$sample <- substr(FPKM_meta$sample_id,1,nchar(FPKM_meta$sample_id)-4)
FPKM_meta <- FPKM_meta %>% 
  left_join(
    meta_2 %>% select(Patient, response_NR = response_NR),
    by = c("sample" = "Patient")
  )
FPKM_meta <- FPKM_meta %>% 
  left_join(
    meta_2 %>% select(Patient, Treatment = Treatment),
    by = c("sample" = "Patient")
  )
SE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fpkm = FPKM_data_1),
  colData = FPKM_meta)

anno_sample_cluster_extended <- readRDS("/cluster3/yflu/STS/cpdb/anno_sample_cluster_extended.rds")
anno_sample_cluster_extended_malig <- subset(anno_sample_cluster_extended,Malignancy == "Malignant")

bulk_malignant <- bulk[,rownames(anno_sample_cluster_extended_malig)]

louvain_markers_tumor <- readRDS("/cluster3/yflu/STS/cpdb/louvain_markers_tumor.rds")
louvain_markers_TNK <- readRDS("/cluster3/yflu/STS/cpdb/louvain_markers_TNK.rds")
louvain_markers_M <- readRDS("/cluster3/yflu/STS/cpdb/louvain_markers_M.rds")
louvain_markers_B <- readRDS("/cluster3/yflu/STS/cpdb/louvain_markers_B.rds")

top_genes <- louvain_markers_tumor %>%
  group_by(cluster) %>%  # Group by cluster column
  slice_head(n = 50) %>%  # Take top 20 rows per cluster
  ungroup()  # Remove grouping
top_genes_1 <- louvain_markers_TNK %>%
  group_by(cluster) %>%  # Group by cluster column
  slice_head(n = 50) %>%  # Take top 20 rows per cluster
  ungroup()  # Remove grouping
top_genes <- rbind(top_genes,top_genes_1)
top_genes_1 <- louvain_markers_M %>%
  group_by(cluster) %>%  # Group by cluster column
  slice_head(n = 50) %>%  # Take top 20 rows per cluster
  ungroup()  # Remove grouping
top_genes <- rbind(top_genes,top_genes_1)
top_genes_1 <- louvain_markers_B %>%
  group_by(cluster) %>%  # Group by cluster column
  slice_head(n = 50) %>%  # Take top 20 rows per cluster
  ungroup()  # Remove grouping
top_genes <- rbind(top_genes,top_genes_1)

top_genes <- unique(top_genes$gene)

set.seed(123)  # 设置随机种子确保可重复性

# 分层抽样50%
FPKM_meta_train <- FPKM_meta %>%
  group_by(response_NR) %>%
  sample_frac(0.5) %>%
  ungroup()
FPKM_data_train <- FPKM_data_1[,FPKM_meta_train$sample_id]
FPKM_meta_valid <- subset(FPKM_meta, !(sample_id %in% FPKM_meta_train$sample_id))
FPKM_data_valid <- FPKM_data_1[,FPKM_meta_valid$sample_id]

SE_train <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fpkm = FPKM_data_train),
  colData = FPKM_meta_train)
SE_valid <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fpkm = FPKM_data_valid),
  colData = FPKM_meta_valid)

train_set <- list(SE_train)
models <- c("NB","SVM","RF","CC","ADB","LGB","LGT")
mymodel <- build_Model(Model='SVM', SE=train_set, response_NR = TRUE,feature_genes = top_genes)

test <- test_Model(mymodel,SE_valid)
test_Model(mymodel,SE_valid)
p1 <- test[[2]] + labs(title="GSE21050 50% train-valid (top 50 DEG)")

SE_meta_STS_test <- as.data.frame(colnames(bulk_malignant))
rownames(SE_meta_STS_test) <- SE_meta_STS_test$`colnames(bulk_malignant)`
colnames(SE_meta_STS_test) <- "sample_id"
SE_meta_STS_test$Louvain <- anno_sample_cluster_extended_malig$Louvain
SE_meta_STS_test$Louvain <- as.character(SE_meta_STS_test$Louvain)
SE_meta_STS_test$`response_NR` <- ifelse(
  SE_meta_STS_test$`Louvain` %in% c("1","2","3"), 
  "R",
  "N"
)
#SE_meta_STS_test$response_NR <- rep("N",nrow(SE_meta_STS_test))
#SE_meta_STS_test$response_NR[c(1)] <- "R"
SE_meta_STS_test$Treatment <- rep("PRE",nrow(SE_meta_STS_test))
bulk_malignant <- bulk_malignant[rownames(FPKM_data_1),]

SE_STS_test <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fpkm = bulk_malignant),
  colData = SE_meta_STS_test)

test_STS <- test_Model(mymodel,SE_STS_test)
test_Model(mymodel,SE_STS_test)
p2 <- test_STS[[2]] + labs(title="scRNA test (top 50 DEG)")
p1+p2

data <- dataProcess(SE_STS_test, top_genes, FALSE, TRUE, FALSE)
predictions <- stats::predict(mymodel, t(data[[1]]), type = "eps-regression")
threshold = median(predictions)
predicted_classes <- ifelse(predictions <= threshold, "N", "R")
SE_meta_STS_test$predicted_classes <- predicted_classes
SE_meta_STS_test$predictions <- predictions

SE_meta_STS_test$Louvain <- factor(SE_meta_STS_test$Louvain)
SE_meta_STS_test$predicted_classes <- factor(SE_meta_STS_test$predicted_classes)

# 创建组合可视化
global_median <- median(SE_meta_STS_test$predictions)
p3 <- ggplot(SE_meta_STS_test, aes(x = Louvain, y = predictions, fill = Louvain)) +
  # 1. 小提琴图显示分布形状
  geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
  
  # 2. 箱线图显示关键统计量
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, 
               color = "black", fill = "white") +
  
  # 3. 散点图显示个体样本（按预测结果着色）
  geom_beeswarm(aes(color = predicted_classes), 
                size = 2, alpha = 0.9, cex = 3, priority = "random") +
  # 添加全局中位线（红色虚线）
  geom_hline(yintercept = global_median, 
             color = "red", 
             linetype = "dashed", 
             size = 0.5,
             alpha = 0.8) +
  # 设置颜色和主题
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_manual(values = c("N" = "#1f77b4", "R" = "#ff7f0e")) +
  labs(title = "Predictions Distribution by Louvain Cluster (50% train-valid)",
       x = "Louvain Cluster",
       y = "Prediction Score",
       color = "Predicted Class") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())
p3
pred_response(SE=SE,method = "Weighted_mean")
(p1 | p2) / p3


models <- c("NB","SVM","RF","CC","ADB","LGB","LGT")
mymodel_all <- build_Model(Model='SVM', SE=SE, response_NR = TRUE,feature_genes = top_genes)

test <- test_Model(mymodel_all,SE_STS_test)
test_Model(mymodel_all,SE_STS_test)

data <- dataProcess(SE_valid, top_genes, FALSE, TRUE, FALSE)
predictions <- stats::predict(mymodel_all, t(data[[1]]), type = "eps-regression")
predictions <- as.data.frame(predictions)
threshold = mean(predictions)
predicted_classes <- ifelse(predictions > threshold, "N", "R")

test_STS <- test_Model(mymodel_all,SE_STS_test)
test_Model(mymodel_all,SE_STS_test)

data <- dataProcess(SE_STS_test, top_genes, FALSE, TRUE, FALSE)
predictions <- stats::predict(mymodel_all, t(data[[1]]), type = "eps-regression")

p4 <- test_STS[[2]] + labs(title="scRNA test (GSE21050 train)")

threshold = median(predictions)
predicted_classes <- ifelse(predictions <= threshold, "N", "R")
SE_meta_STS_test$predicted_classes <- predicted_classes
SE_meta_STS_test$predictions <- predictions

global_median <- median(SE_meta_STS_test$predictions)
p5 <- ggplot(SE_meta_STS_test, aes(x = Louvain, y = predictions, fill = Louvain)) +
  # 1. 小提琴图显示分布形状
  geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
  
  # 2. 箱线图显示关键统计量
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, 
               color = "black", fill = "white") +
  
  # 3. 散点图显示个体样本（按预测结果着色）
  geom_beeswarm(aes(color = predicted_classes), 
                size = 2, alpha = 0.9, cex = 3, priority = "random") +
  # 添加全局中位线（红色虚线）
  geom_hline(yintercept = global_median, 
             color = "red", 
             linetype = "dashed", 
             size = 0.5,
             alpha = 0.8) +
  # 设置颜色和主题
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_manual(values = c("N" = "#1f77b4", "R" = "#ff7f0e")) +
  labs(title = "Predictions Distribution by Louvain Cluster (GSE21050 train)",
       x = "Louvain Cluster",
       y = "Prediction Score",
       color = "Predicted Class") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())
p5
p4+p5
SE_STS_test_1 <- SE_STS_test
SE_STS_test_1$response_NR <- SE_meta_STS_test$predicted_classes
test_STS <- test_Model(mymodel_all,SE_STS_test_1)
p6 <- test_STS[[2]] + labs(title="scRNA (GSE21050 train)")
p6

saveRDS(SE_meta_STS_test,"SE_meta_STS_test.rds")

weights <- t(mymodel_all$coefs) %*% mymodel_all$SV
gene_importance <- data.frame(
  Gene = colnames(weights),
  Weight = as.numeric(weights)
) %>%
  arrange(desc(abs(Weight)))
top_genes_model <- head(gene_importance, 10)

library(kernlab)
library(fastshap)

pred_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "response")
}

data <- dataProcess(SE, top_genes, FALSE, TRUE, FALSE)

shap_values <- fastshap::explain(
  mymodel_all,
  X = as.data.frame(t(data[[1]])),  # 确保与训练数据相同格式
  pred_wrapper = pred_fun,
  nsim = 50  # 增加迭代提高准确性
)

shap_importance <- shap_values %>%
  as.data.frame() %>%
  summarise_all(~ mean(abs(.))) %>%
  t() %>%
  as.data.frame() %>%
  rename(Importance = V1) %>%
  arrange(desc(Importance))

saveRDS(shap_values,"shap_values.rds")


sv_dependence(shap_obj, v = "NRG1") +
  theme_pubr() +  # 专业主题
  font("title", size = 14, face = "bold") +
  font("xlab", size = 12) +
  font("ylab", size = 12) +
  border()  # 添加边框

shap_obj <- shapviz(shap_values, X = as.data.frame(t(data[[1]])))
sv_dependence(shap_obj, v = "NRG1") +
  theme_pubr() +  # 专业主题
  font("title", size = 14, face = "bold") +
  font("xlab", size = 12) +
  font("ylab", size = 12) +
  border()  # 添加边框

p <- sv_dependence(shap_obj, v = "NRG1", color_var = "SCD") 
ggplotly(p)

R_gs <- rownames(shap_importance)[1:10]
R_gs <- list(R_gs)
names(R_gs) <- "Responder"

gsva.res <- gsva(expr = as.matrix(bulk_malignant), 
                gset.idx.list = R_gs,
                method = 'ssgsea',
                kcdf = 'Gaussian',
                verbose = FALSE)
gsva.res <- as.data.frame(t(gsva.res))
gsva.res$Sample <- rownames(gsva.res)
gsva.res <- gsva.res[rownames(anno_sample_cluster_extended_malig),]
SE_meta_STS_test <- readRDS("/cluster3/yflu/STS/TIGER/SE_meta_STS_test.rds")
SE_meta_STS_test <- cbind(SE_meta_STS_test,gsva.res)

P1 <- ggplot(SE_meta_STS_test, mapping=aes(x=Louvain,y=Responder,fill=Louvain))+ ##设置图形的纵坐标横坐标和分组
  stat_boxplot(mapping=aes(x=Louvain,y=Responder),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
  geom_boxplot(aes(fill=Louvain),                             ##分组比较的变量
               position=position_dodge(0.8),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Louvain),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  ggtitle("Responder_score")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none",    
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())
P1
SE_meta_STS_test$Louvain <- as.character(SE_meta_STS_test$Louvain)
my_comparisons <- combn(unique(SE_meta_STS_test$Louvain),2,simplify = F)
my_comparisons_sig <- list()
for (i in 1:length(my_comparisons)) {
  por1 <- subset(SE_meta_STS_test,Louvain == my_comparisons[[i]][1])$Responder
  por2 <- subset(SE_meta_STS_test,Louvain == my_comparisons[[i]][2])$Responder
  
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
                                label.y = seq(from=max(SE_meta_STS_test$Responder), to=max(SE_meta_STS_test$Responder)+32, by=0.4),
                                method="t.test",
                                label="p.signif",hide.ns = T)
P2

P3 <- ggplot(SE_meta_STS_test, mapping=aes(x=predicted_classes,y=Responder,fill=predicted_classes))+ ##设置图形的纵坐标横坐标和分组
  stat_boxplot(mapping=aes(x=predicted_classes,y=Responder),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.8))+     ##bar宽度和组间距
  geom_boxplot(aes(fill=predicted_classes),                             ##分组比较的变量
               position=position_dodge(0.8),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=predicted_classes),shape = 21,size=1.5,position=position_dodge(0.8))+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  ggtitle("Responder_score")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none",    
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())
P3

P4 <- P3 + stat_compare_means(comparisons=combn(unique(SE_meta_STS_test$predicted_classes),2,simplify = F),
                              label.y = seq(from=max(SE_meta_STS_test$Responder), to=max(SE_meta_STS_test$Responder)+32, by=0.4),
                              method="t.test",
                              label="p.signif",hide.ns = T)
P2+P4
