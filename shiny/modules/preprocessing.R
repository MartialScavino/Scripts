setwd("/Users/mscavino/PreprocessingComparison/")


library(shiny)
library(shinyWidgets)
library(shinysky)

library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(plotly)
library(cowplot)
library(viridis)
library(DoubletFinder)
library(cutoff)
library(autothresholdr)

seu <- Read10X("data/filtered_feature_bc_matrix/")
data <- CreateSeuratObject(seu)

# Mets les noms des gènes en majuscule (à cause du cell cycle)
rownames(data@assays$RNA@counts) <- toupper(rownames(data@assays$RNA@counts))
rownames(data@assays$RNA@data) <- toupper(rownames(data@assays$RNA@data))


data <- NormalizeData(data,verbose = FALSE)


# Find and scale variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
data <- ScaleData(data,features = rownames(data))


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["Quality"]] <- "Good"


# Cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


# Dimension reduction
data <- RunPCA(data)
data <- RunTSNE(data)
data <- RunUMAP(data, dims = 1:50)


data <- FindNeighbors(data)
data <- FindClusters(data)

# Stocke les coordonnées des cellules sur la UMAP pour faire des ggplot
data@meta.data["UMAP_1"] <- data[["umap"]]@cell.embeddings[,1]
data@meta.data["UMAP_2"] <- data[["umap"]]@cell.embeddings[,2]

data@meta.data["TSNE_1"] <- data[["tsne"]]@cell.embeddings[,1]
data@meta.data["TSNE_2"] <- data[["tsne"]]@cell.embeddings[,2]

data@meta.data["PCA_1"] <- data[["pca"]]@cell.embeddings[,1]
data@meta.data["PCA_2"] <- data[["pca"]]@cell.embeddings[,2]

ggplot(data@meta.data, aes(TSNE_1, TSNE_2, color = nFeature_RNA)) +
  geom_point() + scale_color_viridis_b()

ggplot(data@meta.data, aes(TSNE_1, TSNE_2, color = percent.mt)) +
  geom_point() + scale_color_viridis_c()

# test <- data[["pca"]]@cell.embeddings
# 
# 
# bizarre <- data@meta.data[names(test[which(test[,1] > 15 & test[,2] < 10), "PC_1"]),]
# 
# 
# table(bizarre$Phase)
# summary(bizarre$nFeature_RNA)
# summary(bizarre$nCount_RNA)
# summary(bizarre$percent.mt)


################## CUTOFF ##########################################

weibull_normal = cutoff::em(data$nFeature_RNA, "weibull", "normal")
normal_normal = cutoff::em(data$nFeature_RNA, "normal", "normal")
gamma_normal = cutoff::em(data$nFeature_RNA, "gamma", "normal")

df = data.frame()
df["fun1", "wei_norm"] <- "dweibull"
df['fun2', "wei_norm"] <- "dnorm"

df["param1_name", "wei_norm"] <- "shape"
df["param2_name", "wei_norm"] <- "scale"
df["param3_name", "wei_norm"] <- "mean"
df["param4_name", "wei_norm"] <- "sd"

df["param1", "wei_norm"] <- weibull_normal$param[1]
df["param2", "wei_norm"] <- weibull_normal$param[2]
df["param3", "wei_norm"] <- weibull_normal$param[3]
df["param4", "wei_norm"] <- weibull_normal$param[4]


df["fun1", "norm_norm"] <- "dnorm"
df['fun2', "norm_norm"] <- "dnorm"

df["param1_name", "norm_norm"] <- "mean"
df["param2_name", "norm_norm"] <- "sd"
df["param3_name", "norm_norm"] <- "mean"
df["param4_name", "norm_norm"] <- "sd"

df["param1", "norm_norm"] <- normal_normal$param[1]
df["param2", "norm_norm"] <- normal_normal$param[2]
df["param3", "norm_norm"] <- normal_normal$param[3]
df["param4", "norm_norm"] <- normal_normal$param[4]

df["fun1", "gamma_norm"] <- "dgamma"
df['fun2', "gamma_norm"] <- "dnorm"

df["param1_name", "gamma_norm"] <- "shape"
df["param2_name", "gamma_norm"] <- "rate"
df["param3_name", "gamma_norm"] <- "mean"
df["param4_name", "gamma_norm"] <- "sd"

df["param1", "gamma_norm"] <- gamma_normal$param[1]
df["param2", "gamma_norm"] <- gamma_normal$param[2]
df["param3", "gamma_norm"] <- gamma_normal$param[3]
df["param4", "gamma_norm"] <- gamma_normal$param[4]


genes <- data.frame("gene" = rownames(data))