# install.packages("babelgene")
# library(babelgene)
# library(data.table)

setwd("/Users/mscavino/PreprocessingComparison/")

library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(cowplot)
library(viridis)
library(sctransform)
library(babelgene)



seu_integrated <- readRDS(file = "Rds/Integration/seu_integrated.rds")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

s.genes_mus.musculus <- orthologs(genes = s.genes, species = "mouse") %>%
  pull(symbol)

# Convert G2M genes
g2m.genes_mus.musculus <- orthologs(genes = g2m.genes, species = "mouse") %>%
  pull(symbol)

# 
# seu_integrated <- CellCycleScoring(seu_integrated, s.features = s.genes_mus.musculus, g2m.features = g2m.genes_mus.musculus, set.ident = TRUE)
# 
# 
# seu_integrated <- ScaleData(seu_integrated,features = rownames(seu_integrated) ,vars.to.regress = c("S.Score", "G2M.Score"))

#saveRDS(seu_integrated, file = "Rds/Integration/seu_integrated_regressed_cell_cycle.rds")

seu_integrated_scaled <- readRDS(file = "Rds/Integration/seu_integrated_regressed_cell_cycle.rds")


seu_integrated <- CellCycleScoring(seu_integrated, s.features = s.genes_mus.musculus, g2m.features = g2m.genes_mus.musculus, set.ident = TRUE)


seu_integrated <- ScaleData(seu_integrated,features = rownames(seu_integrated)) #vars.to.regress = c("S.Score", "G2M.Score")

seu_integrated <- RunPCA(seu_integrated)
seu_integrated_scaled <- RunPCA(seu_integrated_scaled)

seu_integrated <- RunUMAP(seu_integrated, dims = 1:20)
seu_integrated_scaled <- RunUMAP(seu_integrated_scaled, dims = 1:20)


seu_integrated <- FindNeighbors(seu_integrated)
seu_integrated_scaled <- FindNeighbors(seu_integrated_scaled)

seu_integrated<- FindClusters(seu_integrated)
seu_integrated_scaled <- FindClusters(seu_integrated_scaled)



DimPlot(seu_integrated, group.by = "Phase")
DimPlot(seu_integrated_scaled, group.by = "Phase")

DimPlot(seu_integrated, split.by = "Phase")
DimPlot(seu_integrated_scaled, split.by = "Phase")

DimPlot(seu_integrated, reduction = "pca", group.by = "Phase")
DimPlot(seu_integrated_scaled, reduction = "pca", group.by = "Phase")

DimPlot(seu_integrated)
DimPlot(seu_integrated_scaled)


ggplot(seu_integrated@meta.data, aes(x = seurat_clusters, fill = Phase)) +
  geom_bar(position = "fill", width = 0.5) +
  xlab("Phase") +
  ylab("Proportion") +
  theme_light()


ggplot(seu_integrated_scaled@meta.data, aes(x = seurat_clusters, fill = Phase)) +
  geom_bar(position = "fill", width = 0.5) +
  xlab("Phase") +
  ylab("Proportion") +
  theme_light()


ggplot(seu_integrated_scaled@meta.data, aes(x = seurat_clusters, fill = Technology)) +
  geom_bar(position = "fill", width = 0.5) +
  xlab("Phase") +
  ylab("Proportion") +
  theme_light()


# df <- data.frame(list(cluster = unique(seu_integrated$seurat_clusters)))
# df$cluster <- sort(df$cluster)
# 
# freq_G1 <- sapply(1:37, function(i){
#   
#   return(sum(seu_integrated$seurat_clusters == df[i, "cluster"] & seu_integrated$Phase == "G1") / sum(seu_integrated$seurat_clusters == df[i, "cluster"]))
#   
# })
# freq_S <- sapply(1:37, function(i){
#   
#   return(sum(seu_integrated$seurat_clusters == df[i, "cluster"] & seu_integrated$Phase == "S") / sum(seu_integrated$seurat_clusters == df[i, "cluster"]))
#   
# })
# freq_G2M <- sapply(1:37, function(i){
#   
#   return(sum(seu_integrated$seurat_clusters == df[i, "cluster"] & seu_integrated$Phase == "G2M") / sum(seu_integrated$seurat_clusters == df[i, "cluster"]))
#   
# })
# 
# df[["freq_G1"]] <- freq_G1
# df[["freq_S"]] <- freq_S
# df[["freq_G2M"]] <- freq_G2M
# 
# 
# ggplot(data = df, aes(x = cluster, group = 1)) +
#   geom_line(aes(y = freq_G1, color = "G1")) +
#   geom_line(aes(y = freq_S, color = "S")) +
#   geom_line(aes(y = freq_G2M, color = "G2M")) +
#   scale_color_manual(values = c("G1" = "red", "S" = "blue", "G2M" = "green")) +
#   geom_hline(yintercept = 0.6109754, linetype = "dashed", col = "red")+ 
#   geom_hline(yintercept = 0.2824727, linetype = "dashed", col = "blue")+ 
#   geom_hline(yintercept = 0.1065519, linetype = "dashed", col = "green")+ 
#   theme_light()
# 
# test <- table(seu_integrated$Phase)
# test["G1"] / sum(test)
# test["S"] / sum(test)
# test["G2M"] / sum(test)
# 
# 
# 
# #############################@
# df <- data.frame(list(cluster = unique(seu_integrated_scaled$seurat_clusters)))
# df$cluster <- sort(df$cluster)
# 
# freq_G1 <- sapply(1:35, function(i){
#   
#   return(sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"] & seu_integrated_scaled$Phase == "G1") / sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"]))
#   
# })
# freq_S <- sapply(1:35, function(i){
#   
#   return(sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"] & seu_integrated_scaled$Phase == "S") / sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"]))
#   
# })
# freq_G2M <- sapply(1:35, function(i){
#   
#   return(sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"] & seu_integrated_scaled$Phase == "G2M") / sum(seu_integrated_scaled$seurat_clusters == df[i, "cluster"]))
#   
# })
# 
# df[["freq_G1"]] <- freq_G1
# df[["freq_S"]] <- freq_S
# df[["freq_G2M"]] <- freq_G2M
# 
# 
# ggplot(data = df, aes(x = cluster, group = 1)) +
#   geom_line(aes(y = freq_G1, color = "G1")) +
#   geom_line(aes(y = freq_S, color = "S")) +
#   geom_line(aes(y = freq_G2M, color = "G2M")) +
#   scale_color_manual(values = c("G1" = "red", "S" = "blue", "G2M" = "green")) +
#   geom_hline(yintercept = 0.6109754, linetype = "dashed", col = "red")+ 
#   geom_hline(yintercept = 0.2824727, linetype = "dashed", col = "blue")+ 
#   geom_hline(yintercept = 0.1065519, linetype = "dashed", col = "green")+ 
#   theme_light()
# 
# 
# test <- table(seu_integrated_scaled$Phase)
# test["G1"] / sum(test)
# test["S"] / sum(test)
# test["G2M"] / sum(test)

