setwd("/Users/mscavino/PreprocessingComparison/")


library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(PCAtools)
library(cutoff)

seu <- Read10X("data/filtered_feature_bc_matrix/")


data <- CreateSeuratObject(seu)


data <- NormalizeData(data,verbose = FALSE)

# Find and scale variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
data <- ScaleData(data,features = rownames(data))



# n = length(colnames(data))
# x <- sapply(rownames(data), function(i){ 
#   sd(data@assays$RNA@data[i,]) 
#   })
# y <- sapply(rownames(data), function(i){ 
#   sum(data@assays$RNA@data[i,])/n
#   })
# plot(y, x)

############################### QUALITY CONTROL ################################

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")


# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


FeatureScatter(data, "nCount_RNA", "percent.mt")
FeatureScatter(data, "nCount_RNA", "nFeature_RNA")


data[["Quality"]] <- "Good"


data@meta.data[which(data$percent.mt> 20 & data$nFeature_RNA < 600), "Quality"] <- "Bad"


DimPlot(data, reduction = "tsne", group.by = "Quality")
DimPlot(data, reduction = "umap", group.by = "Quality")



##################### CLUSTERING


data <- RunPCA(data)
ElbowPlot(data, 50)

data <- RunTSNE(data, dims = 1:50)
data <- RunUMAP(data, dims = 1:50)

DimPlot(data, reduction = "tsne")
DimPlot(data, reduction = "umap")



data@meta.data["UMAP_1"] <- data[["umap"]]@cell.embeddings[,1]
data@meta.data["UMAP_2"] <- data[["umap"]]@cell.embeddings[,2]

ggplot(data@meta.data, aes(UMAP_1, UMAP_2, color = nFeature_RNA)) +
  geom_point() + scale_color_viridis_b()

ggplot(data@meta.data, aes(UMAP_1, UMAP_2, color = percent.mt)) +
  geom_point() + scale_color_viridis_b()









############ CUTOFF ENTRE 2 DISTRIBUTIONS
# Demande de choisir la distribution qui te semble juste pour les 2 pics
# ça ne marche bien qu'avec weibull (connais pas)
test <- em(data$nFeature_RNA, "log-normal", "normal")

# hist(data$nCount_RNA, breaks = 100)
hist(data$nFeature_RNA, probability = T, breaks = 100)
lines(test)

cut_off <- cutoff(test)


cut_off
# weibull/normal = 1000
# normal/normal = 4560
# gamma/normal = 1200
# log-normal/normal = Probleme

# Visuellement 1000 semble être un bon cutoff. Le cutoff calculé dépend beaucoup 
# des paramètres d'entrée 

# Visiblement lorsqu'on a un mix de 2 modèles dont 1 n'est pas gaussien il n'existe
# pas de méthode toute faite dans R et c'est assez compliqué

# Ce qu'on peut faire en revanche c'est se concentrer sur la vallee plutot que
# sur les pics. C'est une méthode plus naïve et c'est pas forcement ce qui est 
# recherché

library(autothresholdr)

hist(data$nFeature_RNA, probability = T, breaks = 100)
dens <- density(data$nFeature_RNA)
lines(dens)

methods <- c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", 
             "MaxEntropy", "Mean", "MinErrorI", "Minimum",
             "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag",
             "Triangle", "Yen")


auto_thresholds <- sapply(methods,
                          function(x) tryCatch(autothresholdr::auto_thresh(data$nFeature_RNA, x), 
                                               error=function(e) NA))


auto_t


abline(v = auto_thresholds, col = "green")


### Les resultats ne sont pas du tout concluants












############################## DIMENSION REDUCTION #############################



# # PCA
# data <- RunPCA(data, npcs = 50)
# 
# pca <- data[["pca"]]
# eigValues = (pca@stdev)^2  # Variance expliquée par chaque PC
# 
# 
# chosen.elbow <- findElbowPoint(eigValues)
# 
# 
# mat = GetAssayData(data, assay = "RNA", slot = "scale.data")
# horn <- parallelPCA(mat)
# 
# ElbowPlot(data, ndims = 80) + 
#   geom_vline(xintercept = chosen.elbow, col = "red")
#   #geom_vline(xintercept = horn$n, col = "blue")









