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
  geom_point() + scale_color_viridis_c()









############ CUTOFF ENTRE 2 DISTRIBUTIONS
# Demande de choisir la distribution qui te semble juste pour les 2 pics
# ça ne marche bien qu'avec weibull (connais pas)
test <- cutoff::em(data$nFeature_RNA, "log-normal", "normal")

# hist(data$nCount_RNA, breaks = 100)
hist(data$nFeature_RNA, probability = T, breaks = 100)

lines(weibull_normal)
lines(test)
lines(gamma_normal)


  
p = df[, "wei_norm"]
ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(aes(y = ..density..), bins = 100, fill = "#DDA0DD", ) + 
  geom_density(color = "#8B0000") + 
  stat_function(fun = dweibull, n = 101, args = list(shape = weibull_normal$param[1], scale = weibull_normal$param[2]),linetype = "dashed") +
  stat_function(fun = dnorm, n = 101, args = list(mean = weibull_normal$param[3], sd = weibull_normal$param[4]), linetype = "dashed") +
  ylim(0, 0.000578) + theme_light()


cut_off <- cutoff(test)
cut_off["Estimate"]

abline(v = cut_off, col = "red")


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


getPalette = colorRampPalette(brewer.pal(17, "Set1"))
methods_df = data.frame(list(method = methods, x = auto_thresholds))

ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(aes(y = ..density..), bins = 100, fill = "#DDA0DD") + 
  geom_density(color = "#8B0000") +
  geom_vline(data = methods_df, aes(xintercept = x, color = method)) +
  scale_color_manual(name = "Method", values = getPalette(17)) + 
  theme_light()

### Les resultats ne sont pas du tout concluants


########### DOUBLET FINDER

# run DoubletFinder

nExp_poi <- round(0.15*length(colnames(data)))
data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)

colnames(data@meta.data)[length(colnames(data@meta.data))] = "DoubletFinderPrediction"

DimPlot(data, group.by = "DoubletFinderPrediction")
DimPlot(data, reduction = "pca", group.by = "DoubletFinderPrediction")
DimPlot(data, reduction = "tsne", group.by = "DoubletFinderPrediction")


ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = DoubletFinderPrediction)) +
  geom_point() +
  scale_color_manual(values = c("#440154", "#87CEFA"))

ggplot(data@meta.data, aes(nCount_RNA, percent.mt, color = DoubletFinderPrediction)) + 
  geom_point() +
  scale_color_manual(values = c("#440154", "#87CEFA"))

############################## DIMENSION REDUCTION #############################

#data@assays$RNA@data # Normalisé
#data@assays$RNA@counts # Pas normalisé

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


# # Find significant PCs
# stdv <- mouse.sample[["pca"]]@stdev
# sum.stdv <- sum(mouse.sample[["pca"]]@stdev)
# percent.stdv <- (stdv / sum.stdv) * 100
# cumulative <- cumsum(percent.stdv)
# co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
# co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
#                      percent.stdv[2:length(percent.stdv)]) > 0.1), 
#             decreasing = T)[1] + 1
# min.pc <- min(co1, co2)
# min.pc






