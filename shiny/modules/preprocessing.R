setwd("/Users/mscavino/Projet/PreprocessingComparison/")



seu <- Read10X(test)
data <- CreateSeuratObject(seu)


# Mets les noms des gènes en majuscule (à cause du cell cycle)
rownames(data@assays$RNA@counts) <- toupper(rownames(data@assays$RNA@counts))
rownames(data@assays$RNA@data) <- toupper(rownames(data@assays$RNA@data))

# Cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


data <- NormalizeData(data,verbose = FALSE)

# Find and scale variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

data <- ScaleData(data,features = rownames(data)) #vars.to.regress = c("S.Score", "G2M.Score")


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["Quality"]] <- "Good"


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


genes <- data.frame(list(gene = rownames(data)))

#################### AUTOTHRESHOLD
methods <- c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", 
             "MaxEntropy", "Mean", "MinErrorI", "Minimum",
             "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag",
             "Triangle", "Yen")


auto_thresholds <- sapply(methods,
                          function(x) tryCatch(autothresholdr::auto_thresh(data$nFeature_RNA, x), 
                                               error=function(e) NA))

methods_df <- data.frame(method = methods, x = auto_thresholds)
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))




####################### MARKERS
# markers_by_clusters <- FindAllMarkers(data)
# saveRDS(markers_by_clusters, file = "Rds/Allmarkers")

markers_by_clusters <- readRDS("Rds/Shiny/Allmarkers")
top5_markers <- markers_by_clusters %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC)



######################@ Doublet Finder
nExp_poi <- round(0.15*length(colnames(data)))
data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)

colnames(data@meta.data)[length(colnames(data@meta.data))] = "DoubletFinderPrediction"




data$nFeature_RNA <- as.double(data$nFeature_RNA)
keep <- sapply(colnames(data@meta.data), function(names){
  
  check = c("integer", "character")
  
  if (typeof(data@meta.data[, names]) %in% check){ return(TRUE) }
  
  else{ return(FALSE) }
  
})
choice = names(keep[keep == T])

