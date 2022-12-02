# install.packages("babelgene")
# library(babelgene)
# library(data.table)

setwd("/Users/mscavino/PreprocessingComparison/")

seu_integrated <- readRDS(file = "Rds/Integration/seu_integrated.rds")
# 
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# 
# s.genes_mus.musculus <- orthologs(genes = s.genes, species = "mouse") %>% 
#   pull(symbol)
# 
# # Convert G2M genes
# g2m.genes_mus.musculus <- orthologs(genes = g2m.genes, species = "mouse") %>% 
#   pull(symbol)
# 
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


DimPlot(seu_integrated, group.by = "Phase")
DimPlot(seu_integrated_scaled, group.by = "Phase")

DimPlot(seu_integrated, reduction = "pca", group.by = "Phase")
DimPlot(seu_integrated_scaled, reduction = "pca", group.by = "Phase")


