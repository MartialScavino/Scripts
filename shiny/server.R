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


server <- function(input, output, session) {
  
  ######################################## PARTIE QC
  
  # Violin plots
  ViolinSERVER(input, output, session, data)
  
  # Scatter plot du % de mt en fonction du nombre de reads dans la cellule et
  # Scatter plot du nombre de gènes exprimé en fonciton du nombre de read dans la cellule
  ScatterSERVER(input, output, session, data)
  
  
  
  
  
  # Permet d'afficher le nombre de cellule enlevées dans le dataset par les paramètre choisis
  TextSERVER(input, output, session, data)
  
  
  
  # UMAP du dataset groupée par cellules gardées ou non
  output$umap <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = "umap", group.by = "Quality", cols = c("#440154", "#87CEFA"))
    
  })
  
  output$umap_cluster <- renderPlot({
    
    DimPlot(data, reduction = "umap", group.by = "seurat_clusters")
  })
  
  # TSNE du dataset groupé par cellules gardées ou non
  output$tsne <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = "tsne", group.by = "Quality", cols = c("#440154", "#87CEFA"))
    
  })
  
  # PCA du dataset groupée par celulles gardées ou non
  
  output$pca <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = "pca", group.by = "Quality", cols = c("#440154", "#87CEFA"))
    
  })
  
  
  
  
  # Plot 3D qui permet d'afficher 4 variables d'un coup (nCount, nFeatures, %mt et si on garde la cellule ou pas)
  output$test <- renderPlotly({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    plot_ly(data@meta.data, x = ~nCount_RNA, y = ~nFeature_RNA, z = ~percent.mt, color = ~Quality) %>% 
      add_markers() %>% 
      layout(scene = list(xaxis = list(title = 'nCount'),
                          yaxis = list(title = 'nFeature'),
                          zaxis = list(title = 'percent.mt')))
    
  })
  
  
  # Histogramme et density plot des nFeatures
  output$density_feature <- renderPlot({
    
    ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = ..density..), bins = 200, fill = "#DDA0DD", ) + 
      geom_density(color = "#8B0000") + 
      geom_vline(xintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  
  # Histogramme et density plot des % mt
  output$density_mt <- renderPlot({
    
    ggplot(data@meta.data, aes(x = percent.mt)) + 
      geom_histogram(aes(y = ..density..), bins = 200, fill = "#FEB078", ) + 
      geom_density(color = "#f8765c") + 
      geom_vline(xintercept = input$mt, col = "#800000") + theme_light()
    
  })
  
  
  # Histogramme et density plot des % mt avec les cutoff calculés automatiquements
  
  output$cutoff <- renderPlot({
    
    
    p = df[, input$distribution]
    
    param1 = list()
    param1[[p[3]]] = as.numeric(p[7])
    param1[[p[4]]] = as.numeric(p[8])
    
    
    param2 = list()
    param2[[p[5]]] = as.numeric(p[9])
    param2[[p[6]]] = as.numeric(p[10])
    
    
    if (input$distribution == "norm_norm"){ to_cutoff = normal_normal }
    else if (input$distribution == "gamma_norm"){ to_cutoff = gamma_normal}
    else { to_cutoff = weibull_normal}
    
    cut_off = cutoff(to_cutoff)
    
    ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = ..density..), bins = 100, fill = "#DDA0DD") + 
      geom_density(color = "#8B0000") + 
      stat_function(fun = p[1], n = 101, args = param1,linetype = "dashed", aes(color = "Première distribution"), alpha = 0.5, size = 0.7) +
      stat_function(fun = p[2], n = 101, args = param2, linetype = "dashed", aes(color = "Deuxieme distribution"), alpha = 0.5, size = 0.7) +
      ylim(0, 0.000578) + 
      geom_vline(xintercept = cut_off["Estimate"], col = 'black', linetype = 'dashed', size = 1) + 
      geom_text(x = 7500, y = 0.0005, label = paste0("Valeur du cutoff : ", round(cut_off["Estimate"]))) + 
      scale_color_manual(name = "statistics", values = c("Première distribution" = "blue", "Deuxieme distribution" = "red", "Cutoff inféré" = "black")) + 
      theme_light()
    
  })
  
  
  # Feature Plot à partir d'un gène donné
  
  output$featureplot <- renderPlot({
    
    if (input$gene == ""){}
    else {FeaturePlot(data, features = input$gene)}
    
  })
  
  ############################# FIN PARTIE QC 
  
  
}