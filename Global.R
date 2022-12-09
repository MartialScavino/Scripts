setwd("/Users/mscavino/Projet//PreprocessingComparison/")


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

seu <- Read10X("../data/filtered_feature_bc_matrix/")
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

################## PLOTLY FUNCTION #################################

vline <- function(x = 0, color = "red") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color)
  )
}

hline <- function(y = 0, color = "blue") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color)
  )
}



genes <- data.frame("gene" = rownames(data))


################################## UI ##########################################

ui <- fluidPage(
  
  # Titre de l'application
  titlePanel("Comparaison des différents paramètres de preprocessing"),
  navbarPage("My Application", 
             
        tabPanel("QC",      
             
  sidebarLayout(
      
    sidebarPanel(
      chooseSliderSkin(skin = "Shiny"),
      setSliderColor(color = c("#FEB078","#832681"),sliderId =  c(1, 2)),
      
      
      sliderInput("mt", "Choisissez le pourcentage de gène mitochondriaux maximum",
                  min = 0,
                  max = 100,
                  value = 20,
                  step = 1),
      
      sliderInput("features", "Choissisez le nombre de gène exprimé voulu dans chaque cellule",
                  min = 0,
                  max = max(data$nFeature_RNA),
                  value = c(600, 10000), 
                  step = 50),
      
      textOutput("texte")
      
      
      ), #Fin du sidebar panel (partie de gauche)
      
    mainPanel(
      tabsetPanel(
        tabPanel("Umap", plotOutput("umap"), plotOutput("umap_cluster")), 
        tabPanel("Tsne", plotOutput("tsne")),
        tabPanel("Pca", plotOutput("pca")),
        tabPanel("Scatter", plotOutput("scatter_mt"), plotOutput("scatter_feature")),
        tabPanel("Violin", plotOutput("violinplot")),
        tabPanel("Histogramme", plotOutput("density_mt"), plotOutput("density_feature")),
        tabPanel("Automatic Threshold", 
                 selectInput("distribution", label = "Distribution",
                             choices = list(
                                            "Normal/Normal" = "norm_norm",
                                            "Gamma/Normal" = "gamma_norm",
                                            "Weibull/Normal" = "wei_norm")),
                 plotOutput("cutoff")),
        tabPanel("Plotly", plotlyOutput("test")),
        
        tabPanel("FeaturePlot",
                 textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène", 
                                     local = genes, valueKey = "gene",
                                     tokens =  c(1:length(genes$gene)),
                 template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
                 plotOutput("featureplot"))
        ,) #Fin tabsetpanel
        
       )
     )
), #Fin de la partie QC

  tabPanel("Normalisation",
           
           
           # Un autre tab panel set avec d'autres entrée et d'autres plots en sortie
           
           
           
           
           ) # Fin de la partie normalisation

) #Fin du navBarPage

) #Fin du fluid page
################################# SERVER #######################################


server <- function(input, output, session) {
  
  ######################################## PARTIE QC
  
  # Violin plots
  output$violinplot <- renderPlot({
    
    a <- VlnPlot(data, "nCount_RNA") + NoLegend()
    b <- VlnPlot(data, "nFeature_RNA") + 
      geom_hline(yintercept = input$features, col = "#832681") + NoLegend()
    c <- VlnPlot(data, "percent.mt") + 
      geom_hline(yintercept = input$mt, col = '#FEB078') + NoLegend()
    plot_grid(a, b, c, nrow = 1)
    
  })
  
  
  # Scatter plot du % de mt en fonction du nombre de reads dans la cellule
  output$scatter_mt <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, percent.mt, color = nFeature_RNA)) +
      geom_point(data = data@meta.data[data$percent.mt > input$mt,], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$percent.mt <= input$mt,],) + 
      geom_hline(yintercept = input$mt, col = "#FEB078") + 
      scale_color_viridis(option = "B") + theme_light()
    
  })
  
  
  
  # Scatter plot du nombre de gènes exprimé en fonciton du nombre de read dans la cellule
  output$scatter_feature <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
      geom_point(data = data@meta.data[data$nFeature_RNA < input$features[1],], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$nFeature_RNA > input$features[2],], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$nFeature_RNA <= input$features[2] & data$nFeature_RNA >= input$features[1],]) + 
      scale_color_viridis() +
      geom_hline(yintercept = input$features, col = "#832681") + theme_light()
      
    
  })
  
  
  
  # Permet d'afficher le nombre de cellule enlevées dans le dataset par les paramètre choisis
  output$texte <- renderText({
    
    nb_cell_bad <- length(data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"])
    nb_cell_tot <- length(rownames(data@meta.data))
    paste("On enlève", nb_cell_bad, "cellules du jeu de donnée soit", round(100*nb_cell_bad/nb_cell_tot, 2), "% des cellules")
    
  })
  
  
  
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

shinyApp(ui, server)


head(data)


 ########################## TEST ################################################
ggplot(data@meta.data, aes(x = percent.mt)) + 
  geom_histogram(aes(y = ..density..), bins = 200, fill = "#FEB078", ) + 
  geom_density(color = "#f8765c") + 
  geom_vline(xintercept = 50, col = "#EAB000")

ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(aes(y = ..density..), bins = 200, fill = "#DDA0DD", ) + 
  geom_density(color = "#8B0000") + 
  geom_vline(xintercept = 3000, col = "#832681") + 
  theme_light()



# Marche pas de ouf
# data_test <- doubletFinder_v3(seu = data,  PCs = 1:50, pK = 0.09, nExp = round(ncol(data) * 0.4))

# DimPlot(data_test, group.by = "DF.classifications_0.25_0.05_554")

# ggplot(data_test@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = DF.classifications_0.25_0.09_2216)) + 
#   geom_point()


p=df[, "wei_norm"]

