setwd("/Users/mscavino/PreprocessingComparison/")


library(shiny)
library(shinyWidgets)
library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(plotly)
library(cowplot)
library(viridis)
library(DoubletFinder)

seu <- Read10X("data/filtered_feature_bc_matrix/")


data <- CreateSeuratObject(seu)


data <- NormalizeData(data,verbose = FALSE)

# Find and scale variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
data <- ScaleData(data,features = rownames(data))


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["Quality"]] <- "Good"



data <- RunPCA(data)
data <- RunTSNE(data)
data <- RunUMAP(data, dims = 1:50)





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
        tabPanel("Umap", plotOutput("umap")), 
        tabPanel("Tsne", plotOutput("tsne")),
        tabPanel("Scatter", plotOutput("scatter_feature"), plotOutput("scatter_mt")),
        tabPanel("Violin", plotOutput("violinplot")),
        tabPanel("Histogramme", plotOutput("density")),
        tabPanel("Plotly", plotlyOutput("test")),
        tabPanel("FeaturePlot", textInput("gene", NULL,value = "", placeholder = "Veuillez rentrer un gène"),
                 plotOutput("featureplot"))
        
       )
     )
   )
), #Fin de la partie QC

  tabPanel("Normalisation",
           
           
           
           
           ) # Fin de la partie normalisation

) #Fin du navBarPage

) #Fin du fluid page
################################# SERVER #######################################


server <- function(input, output, session) {
  
  
  
  
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
      geom_point() + 
      geom_hline(yintercept = input$mt, col = "#FEB078") + 
      scale_color_viridis()
    
  })
  
  
  
  # Scatter plot du nombre de gènes exprimé en fonciton du nombre de read dans la cellule
  output$scatter_feature <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
      geom_point() + 
      scale_color_viridis() +
      geom_hline(yintercept = input$features, col = "#832681")
      
    
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
    DimPlot(data, reduction = "umap", group.by = "Quality", cols = c("#440154", "#21918c"))
    
  })
  
  # TSNE du dataset groupé par cellules gardées ou non
  output$tsne <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = "tsne", group.by = "Quality", cols = c("#440154", "#21918c"))
    
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
  output$density <- renderPlot({
    
    
    ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = ..density..), bins = 200, fill = "#FFDB6D", ) + 
      geom_density(color = "#440154") + 
      geom_vline(xintercept = input$features, col = "#832681")
    
    
  })
  
  # Feature Plot à partir d'un gène donné
  
  output$featureplot <- renderPlot({
    
    if (input$gene == "")
    {
      
      
    }
    
    else {FeaturePlot(data, features = input$gene)}
    
    })
  
  
}

shinyApp(ui, server)








########################## TEST ################################################
ggplot(data@meta.data, aes(x = percent.mt)) + 
  geom_histogram(aes(y = ..density..), bins = 200, fill = "#febb81", ) + 
  geom_density(color = "#f8765c") + 
  geom_vline(xintercept = 50, col = "#d3436e")

ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(aes(y = ..density..), bins = 200, fill = "#FFDB6D", ) + 
  geom_density(color = "#440154") + 
  geom_vline(xintercept = 50, col = "#832681")



# Marche pas de ouf
# data_test <- doubletFinder_v3(seu = data,  PCs = 1:50, pK = 0.09, nExp = round(ncol(data) * 0.4))

# DimPlot(data_test, group.by = "DF.classifications_0.25_0.05_554")

# ggplot(data_test@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = DF.classifications_0.25_0.09_2216)) + 
#   geom_point()
