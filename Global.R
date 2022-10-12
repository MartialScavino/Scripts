setwd("/Users/mscavino/PreprocessingComparison/")



library(shiny)
library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(plotly)


seu <- Read10X("data/filtered_feature_bc_matrix/")


data <- CreateSeuratObject(seu)


data <- NormalizeData(data,verbose = FALSE)

# Find and scale variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
data <- ScaleData(data,features = rownames(data))


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["Quality"]] <- "Good"

ui <- fluidPage(
  
  # Titre de l'application
  titlePanel("Comparaison des différents paramêtres de preprocessing"),
  
  plotOutput("violinplot"),
  
  
  fluidRow(column(6,
                  
                  
                  plotOutput("scatter_mt"),
                  
                  sliderInput("features", "Choissisez le nombre de gène minimum exprimé dans chaque cellule",
                              min = 0,
                              max = max(data$nFeature_RNA),
                              value = 600, 
                              step = 1),
                  
                  ),
           
           column(6,
                  
                  plotOutput("scatter_feature"),
                  
                  
                  sliderInput("mt", "Choisissez le pourcentage de gène mitochondriaux maximum",
                              min = 0,
                              max = 100,
                              value = 20,
                              step = 1)
                  
                  ),
           ),
      
      plotOutput("umap"),
  
  
    plotlyOutput("test")
      
)


server <- function(input, output, session) {
          
    output$violinplot <- renderPlot({
      
      VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      
    })
    
    
    output$scatter_mt <- renderPlot({
      
      ggplot(data@meta.data, aes(nCount_RNA, percent.mt, color = nFeature_RNA)) +
        geom_point() + 
        scale_fill_viridis_d()
      
    })
    
    
    output$scatter_feature <- renderPlot({
      
      ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
        geom_point() + 
        scale_fill_viridis_d()
      
    })
    
    
    
    
    
    output$umap <- renderPlot({
      
      data$Quality <- "Good"
      data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features), "Quality"] <- "Bad"
      DimPlot(data, reduction = "umap", group.by = "Quality")
      
    })
    
    
    output$test <- renderPlotly({
      
      
      data$Quality <- "Good"
      data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features), "Quality"] <- "Bad"
      
      plot_ly(data@meta.data, x = ~nCount_RNA, y = ~nFeature_RNA, z = ~percent.mt, color = ~Quality) %>% 
        add_markers() %>% 
        layout(scene = list(xaxis = list(title = 'nCount'),
                            yaxis = list(title = 'nFeature'),
                            zaxis = list(title = 'percent.mt')))
      
    })
  
}

shinyApp(ui, server)




