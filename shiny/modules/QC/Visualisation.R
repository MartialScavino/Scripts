VisualisationUI <- tabPanel("Visualisation",
                            column(6,
                                   selectInput("reduction", label = "Plan de visualisation", 
                                               choices = c("umap", "tsne", "pca"),
                                               selected = "umap")), 
                            column(6,
                                   selectInput("group_by", label = "Regroupement des cellules", 
                                               choices = choice,
                                               selected = "Quality")),
                            
                            plotOutput("dimplot"))

# UMAP du dataset groupée par cellules gardées ou non et
# UMAP du dataset groupée par ce qu'on veut
DimplotSERVER <- function(input, output, session, data){
  
  output$dimplot <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = input$reduction, group.by = input$group_by)
    
  })
  
}