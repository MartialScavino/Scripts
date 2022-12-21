
ScatterUI <- tabPanel("Scatter", 
                      plotOutput("scatter_mt"), 
                      plotOutput("scatter_feature"))


# Scatter plot du % de mt en fonction du nombre de reads dans la cellule et
# Scatter plot du nombre de gènes exprimé en fonciton du nombre de read dans la cellule
ScatterSERVER <- function(input, output, session, data) {
  output$scatter_mt <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, percent.mt, color = nFeature_RNA)) +
      geom_point(data = data@meta.data[data$percent.mt > input$mt,], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$percent.mt <= input$mt,],) + 
      geom_hline(yintercept = input$mt, col = "#FEB078") + 
      scale_color_viridis(option = "B") + theme_light()
    
  })
  
  output$scatter_feature <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
      geom_point(data = data@meta.data[data$nFeature_RNA < input$features[1],], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$nFeature_RNA > input$features[2],], alpha = 0.2) + 
      geom_point(data = data@meta.data[data$nFeature_RNA <= input$features[2] & data$nFeature_RNA >= input$features[1],]) + 
      scale_color_viridis() +
      geom_hline(yintercept = input$features, col = "#832681") + theme_light()
    
    
  })
  
}