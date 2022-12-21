HistUI <- tabPanel("Histogramme", 
                   plotOutput("density_mt"), 
                   plotOutput("density_feature"))

# Histogramme et density plot des nFeatures et
# Histogramme et density plot des % mt
HistSERVER <- function(input, output, session, data){
  
  output$density_feature <- renderPlot({
    
    ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = ..density..), bins = 200, fill = "#DDA0DD", ) + 
      geom_density(color = "#8B0000") + 
      geom_vline(xintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  output$density_mt <- renderPlot({
    
    ggplot(data@meta.data, aes(x = percent.mt)) + 
      geom_histogram(aes(y = ..density..), bins = 200, fill = "#FEB078", ) + 
      geom_density(color = "#f8765c") + 
      geom_vline(xintercept = input$mt, col = "#800000") + theme_light()
    
  })
}
