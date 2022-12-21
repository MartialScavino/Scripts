DoubletFinderUI <- tabPanel("Doublet Finder", plotOutput("doublet_finder_mt"),
                            plotOutput("doublet_finder_features"))





DoubletFinderSERVER <- function(input, output, session, data){
  
  output$doublet_finder_mt <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, percent.mt, color = DoubletFinderPrediction)) + 
      geom_point() +
      scale_color_manual(values = c("#440154", "#87CEFA"))
    
  })
  
  
  output$doublet_finder_features <- renderPlot({
    
    ggplot(data@meta.data, aes(nCount_RNA, nFeature_RNA, color = DoubletFinderPrediction)) +
      geom_point() +
      scale_color_manual(values = c("#440154", "#87CEFA"))
    
  })
  
  
}