HeatmapUI <- tabPanel("Heatmap", 
                      plotOutput("gene_heatmap"))



HeatmapSERVER <- function(input, output, session, data){
  
  
  output$gene_heatmap <- renderPlot({
    
    
    DoHeatmap(data, features = top5_markers$gene, group.by = "seurat_clusters")
    
  }, height = 800)
  
}
