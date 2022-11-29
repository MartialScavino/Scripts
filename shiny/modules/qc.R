######################### UI ############################


UmapUI <- tabPanel("Umap", plotOutput("umap"), plotOutput("umap_cluster"))
TsneUI <- tabPanel("Tsne", plotOutput("tsne"))
PcaUI <- tabPanel("Pca", plotOutput("pca"))
ScatterUI <- tabPanel("Scatter", plotOutput("scatter_mt"), plotOutput("scatter_feature"))
ViolinUI <- tabPanel("Violin", plotOutput("violinplot"))
HistUI <- tabPanel("Histogramme", plotOutput("density_mt"), plotOutput("density_feature"))

AutoThresholdUI <- tabPanel("Automatic Threshold", 
                            selectInput("distribution", label = "Distribution",
                                        choices = list(
                                          "Normal/Normal" = "norm_norm",
                                          "Gamma/Normal" = "gamma_norm",
                                          "Weibull/Normal" = "wei_norm")),
                            plotOutput("cutoff"))

PlotlyUI <- tabPanel("Plotly", plotlyOutput("test"))

FeaturePlotUI <- tabPanel("FeaturePlot",
                          textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène", 
                                              local = genes, valueKey = "gene",
                                              tokens =  c(1:length(genes$gene)),
                                              template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
                          plotOutput("featureplot"))














######################### SERVER ############################

### Violin Plot
ViolinSERVER <- function(input, output, session, data){
  
  output$violinplot <- renderPlot({
    
    a <- VlnPlot(data, "nCount_RNA") + NoLegend()
    b <- VlnPlot(data, "nFeature_RNA") + 
      geom_hline(yintercept = input$features, col = "#832681") + NoLegend()
    c <- VlnPlot(data, "percent.mt") + 
      geom_hline(yintercept = input$mt, col = '#FEB078') + NoLegend()
    plot_grid(a, b, c, nrow = 1)
    
  })
  
}

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



# Permet d'afficher le nombre de cellule enlevées dans le dataset par les paramètre choisis
TextSERVER <- function(input, output, session, data) {
  output$texte <- renderText({
  
  nb_cell_bad <- length(data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"])
  nb_cell_tot <- length(rownames(data@meta.data))
  paste("On enlève", nb_cell_bad, "cellules du jeu de donnée soit", round(100*nb_cell_bad/nb_cell_tot, 2), "% des cellules")
  
})
}


