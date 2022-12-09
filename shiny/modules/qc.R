######################### UI ############################

data$nFeature_RNA <- as.double(data$nFeature_RNA)
keep <- sapply(colnames(data@meta.data), function(names){
  
  check = c("integer", "character")
  
  if (typeof(data@meta.data[, names]) %in% check){ return(TRUE) }
  
  else{ return(FALSE) }
  
})
choice = names(keep[keep == T])


VisualisationUI <- tabPanel("Visualisation",
                            column(6,
                                   selectInput("reduction", label = "Plan de visualisation", 
                                               choices = c("umap", "tsne", "pca"),
                                               selected = "umap")), 
                            column(6,
                                   selectInput("group_by", label = "regroupement des cellules", 
                                               choices = choice,
                                               selected = "Quality")),
                            
                            plotOutput("dimplot"))



ScatterUI <- tabPanel("Scatter", plotOutput("scatter_mt"), plotOutput("scatter_feature"))
HeatmapUI <- tabPanel("Heatmap", plotOutput("gene_heatmap"))
ViolinUI <- tabPanel("Violin", plotOutput("violinplot"))
HistUI <- tabPanel("Histogramme", plotOutput("density_mt"), plotOutput("density_feature"))


######## AUTO THRESHOLD

CutOffUI <- tabPanel("Cutoff", 
                     selectInput("distribution", label = "Distribution",
                        choices = list(
                          "Normal/Normal" = "norm_norm",
                          "Gamma/Normal" = "gamma_norm",
                          "Weibull/Normal" = "wei_norm")),
                     plotOutput("cutoff"))


AutoThresholdUI <- tabPanel("AutothresholdR", plotOutput('autothreshold'))

ComputedThresholdUI <- tabPanel("Automatic Threshold",
                            tabsetPanel(CutOffUI, AutoThresholdUI))
###################


DoubletFinderUI <- tabPanel("Doublet Finder", plotOutput("doublet_finder_mt"),
                            plotOutput("doublet_finder_features"))



PlotlyUI <- tabPanel("Plotly", plotlyOutput("test"))

# FeaturePlotUI <- tabPanel("FeaturePlot",
#                           textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
#                                               local = genes, valueKey = "gene",
#                                               tokens =  c(1:length(genes$gene)),
#                                               template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
#                           plotOutput("featureplot"))

# FeaturePlotUI <- tabPanel("FeaturePlot",
#                           selectizeInput(
#                             inputId = 'gene', label = 'Select Something',
#                             choices = NULL,
#                             selected = 1,
#                             multiple = T
#                           ),
#                           textOutput("featureplot"))


FeaturePlotSideBarUI <- sidebarPanel(
    actionButton("reset", "Reinitialiser la signature", "danger"),
    p(markdown('---')),
    p("Gènes dans la signature"),
    p("(cliquez sur un gène pour le retirer)"),
    p(markdown('---')),
    uiOutput("list_button"),
    p(HTML("<br><br><br>")),
    actionButton('plot_button', 'Calculer la figure',"primary"))
  
FeaturePlotMainUI <-mainPanel(
  textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
                      local = genes, 
                      valueKey = "gene",
                      tokens =  c(1:length(genes$gene)),
                      template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")
  ),
  
    plotOutput(outputId = "test_features"))


FeaturePlotUI <- tabPanel("FeaturePlot", FeaturePlotSideBarUI, FeaturePlotMainUI)


### Slider MT
SliderMtUI <-  sliderInput("mt", "Choisissez le pourcentage de gène mitochondriaux maximum",
                min = 0,
                max = 100,
                value = 20,
                step = 1)

SliderFeaturesUI <- sliderInput("features", "Choissisez le nombre de gène exprimé voulu dans chaque cellule",
            min = 0,
            max = max(data$nFeature_RNA),
            value = c(600, 10000), 
            step = 50)

NumberCellsOutUI <- textOutput("texte")


SideBarPanelUI <- sidebarPanel(
  
  # Options
  chooseSliderSkin(skin = "Shiny"),
  setSliderColor(color = c("#FEB078","#832681"),sliderId =  c(1, 2)),
  
  # Sidebar
  SliderMtUI,
  SliderFeaturesUI,
  NumberCellsOutUI
  
, width = 3)













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



# UMAP du dataset groupée par cellules gardées ou non et
# UMAP du dataset groupée par ce qu'on veut
DimplotSERVER <- function(input, output, session, data){
  
  output$dimplot <- renderPlot({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    DimPlot(data, reduction = input$reduction, group.by = input$group_by)
    
  })
  
}

# Plot 3D qui permet d'afficher 4 variables d'un coup (nCount, nFeatures, %mt et si on garde la cellule ou pas)
PlotlySERVER <- function(input, output, session, data){
  
  output$test <- renderPlotly({
    
    data$Quality <- "Good"
    data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"] <- "Bad"
    plot_ly(data@meta.data, x = ~nCount_RNA, y = ~nFeature_RNA, z = ~percent.mt, color = ~Quality) %>% 
      add_markers() %>% 
      layout(scene = list(xaxis = list(title = 'nCount'),
                          yaxis = list(title = 'nFeature'),
                          zaxis = list(title = 'percent.mt')))
    
  })
  
}

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


# Histogramme et density plot des % mt avec les cutoff calculés automatiquements
HistCutoffSERVER <- function(input, output, session, data){
  
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
      geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "#DDA0DD") + 
      geom_density(color = "#8B0000") + 
      stat_function(fun = c(p[1]), n = 101, args = param1,linetype = "dashed", aes(color = "Première distribution"), alpha = 0.5, linewidth = 0.7) +
      stat_function(fun = c(p[2]), n = 101, args = param2, linetype = "dashed", aes(color = "Deuxieme distribution"), alpha = 0.5, linewidth = 0.7) +
      ylim(0, 0.000578) + 
      geom_vline(xintercept = cut_off["Estimate"], col = 'black', linetype = 'dashed', linewidth = 1) + 
      geom_text(x = 7500, y = 0.0005, label = paste0("Valeur du cutoff : ", round(cut_off["Estimate"]))) + 
      scale_color_manual(name = "statistics", values = c("Première distribution" = "blue", "Deuxieme distribution" = "red", "Cutoff inféré" = "black")) + 
      theme_light()
    
  })
  
}


HistAutoThresholdSERVER <- function(input, output, session, data){
  output$autothreshold <- renderPlot({
    
    ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = ..density..), bins = 100, fill = "#DDA0DD", alpha = 0.3) + 
      geom_density(color = "#8B0000") +
      geom_vline(data = methods_df, aes(xintercept = x, color = method), size = 0.8) +
      scale_color_manual(name = "Method", values = getPalette(17)) + 
      theme_light() + 
      theme(legend.position = "bottom")
    
  }, height = 600)
  
  
}


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

#Feature Plot à partir d'un gène donné
# FeaturePlotSERVER <- function(input, output, session, data){
# 
#   output$featureplot <- renderPlot({
# 
#     if (input$gene == ""){}
# 
#     else if (length(input$gene < 2)){FeaturePlot(data, features = input$gene)}
# 
#      else{
# 
#        data <- AddModuleScore(data, features = input$gene)
#        FeaturePlot(data, features = input$Cluster1)
# 
#        }
# 
#   })
# 
# }

HeatmapSERVER <- function(input, output, session, data){
  
  
  output$gene_heatmap <- renderPlot({
    
    
    DoHeatmap(data, features = top5_markers$gene, group.by = "seurat_clusters")
    
  }, height = 800)
  
}



FeaturePlotSERVER <- function(input, output, session, data) {
  
  liste = reactiveValues(gene_name = character(0))
  
  # Ajoute un gène à la liste lorsqu'on écrit un nouveau gène
  observeEvent(input$gene,{
               
               if (input$gene != "" & input$gene %in% liste$gene_name == FALSE){
                 
                 liste$gene_name <- c(liste$gene_name, input$gene)
               }
    })
  
  
  # Créer un bouton pour chaque élément dans la liste de gène signature
  observeEvent(liste$gene_name, output$list_button <- renderUI({
    
    if (length(liste$gene_name) > 0){
      print("Création des boutons")
      print(liste$gene_name)
      tagList(
        lapply(1:length(liste$gene_name), function(i){
          
          id1 <- paste0('slider_',liste$gene_name[i])
          actionButton(id1, liste$gene_name[i], "info")
        })
      )
    }
    
    else{
      
      output$list_button <- renderUI(actionButton('bla', "Aucun gène sélectionné", "inverse"))
      
    }
  }))
  
  # Reset list quand on appuie sur le bouton
  observeEvent(input$reset, liste$gene_name <- character(0))
  
  
  # Make button removable when clicking on it
  observe({
    
    print("Début de suppression des boutons")
    print(liste$gene_name)
    lapply(1:length(liste$gene_name), function(x){
      
      id <- paste0('slider_',liste$gene_name[x])
      operator <- "$"
      real_id <- do.call(operator, list(input, id))
      
      observeEvent(real_id, {
        
        if (liste$gene_name[x] != input$gene){ # Probleme à ce niveau
        
        print("Bouton pressé")
        liste$gene_name <- liste$gene_name[-x]
        print("Fin de suppression des boutons")
        print(liste$gene_name)}})
      
    })
  })
  
  # Calcule le plot quand le bouton est sélectionné 
  # et quand le bouton change ou quand on ajoute ou enlève un gène de la liste
  observeEvent(input$plot_button,{
    
    output$test_features <- renderPlot({ 
      isolate({
        if (length(liste$gene_name) == 1){
          
          p <- FeaturePlot(data, features = liste$gene_name[1])
          return(p)
        }
        
        else if (length(liste$gene_name) > 1){
          
          data <- AddModuleScore_UCell(obj = data, features = list(Signature = liste$gene_name))
          p <- FeaturePlot(object = data, features = "Signature_UCell") + ggtitle("Signature")
          return(p)
        }
        
        else { return() }
        
      })
    })
    
    
  })
  
}   