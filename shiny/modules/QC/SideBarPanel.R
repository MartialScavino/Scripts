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


# Permet d'afficher le nombre de cellule enlevées dans le dataset par les paramètre choisis
TextSERVER <- function(input, output, session, data) {
  output$texte <- renderText({
    
    nb_cell_bad <- length(data@meta.data[which(data$percent.mt> input$mt | data$nFeature_RNA < input$features[1] | data$nFeature_RNA > input$features[2]), "Quality"])
    nb_cell_tot <- length(rownames(data@meta.data))
    paste("On enlève", nb_cell_bad, "cellules du jeu de donnée soit", round(100*nb_cell_bad/nb_cell_tot, 2), "% des cellules")
    
  })
}