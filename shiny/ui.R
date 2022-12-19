setwd("/Users/mscavino/Projet/PreprocessingComparison/")


ui <- fluidPage(
                # Sert à faire une animation sur le bouton du feature plot pendant que ça charge
                useShinyjs(), 
                headUI,
  
  # Titre de l'application
  titlePanel("Comparaison des différents paramètres de preprocessing"),
  navbarPage("My Application",
             tabPanel("Données",
                      
                      shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE),
                      
                      actionButton("prepro", "Lancer prepro")
                      
                      ),
             if (exists("test")){
             tabPanel("QC",      
                      
                      SideBarPanelUI,
                        
                        mainPanel(
                          tabsetPanel(
                            
                            VisualisationUI,
                            ScatterUI,
                            HeatmapUI,
                            ViolinUI,
                            HistUI,
                            ComputedThresholdUI,
                            DoubletFinderUI,
                            PlotlyUI,
                            FeaturePlotUI) #Fin tabsetpanel
                          
                        )
                      ) #Fin de la partie QC
             
             tabPanel("Normalisation",
                      
                      
                      # Un autre tab panel set avec d'autres entrée et d'autres plots en sortie
                      
                      
                      
                      
             ) # Fin de la partie normalisation
             }
  ) #Fin du navBarPage
  
) #Fin du fluid page