setwd("/Users/mscavino/PreprocessingComparison/")


ui <- fluidPage(
  
  # Titre de l'application
  titlePanel("Comparaison des différents paramètres de preprocessing"),
  navbarPage("My Application", 
             
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
                      ), #Fin de la partie QC
             
             tabPanel("Normalisation",
                      
                      
                      # Un autre tab panel set avec d'autres entrée et d'autres plots en sortie
                      
                      
                      
                      
             ) # Fin de la partie normalisation
             
  ) #Fin du navBarPage
  
) #Fin du fluid page