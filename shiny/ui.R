setwd("/Users/mscavino/PreprocessingComparison/")


library(shiny)
library(shinyWidgets)
library(shinysky)

library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(plotly)
library(cowplot)
library(viridis)
library(DoubletFinder)
library(cutoff)
library(autothresholdr)

ui <- fluidPage(
  
  # Titre de l'application
  titlePanel("Comparaison des différents paramètres de preprocessing"),
  navbarPage("My Application", 
             
             tabPanel("QC",      
                      
                      sidebarLayout(
                        
                        sidebarPanel(
                          chooseSliderSkin(skin = "Shiny"),
                          setSliderColor(color = c("#FEB078","#832681"),sliderId =  c(1, 2)),
                          
                          
                          sliderInput("mt", "Choisissez le pourcentage de gène mitochondriaux maximum",
                                      min = 0,
                                      max = 100,
                                      value = 20,
                                      step = 1),
                          
                          sliderInput("features", "Choissisez le nombre de gène exprimé voulu dans chaque cellule",
                                      min = 0,
                                      max = max(data$nFeature_RNA),
                                      value = c(600, 10000), 
                                      step = 50),
                          
                          textOutput("texte")
                          
                          
                        ), #Fin du sidebar panel (partie de gauche)
                        
                        mainPanel(
                          tabsetPanel(
                            UmapUI, 
                            TsneUI,
                            PcaUI,
                            ScatterUI,
                            ViolinUI,
                            HistUI,
                            AutoThresholdUI,
                            PlotlyUI,
                            FeaturePlotUI) #Fin tabsetpanel
                          
                        )
                      )
             ), #Fin de la partie QC
             
             tabPanel("Normalisation",
                      
                      
                      # Un autre tab panel set avec d'autres entrée et d'autres plots en sortie
                      
                      
                      
                      
             ) # Fin de la partie normalisation
             
  ) #Fin du navBarPage
  
) #Fin du fluid page