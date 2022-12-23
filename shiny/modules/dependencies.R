if (!require("pacman")) install.packages("pacman")



if (!require("shinysky")) devtools::install_github("AnalytixWare/ShinySky")
if (!require("DoubletFinder")) remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

packages <- c("BiocManager","shiny", "shinyWidgets", "shiny.fluent", "shinyjs", "shinyFiles",
              "Seurat", "tidyverse", "rstudioapi", "dplyr", "stringr", "plotly", "cowplot",
              "viridis", "cutoff", "autothresholdr", "RColorBrewer", "data.table",
              "UCell", "svDialogs", "SCpubr", "ggplotify")

pacman::p_load(char = packages)

df <- installed.packages()
poubelle <- sapply(packages, function(x){
  
  if (!x %in% rownames(df)){
    
    print(paste0("Le Package ", x, " n'a pas été installé correctement"))
    return(0) 
  }
  
  return(1)
  
  
})
rm(poubelle)
