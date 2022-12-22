if (!require("pacman")) install.packages("pacman")



if (!require("shinysky")) devtools::install_github("AnalytixWare/ShinySky")
if (!require("DoubletFinder")) remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

packages <- c("BiocManager","shiny", "shinyWidgets", "shiny.fluent", "shinyjs", "shinyFiles",
              "Seurat", "tidyverse", "rstudioapi", "dplyr", "stringr", "plotly", "cowplot",
              "viridis", "cutoff", "autothresholdr", "RColorBrewer", "data.table",
              "UCell", "svDialogs", "SCpubr", "ggplotify")

pacman::p_load(char = packages)

