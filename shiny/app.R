setwd("/Users/mscavino/Projet/PreprocessingComparison/")


# Libraries
library(shiny)
library(shinyWidgets)
library(shinysky)
library(shiny.fluent)
library(shinyjs)
library(shinyFiles)

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
library(RColorBrewer)
library(data.table)
library(UCell)


# Sourcing modules


source("Scripts/shiny/modules/preprocessing.R")
source("Scripts/shiny/modules/qc.R")

runApp("Scripts/shiny")

