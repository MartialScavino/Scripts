if (!require("rstudioapi")) install.packages("rstudioapi")

# Setting working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Installing packages
source("modules/dependencies.R")


# Sourcing modules
source("modules/preprocessing.R")
source("modules/qc.R")


runApp(".")
