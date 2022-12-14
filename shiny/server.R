server <- function(input, output, session) {
  
  
  # Violin plots
  ViolinSERVER(input, output, session, data)
  
  # Scatter plot du % de mt en fonction du nombre de reads dans la cellule et
  # Scatter plot du nombre de gènes exprimé en fonciton du nombre de read dans la cellule
  ScatterSERVER(input, output, session, data)
  
  #DimHeatMap
  HeatmapSERVER(input, output, session, data)
  
  # Permet d'afficher le nombre de cellule enlevées dans le dataset par les paramètre choisis
  TextSERVER(input, output, session, data)
  
  # Permet d'afficher les Dimplots avec les paramètres choisis
  DimplotSERVER(input, output, session, data)
  
  # Plot 3D qui permet d'afficher 4 variables d'un coup (nCount, nFeatures, %mt et si on garde la cellule ou pas)
  PlotlySERVER(input, output, session, data)
  
  
  # Histogramme et density plot des nFeatures et
  # Histogramme et density plot des % mt
  HistSERVER(input, output, session, data)
  
  # Histogramme et density plot des % mt avec les cutoff calculés automatiquements avec cutoff
  # HistCutoffSERVER(input, output, session, data)
  
  # Histogramme et density plot des % mt avec les cutoff calculés automatiquements avec autothresholdr
  # HistAutoThresholdSERVER(input, output, session, data)
  
  
  # Scatter plot pout voir les Doublets prédits par doublet finder 
  # DoubletFinderSERVER(input, output, session, data)
  
  # Feature Plot à partir d'un gène donné
  FeaturePlotSERVER(input, output, session, data)
  
  ############################# FIN PARTIE QC 
  
}