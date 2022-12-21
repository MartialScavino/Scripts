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
