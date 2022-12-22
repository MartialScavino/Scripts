FeaturePlotSideBarUI <- sidebarPanel(
  actionButton("reset", "Reinitialiser la signature", styleclass = "danger"),
  p(markdown('---')),
  p("Gènes dans la signature"),
  p("(cliquez sur un gène pour le retirer)"),
  uiOutput("list_button"),
  p(markdown('---')),
  p(HTML("<br><br><br>")),
  actionButton('plot_button', span('Calculer la figure', id="UpdateAnimate", class=""), styleclass = "primary"))

FeaturePlotMainUI <-mainPanel(
  textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
                      local = genes, 
                      valueKey = "gene",
                      tokens =  c(1:length(genes$gene)),
                      template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")
  ),
  plotOutput(outputId = "test_features"))


FeaturePlotUI <- tabPanel("FeaturePlot", FeaturePlotSideBarUI, FeaturePlotMainUI)






FeaturePlotSERVER <- function(input, output, session, data) {
  
  liste = reactiveValues(df = data.frame(list(gene_name = character(0), show_button = character(0))))
  
  # Ajoute un gène à la liste lorsqu'on écrit un nouveau gène
  observeEvent(input$gene,{
    
    if (input$gene != "" & input$gene %in% liste$df$gene_name == FALSE){
      liste$df[nrow(liste$df) + 1,] <- c(input$gene, TRUE)
    }
    
    else if (input$gene != "" & input$gene %in% liste$df$gene_name == TRUE){
      if (liste$df[which(liste$df$gene_name == input$gene), "show_button"] == FALSE){
        
        liste$df[which(liste$df$gene_name == input$gene), "show_button"] <- TRUE
        
      }
    }
  })
  
  
  # Créer un bouton pour chaque élément dans la liste de gène signature
  observeEvent(liste$df$gene_name, output$list_button <- renderUI({
    
    if (length(liste$df$gene_name) > 0){
      tagList(
        lapply(1:length(liste$df$gene_name), function(i){
          
          if (liste$df$show_button[i] == TRUE){
            
            id1 <- paste0('slider_',liste$df$gene_name[i])
            actionButton(id1, liste$df$gene_name[i], "info", 
                         onclick = "Shiny.onInputChange('myclick', {id : this.id, val : this})")
          }
        })
      )
    }
    
    else{
      
      output$list_button <- renderUI(actionButton('bla', "Aucun gène sélectionné", "inverse"))
      
    }
  }))
  
  # Reset list quand on appuie sur le bouton
  observeEvent(input$reset, {
    
    liste$df <- data.frame(list(gene_name = character(0), show_button = character(0)))
    
  })
  
  # Make button removable when clicking on it
  observeEvent(input$myclick, {
    
    id <- strsplit(input$myclick$id, "_")[[1]][2]
    liste$df[which(liste$df$gene_name == id), "show_button"] <- FALSE
    
  })
  
  # Calcule le plot quand le bouton est sélectionné 
  observeEvent(input$plot_button,{
    
    output$test_features <- renderPlot({
      
      isolate({
        
        liste_gene_plot <- liste$df[which(liste$df$show_button == TRUE), "gene_name"]
        
        if (length(liste_gene_plot) == 1){
          
          p <- FeaturePlot(data, features = liste_gene_plot)
          return(p)
        }
        
        else if (length(liste_gene_plot) > 1){
          
          shinyjs::addClass(id = "UpdateAnimate", class = "loading dots")
          shinyjs::disable("plot_button")
          
          data <- AddModuleScore_UCell(obj = data, features = list(Signature = liste_gene_plot))
          p <- FeaturePlot(object = data, features = "Signature_UCell") + ggtitle("Signature")
          
          shinyjs::enable("plot_button")
          shinyjs::removeClass(id = "UpdateAnimate", class = "loading dots")
          return(p)
          
        }
      
        else { return() }
      
      })
    })
  })
} 
