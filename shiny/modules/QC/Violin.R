ViolinUI <- tabPanel("Violin", 
                     selectInput("group", label = "Axe x", 
                                 choices = choice[choice != "Quality"], 
                                 selected = "orig.ident"),
                     plotOutput("violinplot"))


### Violin Plot
ViolinSERVER <- function(input, output, session, data){
  
  output$violinplot <- renderPlot({
    
    a <- VlnPlot(data, "nCount_RNA", group.by = input$group) + NoLegend()
    b <- VlnPlot(data, "nFeature_RNA", group.by = input$group) + 
      geom_hline(yintercept = input$features, col = "#832681") + NoLegend()
    c <- VlnPlot(data, "percent.mt", group.by = input$group) + 
      geom_hline(yintercept = input$mt, col = '#FEB078') + NoLegend()
    plot_grid(a, b, c, nrow = 1)
    
  })
  
}