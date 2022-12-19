
icon_find("rocket")
UI <- fluidPage(
  
  sidebarPanel(
  
                textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
                                    local = genes, 
                                    valueKey = "gene",
                                    tokens =  c(1:length(genes$gene)),
                                    template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
                
                actionButton("reset", "Reinitialiser", "danger"),
          
          p(HTML("<br>")),
          p("Gènes dans la signature (cliquez sur un gène pour le retirer)"),
          p(markdown('---')),
          uiOutput("list_button"),
          p(HTML("<br><br><br>")),
          actionButton('plot_button', 'Calculer la figure',"primary")),
  
  mainPanel(
          
          plotOutput(outputId = "test_features")))



server <- function(input, output, session) {
  
  liste = reactiveValues(gene_name = character(0))
  
    # Ajoute un gène à la liste lorsqu'on écrit un nouveau gène
    observeEvent(input$gene,
      
      if (input$gene != "" & input$gene %in% liste$gene_name == FALSE){
        
        liste$gene_name <- c(liste$gene_name, input$gene)
    })
    
    
    # Créer un bouton pour chaque élément dans la liste de gène signature
    observeEvent(liste$gene_name, output$list_button <- renderUI({
      
      if (length(liste$gene_name) > 0){
      tagList(
       lapply(1:length(liste$gene_name), function(i){

         id <- paste0('slider_',liste$gene_name[i])
         actionButton(id, liste$gene_name[i], "info")
            })
          )
      }
      
      else{
        
        output$list_button <- renderUI(actionButton('bla', "Aucun gène sélectionné", "inverse"))
        
      }
      }))
    
    # Reset list quand on appuie sur le bouton
    observeEvent(input$reset, liste$gene_name <- character(0))
    
    # Make button removable when clicking on it
    observe({
    lapply(1:length(liste$gene_name), function(i){
        
        id <- paste0('slider_',liste$gene_name[i])
        operator <- "$"
        real_id <- do.call(operator, list(input, id))
        observeEvent(real_id, liste$gene_name <- liste$gene_name[-i])
      
        })
      })
    
    # Calcule le plot quand le bouton est sélectionné 
    # et quand le bouton change ou quand on ajoute ou enlève un gène de la liste
    observeEvent(input$plot_button,{
        
        output$test_features <- renderPlot({ 
          isolate({
            if (length(liste$gene_name) == 1){
            
            p <- FeaturePlot(data, features = liste$gene_name[1])
            return(p)
          }
          
            else if (length(liste$gene_name) > 1){
            
            data <- AddModuleScore_UCell(obj = data, features = list(Signature = liste$gene_name))
            p <- FeaturePlot(object = data, features = "Signature_UCell") + ggtitle("Signature")
            return(p)
          }
          
          else { return() }
          
          })
        })
      
        
    })

}       
shinyApp(UI, server)



#################################################

if (interactive()) {
  
  library(shiny)
  library(shinyWidgets)
  
  fruits <- as.vector(genes$gene)
  
  ui <- fluidPage(
    
    
    textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
                        local = genes, 
                        valueKey = "gene",
                        tokens =  c(1:length(genes$gene)),
                        template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
    
    tags$h2("Multi update"),
    multiInput(
      inputId = "my_multi",
      label = "Fruits :",
      choices = character(0),
      width = "350px"
    ),
    verbatimTextOutput(outputId = "res"),
    selectInput(
      inputId = "selected",
      label = "Update selected:",
      choices = character(0),
      multiple = TRUE
    ),
    textInput(inputId = "label", label = "Update label:")
  )
  
  server <- function(input, output, session) {
    
    output$res <- renderPrint(input$my_multi)
    
    observeEvent(input$gene, {
      updateMultiInput(
        session = session,
        inputId = "my_multi",
        selected = input$selected,
        choices = genes$gene[genes$gene %like% input$gene]
      )
    })
    
    observeEvent(input$label, {
      updateMultiInput(
        session = session,
        inputId = "my_multi",
        label = input$label,
        choices = genes$gene[genes$gene %like% input$gene]
      )
    }, ignoreInit = TRUE)
  }
  
  shinyApp(ui, server)
  
}

























library("shiny")
# mylist
genes <- as.vector(rownames(data))
# ui
ui <- fluidPage(
  
  
  searchInput("search", btnSearch = icon("magnifying-glass")),
  selectizeInput(
    inputId = 'mylist', label = 'Select Something',
    choices = NULL,
    selected = 1,
    multiple = T
  )
)


# server
server <- function(input, output, session) {
  observeEvent(input$search, 
               updateSelectizeInput(session = session, 
                                    inputId = 'mylist', 
                                    choices = genes[genes %like% input$search], server = TRUE))
}
# app
shinyApp(ui = ui, server = server)


liste <- (test = character(0))







data <- AddModuleScore(data, features = c("S100A8", "S100A1"), name = "Cluster")
data <- AddModuleScore(data, features = c("S100A13", "S100A8", "S100A1"), name = "Cluster")
head(data)
FeaturePlot(data, "Cluster1")



## Only run examples in interactive R sessions
if (interactive()) {
  
  ui <- fluidPage(
    uiOutput("moreControls")
  )
  
  server <- function(input, output) {
    output$moreControls <- renderUI({
      tagList(
        sliderInput("n", "N", 1, 1000, 500),
        textInput("label", "Label")
      )
    })
  }
  shinyApp(ui, server)
}







library(shiny)
actionButton <- function(inputId, label, btn.style = "" , css.class = "") {
  if ( btn.style %in% c("primary","info","success","warning","danger","inverse","link")) {
    btn.css.class <- paste("btn",btn.style,sep="-")
  } else btn.css.class = ""
  
  tags$button(id=inputId, type="button", class=paste("btn action-button",btn.css.class,css.class,collapse=" "), label)
}

ui <- fluidPage(shinyUI(basicPage(
  actionButton("primary","primary","primary",c("btn-large", "meheh")),
  actionButton("info","info","info"),
  actionButton("success","success","success"),
  actionButton("warning","warning","warning"),
  actionButton("danger","danger","danger"),
  actionButton("inverse","inverse","inverse"),
  actionButton("link","link","link")
)
)
)


server <- function(input,output){
  
}


shinyApp(ui, server)




data <- AddModuleScore_UCell(data, features = list(Signatures = c("S100A8", "S100A13")), name = NULL)
head(data)


if (interactive()) {
  
  ui <- fluidPage(
    actionButton("update", "Update other buttons"),
    br(),
    actionButton("goButton", "Go"),
    br(),
    actionButton("goButton2", "Go 2", icon = icon("area-chart")),
    br(),
    actionButton("goButton3", "Go 3")
  )
  
  server <- function(input, output, session) {
    observe({
      req(input$update)
      
      # Updates goButton's label and icon
      updateActionButton(session, "goButton",
                         label = "New label",
                         icon = icon("calendar"))
      
      # Leaves goButton2's label unchaged and
      # removes its icon
      updateActionButton(session, "goButton2",
                         icon = character(0))
      
      # Leaves goButton3's icon, if it exists,
      # unchaged and changes its label
      updateActionButton(session, "goButton3",
                         label = "New label 3")
    })
  }
  
  shinyApp(ui, server)
}





















ui <- fluidPage(
  useShinyjs(), 
  titlePanel("Old Faithful Geyser Data"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 30)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tags$head(tags$style(type="text/css", '
            .loading {
                display: inline-block;
                overflow: hidden;
                height: 1.3em;
                margin-top: -0.3em;
                line-height: 1.5em;
                vertical-align: text-bottom;
                box-sizing: border-box;
            }
            .loading.dots::after {
                text-rendering: geometricPrecision;
                content: "⠋\\A⠙\\A⠹\\A⠸\\A⠼\\A⠴\\A⠦\\A⠧\\A⠇\\A⠏";
                animation: spin10 1s steps(10) infinite;
                animation-duration: 1s;
                animation-timing-function: steps(10);
                animation-delay: 0s;
                animation-iteration-count: infinite;
                animation-direction: normal;
                animation-fill-mode: none;
                animation-play-state: running;
                animation-name: spin10;
            }
            .loading::after {
                display: inline-table;
                white-space: pre;
                text-align: left;
            }
            @keyframes spin10 { to { transform: translateY(-15.0em); } }
            ')),
      plotOutput("distPlot"),
      actionButton("btnUpdate", span("Update", id = "UpdateAnimate", class = "loading dots"))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  UpdatePlot <- reactiveVal(1)
  
  output$distPlot <- renderPlot({
    req(UpdatePlot())
    
    Sys.sleep(1) # just for show
    # Button settings        
    shinyjs::enable("btnUpdate")
    shinyjs::removeClass(id = "UpdateAnimate", class = "loading dots")
    
    x    <- faithful[, 2]
    isolate( # update chart ONLY after button press
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
    )
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })
  
  
  
  observeEvent(input$btnUpdate, { # User requests update
    UpdatePlot(NULL)
    
    shinyjs::addClass(id = "UpdateAnimate", class = "loading dots")
    shinyjs::disable("btnUpdate")
    
    UpdatePlot(TRUE)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


















if (interactive()) {
  ui <- fluidPage(
    tags$h1("Search Input"),
    br(),
    searchInput(
      inputId = "search", label = "Enter your text",
      placeholder = "A placeholder",
      btnSearch = icon("magnifying-glass"),
      btnReset = icon("xmark"),
      width = "450px"
    ),
    br(),
    verbatimTextOutput(outputId = "res")
  )
  
  server <- function(input, output, session) {
    output$res <- renderPrint({
      input$search
    })
  }
  
  shinyApp(ui = ui, server = server)
}
