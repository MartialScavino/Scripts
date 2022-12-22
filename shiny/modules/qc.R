

# source("modules/QC/ComputedThreshold.R")
# source("modules/QC/DoubletFinder.R")
source("modules/QC/FeaturePlot.R")
source("modules/QC/Heatmap.R")
source("modules/QC/Hist.R")
source("modules/QC/Plotly.R")
source("modules/QC/Scatter.R")
source("modules/QC/SideBarPanel.R")
source("modules/QC/Violin.R")
source("modules/QC/Visualisation.R")




headUI <- tags$head(tags$style(type="text/css", '
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
            '))




# FeaturePlotUI <- tabPanel("FeaturePlot",
#                           textInput.typeahead("gene", placeholder = "Veuillez rentrer un gène",
#                                               local = genes, valueKey = "gene",
#                                               tokens =  c(1:length(genes$gene)),
#                                               template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")),
#                           plotOutput("featureplot"))

# FeaturePlotUI <- tabPanel("FeaturePlot",
#                           selectizeInput(
#                             inputId = 'gene', label = 'Select Something',
#                             choices = NULL,
#                             selected = 1,
#                             multiple = T
#                           ),
#                           textOutput("featureplot"))










######################### SERVER ############################







#Feature Plot à partir d'un gène donné
# FeaturePlotSERVER <- function(input, output, session, data){
# 
#   output$featureplot <- renderPlot({
# 
#     if (input$gene == ""){}
# 
#     else if (length(input$gene < 2)){FeaturePlot(data, features = input$gene)}
# 
#      else{
# 
#        data <- AddModuleScore(data, features = input$gene)
#        FeaturePlot(data, features = input$Cluster1)
# 
#        }
# 
#   })
# 
# }


