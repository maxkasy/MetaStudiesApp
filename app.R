# to upload, run rsconnect::deployApp('MetaStudiesApp')
library(shiny)
rm(list = ls())
source("metastudiesfunctions.r")
source("metastudiesplots.r")

ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File of estimates and standard errors",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      hr(),
      h3("Model parameters"),
      fluidRow(
        column(6, checkboxInput("symmetric", "Symmetric p(.)", value = FALSE),
                  checkboxGroupInput("cutoffs", "Cutoffs for p(.)",
                                        choiceNames = c("1.65", "1.96","2.58"),
                                        choiceValues = c(1.645, 1.960, 2.576),
                                  selected = 1.960)),
        column(6, radioButtons("modelmu", "Model for the distribution of effects",
                               choices = c("Normal" = "normal",
                                           "Student-t" = "t" #,
                                              #"nonparametric" = 3
                                           ),
                                           selected = "normal"))
      ),
      #hr(),
      actionButton(inputId = "estimatebutton", label = "Estimate model")
    ),
    
    mainPanel(
      h2("Funnel plot, histogram of z-stats", align = "center"),
      fluidRow(splitLayout(cellWidths = c("50%", "50%"), 
                    plotOutput("funnel"),
                    plotOutput("hist")))
    )
  ),
  hr(),
  h2("Model estimates", align = "center"),
  h4("Distribution of true effects, conditional publication proabilities", align = "center"),
  column(12, align="center",
         tableOutput("estimatestable"),
         plotOutput("estplot", width = "70%")
  )
)



server <- function(input, output, session) {

  #object to store data and estimation results
  v = reactiveValues()

  # read data, generate funnel plot
  output$funnel=renderPlot({
    req(input$file1)
    metadata=read.table(input$file1$datapath,
                        sep=",")
    v$X=metadata[,1]
    v$sigma=metadata[,2]
    metastudies_plot(v$X,v$sigma)
  })
  
  # generate histogram
  output$hist=renderPlot({
    req(v$X)
    z_histogram(v$X, v$sigma)
  })
  
  #estimation
  observeEvent(input$estimatebutton,{
    v$cutoffs=as.numeric(unlist(input$cutoffs))
    v$symmetric=input$symmetric
    v$modelmu=input$modelmu
    if (!v$symmetric) v$cutoffs= c(-rev(v$cutoffs), 0 ,v$cutoffs)
    
    v$estimates=metastudies_estimation(v$X,v$sigma,v$cutoffs, v$symmetric, model= v$modelmu)
    
    output$estplot = renderPlot({
        estimates_plot(v$X,v$sigma,v$cutoffs, v$symmetric,v$estimates, model= v$modelmu)
    })
  })
  
  #render estimates to table
  output$estimatestable =  renderTable(rownames=TRUE,hover=TRUE,digits=3,
   {req(v$estimates)
    estimatestable(v$estimates$Psihat, v$estimates$SE, v$cutoffs, v$symmetric, v$modelmu)
   })


}

# Run the app ----
shinyApp(ui = ui, server = server)