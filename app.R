# to upload, run rsconnect::deployApp('MetaStudiesApp')

library(shiny)
rm(list = ls())
source("metastudiesfunctions.r")

ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  sidebarLayout(
    sidebarPanel(
      width=5,
      fileInput("file1", "Choose CSV File of estimates and standard errors",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      hr(),
      h3("Model parameters"),
      fluidRow(
        column(6, checkboxInput("symmetric", "Symmetric p(.)", value = TRUE),
                  checkboxGroupInput("cutoffs", "Cutoffs for p(.)",
                                  choices=c("1.65" = 1.645,
                                            "1.96" = 1.960,
                                            "2.33" = 2.326),
                                  selected = 1.960)),
        column(6, radioButtons("modelmu", "Model for the distribution of effects",
                               choices = c("Normal" = 1,
                                              "Student-t" = 2,
                                              "nonparametric" = 3),selected = 1))
      ),
      hr(),
      actionButton(inputId = "estimatebutton", label = "Estimate model")
    ),
    
    mainPanel(
      width=7,
      # tableOutput("estimatestable"),
      hr(),
      plotOutput("funnel", width = "100%")
    )
  )
)



server <- function(input, output, session) {

  #object to store estimation results
  v = reactiveValues()

  
  observeEvent(input$estimatebutton,{
    req(input$file1)
    
    #loading covariate file
    metadata=read.table(input$file1$datapath,
                          sep=",")
    X=metadata[,1]
    sigma=metadata[,2]
    
    # generate funnel plot
    output$funnel=renderPlot({
      metastudies_plot(X,sigma)
    })
    
    #invoke estimation commands here
    #v$estimates=metastudies_estimation(X,sigma, input$cutoffs,input$symmetric)
    
  })
  
  
  # render estimates to table
  
  # output$estimatestable =  renderTable({
  #   v$estimates
  #  })


}

# Run the app ----
shinyApp(ui = ui, server = server)