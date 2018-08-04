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
        column(6, checkboxInput("symmetric", "Symmetric p(.)", value = TRUE),
                  checkboxGroupInput("cutoffs", "Cutoffs for p(.)",
                                        choiceNames = c("1.65", "1.96","2.33"),
                                        choiceValues = c(1.645, 1.960, 2.326),
                                  selected = 1.960)),
        column(6, radioButtons("modelmu", "Model for the distribution of effects",
                               choices = c("Normal" = 1 #,
                                              #"Student-t" = 2,
                                              #"nonparametric" = 3
                                           ),
                                           selected = 1))
      ),
      #hr(),
      actionButton(inputId = "estimatebutton", label = "Estimate model")
    ),
    
    mainPanel(
      h2("Funnel plot, histogram", align = "center"),
      fluidRow(splitLayout(cellWidths = c("50%", "50%"), 
                    plotOutput("funnel"),
                    plotOutput("hist")))
    )
  ),
  hr(),
  h2("Model estimates", align = "center"),
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
    if (!v$symmetric) v$cutoffs= c(-rev(v$cutoffs), 0 ,v$cutoffs)
    
    v$estimates=metastudies_estimation(v$X,v$sigma,v$cutoffs, v$symmetric)
    
    output$estplot = renderPlot({
        estimates_plot(v$X,v$sigma,v$cutoffs, v$symmetric,v$estimates)
    })
  })
  
  #render estimates to table
  output$estimatestable =  renderTable(rownames=TRUE,hover=TRUE,digits=3,
   {req(v$estimates)
    l=length(v$estimates$Psihat)
    estimates=matrix(0,2,l)
    estimates[1,]=v$estimates$Psihat
    estimates[2,]=v$estimates$SE
    rownames(estimates)=c("estimate", "standard error")
    colnames(estimates)=rep(" ",l)
    colnames(estimates)[1]=intToUtf8(956) #mu
    colnames(estimates)[2]=intToUtf8(964) #tau
    if (v$symmetric){
      colnames(estimates)[3]=paste("[0,", v$cutoffs[1],"]")
      for (i in seq(2, length(v$cutoffs), length=max(0, length(v$cutoffs) - 1))) {
        colnames(estimates)[2+i]=paste("(", v$cutoffs[i-1], ",", v$cutoffs[i],"]")
      }
    } else {
      colnames(estimates)[3]=paste("(-", intToUtf8(8734), ",", v$cutoffs[1],"]")
      for (i in seq(2, length(v$cutoffs), length=max(0, length(v$cutoffs) - 1))) {
        colnames(estimates)[2+i]=paste("(", v$cutoffs[i-1], ",", v$cutoffs[i],"]")
      }
    }
    #colnames(estimates)[3]=expression(beta_p)
    estimates
   })


}

# Run the app ----
shinyApp(ui = ui, server = server)