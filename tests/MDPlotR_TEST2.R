
library(shiny)
library(shinythemes)
library(dplyr)
library(plotly)
library(crosstalk)
library(DBI)
library(duckdb)

ui <- shiny::navbarPage(
  "MDPlotR: interactive mass defect plots",
  theme = shinythemes::shinytheme('spacelab'),
  shiny::tabPanel("DataAnalysis",
                  shiny::fluidPage(shiny::sidebarLayout(
                    shiny::sidebarPanel(
                      shiny::fileInput(
                        'file1',
                        'Choose Parquet File',
                        accept = c('.parquet')
                      ),
                      shiny::fluidRow(shiny::column(
                        12,
                        shiny::textInput("cus1", "MD formula 1", value = "CH2,O")
                      )),
                      shiny::fluidRow(shiny::column(
                        12,
                        shiny::textInput("cus2", "MD formula 2", value = "Cl-H")
                      )),
                      shiny::actionButton('go', 'Plot', width = "100%"),
                      shiny::tags$br(),
                      shiny::radioButtons(
                        inputId = "rounding",
                        label = "Rounding",
                        choices = c("round", "ceiling", "floor"),
                        selected = "round",
                        inline = TRUE
                      ),
                      shiny::checkboxInput('ins', 'Show intensity as size', FALSE),
                      shiny::checkboxInput("show_leg", "Show plot legends", TRUE),
                      # Plot controls
                      shiny::uiOutput("slide1"),
                      shiny::uiOutput("slide2"),
                      shiny::uiOutput("slide3"),
                      shiny::uiOutput("plotctr"),
                      shiny::uiOutput("plotctr2"),
                      width = 3
                    ),
                    shiny::mainPanel(
                      shiny::uiOutput("plot"),
                      shiny::tags$br(),
                      DT::DTOutput("x1"),
                      plotly::plotlyOutput("barplot"),
                      shiny::fluidRow(shiny::column(
                        3, shiny::downloadButton("x3", "Export Data")
                      )),
                      shiny::tags$br()
                    )
                  ))
  )
)

#---------------- Shiny Server ----------------#

server <- function(input, output, session) {
  
  MD_data <- reactive({
    req(input$file1)
    
    # Connect to DuckDB in-memory
    con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
    
    # Register Parquet file
    dbExecute(con, paste0("CREATE TABLE md AS SELECT * FROM parquet_scan('", input$file1$datapath, "')"))
    
    # Query all data
    df <- dbGetQuery(con, "SELECT * FROM md")
    
    # Mass defect calculations
    df$RMD <- round((round(df$mz) - df$mz) / df$mz * 10^6)
    df$OMD <- round((round(df$mz) - df$mz) * 10^3)
    
    # Higher-order mass defect
    mdh1 <- getmdh(df$mz, cus = input$cus1, method = input$rounding)
    mdh2 <- getmdh(df$mz, cus = input$cus2, method = input$rounding)
    
    name1 <- paste0("Formula1_", colnames(mdh1))
    name2 <- paste0("Formula2_", colnames(mdh2))
    
    mdh <- cbind(mdh1[,-1], mdh2[,-1])
    colnames(mdh) <- c(name1[-1], name2[-1])
    
    df <- cbind(df, mdh)
    
    dbDisconnect(con, shutdown = TRUE)
    
    df
  })
  
  # Filtering sliders
  output$slide1 <- renderUI({
    minZ <- min(MD_data()$intensity)
    maxZ <- max(MD_data()$intensity)
    sliderInput("slide1", "Intensity range filter", min = minZ, max = maxZ, value = c(minZ, maxZ))
  })
  
  output$slide2 <- renderUI({
    minZ <- min(MD_data()$mz)
    maxZ <- max(MD_data()$mz)
    sliderInput("slide2", "m/z range", min = minZ, max = maxZ, value = c(minZ, maxZ))
  })
  
  output$slide3 <- renderUI({
    minZ <- min(MD_data()$rt)
    maxZ <- max(MD_data()$rt)
    sliderInput("slide3", "retention time range", min = minZ, max = maxZ, value = c(minZ, maxZ))
  })
  
  # Plot UI
  output$plot <- renderUI({
    fluidRow(column(6, plotly::plotlyOutput("DTPlot1")),
             column(6, plotly::plotlyOutput("DTPlot2")))
  })
  
  output$plotctr <- renderUI({
    fluidRow(
      h4("Plot controls"),
      tags$br(),
      column(6,
             selectInput(
               inputId = 'selectintensity',
               label = 'Variable for intensity',
               choices = names(dplyr::select(MD_data(), where(is.numeric))),
               selected = "intensity"
             )
      ),
      tags$br(), tags$br(), tags$br(), tags$br(), tags$br(),
      column(6,
             selectInput('xvar1', 'X variable for Plot 1', choices = names(MD_data()), selected = "rt")
      ),
      column(6,
             selectInput('yvar1', 'Y variable for Plot 1', choices = names(MD_data()), selected = "mz")
      ),
      column(6,
             selectInput('xvar2', 'X variable for Plot 2', choices = names(MD_data()), selected = "RMD")
      ),
      column(6,
             selectInput('yvar2', 'Y variable for Plot 2', choices = names(MD_data()), selected = "mz")
      )
    )
  })
  
  output$plotctr2 <- renderUI({
    fluidRow(
      tags$br(),
      textInput('x1', 'x axis label for plot 1', input$xvar1),
      textInput('y1', 'y axis label for plot 1', input$yvar1),
      textInput('x2', 'x axis label for plot 2', input$xvar2),
      textInput('y2', 'y axis label for plot 2', input$yvar2)
    )
  })
  
  # Observe plot button
  observeEvent(input$go, {
    m <- MD_data()
    m <- m[m$intensity >= input$slide1[1] &
             m$intensity <= input$slide1[2] &
             m$mz >= input$slide2[1] &
             m$mz <= input$slide2[2] &
             m$rt >= input$slide3[1] &
             m$rt <= input$slide3[2], ]
    
    d <- crosstalk::SharedData$new(m)
    
    MDplot_y1 <- m[, input$yvar1]
    MDplot_x1 <- m[, input$xvar1]
    
    intensity <- if (input$ins) m[, input$selectintensity] else NULL
    
    MDplot_x2 <- m[, input$xvar2]
    MDplot_y2 <- m[, input$yvar2]
    
    # Plot 1
    output$DTPlot1 <- plotly::renderPlotly({
      plotly::plot_ly(x = MDplot_x1, y = MDplot_y1, type = "scatter", mode = "markers",
                      size = intensity, marker = list(color = intensity, colorscale = "Hot"),
                      showlegend = input$show_leg) %>%
        plotly::layout(xaxis = list(title = input$x1), yaxis = list(title = input$y1))
    })
    
    # Plot 2
    output$DTPlot2 <- plotly::renderPlotly({
      plotly::plot_ly(x = MDplot_x2, y = MDplot_y2, type = "scatter", mode = "markers",
                      size = intensity, marker = list(color = intensity, colorscale = "Hot"),
                      showlegend = input$show_leg) %>%
        plotly::layout(xaxis = list(title = input$x2), yaxis = list(title = input$y2))
    })
    
    # DataTable
    output$x1 <- DT::renderDT({
      DT::datatable(m, editable = TRUE, rownames = FALSE, selection = "none", filter = "top",
                    options = list(scrollX = TRUE))
    })
    
    # Barplot
    nperc <- function(x) round(x / max(x) * 100, 1)
    output$barplot <- plotly::renderPlotly({
      bar_out <- m[d$selection(), ]
      selectInt <- m[d$selection(), input$selectintensity]
      plotly::plot_ly(data = bar_out, x = bar_out$mz, y = nperc(selectInt), type = "bar") %>%
        plotly::layout(xaxis = list(title = "m/z"), yaxis = list(title = "Relative intensity (%)"))
    })
    
    # Export
    output$x3 <- downloadHandler(
      'MDplot_annotated_export.csv',
      content = function(file) {
        write.csv(m[d$selection(), ], file)
      }
    )
  })
  
  # Close app when session ends
  if (!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
}

shiny::shinyApp(ui, server)
