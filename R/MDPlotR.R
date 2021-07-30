#' MDPlotR
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
MDPlotR <- function() {

ui <- shiny::navbarPage(
  "MDPlotR: interactive mass defect plots",
  theme = shinythemes::shinytheme('spacelab'),
  shiny::tabPanel("DataAnalysis",
          shiny::fluidPage(shiny::sidebarLayout(
             shiny::sidebarPanel(
              shiny::fileInput(
                 'file1',
                 'Choose CSV File',
                 accept = c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv')
               ),
               shiny::fluidRow(shiny::column(
                 12,
                 shiny::textInput("cus1", "MD formula 1", value = 'CH2,O')
               )
               
               ),
               
              shiny::fluidRow(shiny::column(
                 12,
                 shiny::textInput("cus2", "MD formula 2", value = 'Cl-H')
               )),
              shiny::actionButton('go', 'Plot', width = '100%'),
              shiny::tags$br(),
              shiny::radioButtons(
                inputId = "rounding",
                label = "Rounding",
                choices = c("round", "ceiling", "floor"),
                selected = "round",
                inline = TRUE
              ),
              shiny::checkboxInput('ins', 'Show intensity as size', F),
              shiny::checkboxInput("show_leg", "Show plot legends", T),
              #Plot controls
              shiny::uiOutput("plotctr"),
              shiny::uiOutput("plotctr2"),
              shiny::uiOutput("slide1"),
              shiny::uiOutput("slide2"),
              shiny::uiOutput("slide3"),
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
               shiny::tags$br(),
               )
             )
             )),
  shiny::tabPanel(
    "Instructions",
    shiny::sidebarLayout(
      shiny::sidebarPanel(shiny::h3("Table of content"),
                          shiny::h4("File input"),
                          shiny::h4("Equation"),
                   width = 3),
      shiny::mainPanel(
        shiny::withMathJax(shiny::includeMarkdown("data/MDPlotR_instructions.md"))
      )
    )
  )
)



#-----------------------Shiny Server function----------#



server = function(input, output, session) {
  MD_data <- reactive({
    req(input$file1) #requires that the input is available
    df <- vroom::vroom(input$file1$datapath) # use vroom for faster loading of large files
    df$RMD <- round((round(df$mz) - df$mz) / df$mz * 10 ^ 6)
    df$OMD <- round((round(df$mz) - df$mz) * 10 ^ 3)
    
    # higher-order mass defect calculation based on: doi.org/10.1021/ac200654j
    
    mdh1 <- getmdh(df$mz,cus = input$cus1, method = input$rounding)
    mdh2 <- getmdh(df$mz,cus = input$cus2, method = input$rounding)
    name1 <- paste0("Formula1_", colnames(mdh1))
    name2 <- paste0("Formula2_", colnames(mdh2))
    mdh <- cbind(mdh1[,-1],mdh2[,-1])
    colnames(mdh) <- c(name1[-1],name2[-1])
    
    df <- cbind(df,mdh)
    
  })
  # Filtering the intensity, mz, and rt
  output$slide1 <- renderUI({
    minZ <- min(MD_data()$intensity)
    maxZ <- max(MD_data()$intensity)
    
    sliderInput(
      "slide1",
      "Intensity range filter",
      min = minZ,
      max = maxZ,
      value = c(minZ, maxZ)
    )
  })
  output$slide2 <- shiny::renderUI({
    minZ <- min(MD_data()$mz)
    maxZ <- max(MD_data()$mz)
    
    sliderInput(
      "slide2",
      "mass to charge ratio range",
      min = minZ,
      max = maxZ,
      value = c(minZ, maxZ)
    )
  })
  output$slide3 <- shiny::renderUI({
    minZ <- min(MD_data()$rt)
    maxZ <- max(MD_data()$rt)
    
    shiny::sliderInput(
      "slide3",
      "retention time range",
      min = minZ,
      max = maxZ,
      value = c(minZ, maxZ)
    )
  })
  
  
  ## for plot control ##
  output$plot <- shiny::renderUI({
    shiny::fluidRow(shiny::column(6, plotly::plotlyOutput("DTPlot1")),
                      shiny::column(6, plotly::plotlyOutput("DTPlot2")))
    })
  
  output$plotctr <- shiny::renderUI({
        shiny::fluidRow(
        shiny::h4("Plot controls"),
        shiny::tags$br(),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'selectintensity', #select which variable to use as intensity input for bar plot (not working yet for scatterplot?)
            label = 'Variable for intensity',
            choices = names(dplyr::select(MD_data(), where(is.numeric))), # select only numeric variables
            selected = names(MD_data()["intensity"])
          )
        ),
        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'xvar1',
            label = 'X variable for Plot 1',
            choices = names(MD_data()),
            selected = names(MD_data()["rt"])
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'yvar1',
            label = 'Y variable for Plot 1',
            choices = names(MD_data()),
            selected = names(MD_data()["mz"])
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'xvar2',
            label = 'X variable for Plot 2',
            choices = names(MD_data()),
            selected = names(MD_data()["RMD"])
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'yvar2',
            label = 'Y variable for Plot 2',
            choices = names(MD_data()),
            selected = names(MD_data()["mz"])
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'zvar1',
            label = 'Symbol variable for plot',
            choices = list(`NULL` = 'NA',`Variable` = names(MD_data())),
            selected = 'NULL'
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'zvar2',
            label = 'Symbol variable for plot 2',
            choices = list(`NULL` = 'NA',`Variable` = names(MD_data())),
            selected = 'NULL'
          )
        )
      )
    })
  output$plotctr2 <- shiny::renderUI({
        shiny::fluidRow(
        shiny::tags$br(),
        shiny::textInput(
          'x1',
          'x axis label for plot 1',
          input$xvar1
        ),
        shiny::textInput(
          'y1',
          'y axis label for plot 1',
          input$yvar1
        ),
        shiny::textInput(
          'x2',
          'x axis label for plot 2',
          input$xvar2
        ),
        shiny::textInput(
          'y2',
          'y axis label for plot 2',
          input$yvar2
        )
      )
    })
 
#------For MD Plot Panel-----
  
  #OE#
  shiny::observeEvent(input$go, {
    m <- MD_data()
    m <-
      m[m$intensity >= input$slide1[1] &
          m$intensity <= input$slide1[2] &
          m$mz >= input$slide2[1] &
          m$mz <= input$slide2[2] &
          m$rt >= input$slide3[1] &
          m$rt <= input$slide3[2],]
    
    d <- crosstalk::SharedData$new(m)
    
    
    
    MDplot_y1 <- m[, input$yvar1]
    
    MDplot_x1 <- m[, input$xvar1]
    
    if(input$zvar1 == 'NA'){
      MDplot_z1 <- 1
    }else{
      MDplot_z1 <- m[, input$zvar1]
    }
    
    # Checkbox option for size of markers by intensity variable
    if (input$ins) {
      intensity <- m[, input$selectintensity]
    } else{
      intensity <- NULL
    }
    
    MDplot_x2 <-
        m[, input$xvar2]
      
      MDplot_y2 <-
        m[, input$yvar2]
      
      if(input$zvar2 == 'NA'){
        MDplot_z2 <- 1
      }else{
        MDplot_z2 <- m[, input$zvar2]
      }
      

    
#-----Plot 1-------
    output$DTPlot1 <- plotly::renderPlotly({
      s <- input$x1_rows_selected
      if (!length(s)) {
        p <- d %>%
          plotly::plot_ly(
            x = MDplot_x1,
            y = MDplot_y1,
            symbol = MDplot_z1,
            showlegend = input$show_leg,
            type = "scatter",
            size = intensity,
            mode = "markers",
            marker = list(
              line = list(
                width = 1,
                color = '#FFFFFF'
              )
            ),
            color = I('black'),
            #color = ~intensity,
            name = 'Unfiltered'
          ) %>%
          plotly::layout(
            legend = list(
              orientation = "h",
              xanchor = "center",
              x = 0.5,
              y = 100
            ),
            showlegend = T,
            xaxis = list(title = input$x1),
            yaxis = list(title = input$y1)
          ) %>%
          plotly::highlight(
            "plotly_selected",
            color = I('red'),
            selected = plotly::attrs_selected(name = 'Filtered')
          )
      } else if (length(s)) {
        pp <- m %>%
          plotly::plot_ly() %>%
          plotly::add_trace(
            x = MDplot_x1,
            y = MDplot_y1,
            symbol = MDplot_z1,
            type = "scatter",
            size = intensity,
            mode = "markers",
            marker = list(
              line = list(
                width = 1,
                color = '#FFFFFF'
              )
            ),
            color = I('black'),
            #color = ~intensity,
            name = 'Unfiltered'
          ) %>%
          plotly::layout(
            legend = list(
              orientation = "h",
              xanchor = "center",
              x = 0.5,
              y = 100
            ),
            showlegend = T,
            xaxis = list(title = input$x1),
            yaxis = list(title = input$y1)
          )
        
        # selected data
        pp <-
          plotly::add_trace(
            pp,
            data = m[s, , drop = F],
            x = MDplot_x1[s],
            y = MDplot_y1[s],
            type = "scatter",
            size = intensity[s],
            mode = "markers",
            marker = list(
              line = list(
                width = 1,
                color = '#FFFFFF'
              )
            ),
            color = I('red'),
            name = 'Filtered'
          )
      }
      

      
    })
    
#-----Plot 2-------
    output$DTPlot2 <- plotly::renderPlotly({
        t <- input$x1_rows_selected
        
        if (!length(t)) {
          p <- d %>%
            plotly::plot_ly(
              x = MDplot_x2,
              y = MDplot_y2,
              symbol = MDplot_z2,
              showlegend = input$show_leg,
              type = "scatter",
              size = intensity,
              mode = "markers",
              marker = list(
                line = list(
                  width = 1,
                  color = '#FFFFFF'
                )
              ),
              color = I('black'),
              #color = ~intensity,
              name = 'Unfiltered'
            ) %>%
            plotly::layout(
              legend = list(
                orientation = "h",
                xanchor = "center",
                x = 0.5,
                y = 100
              ),
              showlegend = T,
              xaxis = list(title = input$x2),
              yaxis = list(title = input$y2)
            ) %>%
            plotly::highlight(
              "plotly_selected",
              color = I('red'),
              selected = plotly::attrs_selected(name = 'Filtered')
            )
        } else if (length(t)) {
          pp <- m %>%
            plotly::plot_ly() %>%
            plotly::add_trace(
              x = MDplot_x2,
              y = MDplot_y2,
              symbol = MDplot_z2,
              type = "scatter",
              size = intensity,
              mode = "markers",
              marker = list(
                line = list(
                  width = 1,
                  color = '#FFFFFF'
                )
              ),
              color = I('black'),
              #color = ~intensity,
              name = 'Unfiltered'
            ) %>%
            plotly::layout(
              legend = list(
                orientation = "h",
                xanchor = "center",
                x = 0.5,
                y = 100
              ),
              showlegend = T,
              xaxis = list(title = input$x2),
              yaxis = list(title = input$y2)
            )
          
          # selected data
          pp <-
            plotly::add_trace(
              pp,
              data = m[t, , drop = F],
              x = MDplot_x2[t],
              y = MDplot_y2[t],
              type = "scatter",
              size = intensity[t],
              mode = "markers",
              marker = list(
                line = list(
                  width = 1,
                  color = '#FFFFFF'
                )
              ),
              color = I('red'),
              name = 'Filtered'
            )
        }
        
      })
    

#----- Datatable of selected rows--------
    
    ## TO FIX: if highlighted twice then the DT does not index correctly
    
    output$x1 <- DT::renderDT({
      T_out1 <- DT::datatable(
        m[d$selection(), ],
        editable = TRUE,
        rownames = FALSE,
        filter = "top",
        options = list(
          scrollX = TRUE)
      )
      
      dt <-
        DT::datatable(
          m,
          editable = TRUE,
          rownames = FALSE,
          filter = "top",
          options = list(
            scrollX = TRUE)
        )
      
      if (NROW(T_out1) == 0) {
        dt
      } else {
        T_out1
      }
    })
    
#-----Barplot------
    
    # Function to normalize to the highest intensity of selected features
    nperc <- function(x) {
      norm <- round(x/max(x) * 100, 1)
      return(norm)
    }
    # generate the barplot only when selecting data
    output$barplot <- plotly::renderPlotly({
      
      bar_out <- m[d$selection(), ]
      selectInt <- m[d$selection(), input$selectintensity] # intensity based on the selectintensity selectinput
      
      
      pbar <-  plotly::plot_ly() %>%
          plotly::add_trace(
            data = bar_out,
            x = bar_out$mz,
            y = nperc(selectInt), 
            type = "bar") %>%
        plotly::layout(
          xaxis = list(title = "m/z"),
          yaxis = list(title = "Relative intensity (%)")
        )
      })
    
#-----Exporting the annotated data------
    
    ##NOT FULLY WORKING YET?###
    
    output$x3 <- shiny::downloadHandler(
      'MDplot_annotated_export.csv',
      content = function(file) {
        s <- input$x1_rows_selected
        if (length(s)) {
          write.csv(m[s, , drop = FALSE], file)
        } else if (!length(s)) {
          write.csv(m[d$selection(), ], file)
          
          
        }
      }
    )
  })
  
  # Close the app when the session completes
  if(!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
  
}

shiny::shinyApp(ui, server)
}
