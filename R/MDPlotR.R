


#------ShinyApp------------------------------

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
#require(shiny, dplyr, vroom, shinythemes, DT, plotly, crosstalk, rcdk)

ui <- shiny::navbarPage(
  "MDPlotR: interactive mass defect plots",
  theme = shinythemes::shinytheme('spacelab'),
  shiny::tabPanel("Get Started",
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
              shiny::br(),
              shiny::radioButtons(
                 inputId = "single",
                 label = "Single or Double plots",
                 choices = c("Single", "Double"),
                 selected = "Single",
                 inline = TRUE
               ),
              shiny::checkboxInput('ins', 'Show intensity as size', F),
              shiny::checkboxInput("show_leg", "Show plot legends", T),
              shiny::uiOutput("plotctr"),
              shiny::uiOutput("plotctr2"),
              shiny::uiOutput("slide1"),
              shiny::uiOutput("slide2"),
              shiny::uiOutput("slide3"),
               width = 3
             ),
             shiny::mainPanel(
               shiny::uiOutput("plot"),
               DT::DTOutput("x1"),
               shiny::fluidRow(shiny::column(
                 3, shiny::downloadButton("x3", "Export Data")
               )),
               shiny::tags$br(),
               )
             )
             ))
  # shiny::tabPanel(
  #   "Instructions",
  #   shiny::sidebarLayout(
  #     shiny::sidebarPanel(shiny::h3("Table of content"),
  #                         shiny::h4("File input"),
  #                         shiny::h4("Equation"),
  #                  width = 3),
  #     shiny::mainPanel(
  #       shiny::withMathJax(shiny::includeMarkdown("instructions.md"))
  #     )
  #   )
  # )
)



#-----------------------Shiny Server functon----------#


#--FUNCTIONS---------------------------

# Initiate functions parsed from EnviGCMS #
getmass <- function(data) {
  if (grepl('-', data)) {
    name <- unlist(strsplit(data, '-'))
    iso1 <- rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
    iso2 <- rcdk::get.isotopes.pattern(rcdk::get.formula(name[2]))
    cus <- as.numeric(iso1[max(iso1[, 2]), 1]) - as.numeric(iso2[max(iso2[, 2]), 1])
  } else{
    iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(data))
    cus <- as.numeric(iso[max(iso[, 2]), 1])
  }
  return(cus)
}

getmdh <- function(mz, cus = c("CH2,H2"), method = "round"){
  getorder <- function(input) {
    if (grepl(',', input)) {
      name <- unlist(strsplit(input, ','))
    } else{
      name <- input
    }
    return(name)
  }
  temp <- getorder(cus)
  cus <- NULL
  for (i in 1:length(temp)) {
    cus <- c(cus, getmass(temp[i]))
  }
  if (length(cus) == 2) {
    omd <- mz * round(cus[1]) / cus[1]
    sumd <- cus[2] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      
      MD1 <-
        round(round(omd) - omd,
              digits = 6)
      md2 <- round(round(sumd) - sumd,
                   digits = 6)
      smd <-  MD1 / md2
      MD2 <-
        round(round(smd) - smd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2)
      
    } else if (method == 'floor') {
      MD1 <-
        round(floor(omd) - omd,
              digits = 6)
      md2 <- round(floor(sumd) - sumd,
                   digits = 6)
      smd <-  MD1 / md2
      MD2 <-
        round(floor(smd) - smd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2)
      
    } else if (method == 'ceiling'){
      MD1 <-
        round(ceiling(signif(omd)) - omd,
              digits = 6)
      md2 <- round(ceiling(signif(sumd))- sumd,
                   digits = 6)
      smd <-  MD1 / md2
      MD2 <-
        round(ceiling(signif(smd)) - smd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD1)
    }
  } else if (length(cus) == 3) {
    omd <- mz * round(cus[1]) / cus[1]
    sumd <- cus[2] * round(cus[1]) / cus[1]
    tumd <- cus[3] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      MD1 <-
        round(round(omd) - omd,
              digits = 6)
      md2 <- round(round(sumd) - sumd,
                   digits = 6)
      md3 <- round(round(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(round(smd) - smd,
              digits = 6)
      md3 <- round(round(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD3 <-
        round(round(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    } else if (method == 'floor') {
      MD1 <-
        round(floor(omd) - omd,
              digits = 6)
      md2 <- round(floor(sumd) - sumd,
                   digits = 6)
      md3 <- round(floor(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(floor(smd) - smd,
              digits = 6)
      md3 <- round(floor(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD3 <-
        round(floor(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    } else{
      MD1 <-
        round(ceiling(omd) - omd,
              digits = 6)
      md2 <- round(ceiling(sumd) - sumd,
                   digits = 6)
      md3 <- round(ceiling(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(ceiling(smd) - smd,
              digits = 6)
      md3 <- round(ceiling(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD3 <-
        round(ceiling(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    }
    
  } else if (length(cus) > 3) {
    message("Sorry, only the first three unit would be used.")
    omd <- mz * round(cus[1]) / cus[1]
    sumd <- cus[2] * round(cus[1]) / cus[1]
    tumd <- cus[3] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      MD1 <-
        round(round(omd) - omd,
              digits = 6)
      md2 <- round(round(sumd) - sumd,
                   digits = 6)
      md3 <- round(round(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(round(smd) - smd,
              digits = 6)
      md3 <- round(round(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD3 <-
        round(round(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    } else if (method == 'floor') {
      MD1 <-
        round(floor(omd) - omd,
              digits = 6)
      md2 <- round(floor(sumd) - sumd,
                   digits = 6)
      md3 <- round(floor(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(floor(smd) - smd,
              digits = 6)
      md3 <- round(floor(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD1_3 <-
        round(floor(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    } else{
      MD1 <-
        round(ceiling(omd) - omd,
              digits = 6)
      md2 <- round(ceiling(sumd) - sumd,
                   digits = 6)
      md3 <- round(ceiling(tumd) - tumd,
                   digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <-
        round(ceiling(smd) - smd,
              digits = 6)
      md3 <- round(ceiling(tsmd) - tsmd,
                   digits = 6)
      tmd <- MD2 / md3
      MD3 <-
        round(ceiling(tmd) - tmd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    }
  } else{
    
    if (method == 'round') {
      omd <- mz * round(cus) / cus
      MD1 <-
        round(round(omd) - omd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1)
    } else if (method == 'floor') {
      omd <- mz * floor(cus) / cus
      MD1 <-
        round(floor(omd) - omd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1)
    } else{
      omd <- mz * ceiling(cus) / cus
      MD1 <-
        round(ceiling(omd) - omd,
              digits = 6)
      re <- cbind.data.frame(mz,MD1)
    }
  }
  return(re)
}



server = function(input, output, session) {
  MD_data <- reactive({
    #  require that the input is available
    req(input$file1)
    df <- vroom::vroom(input$file1$datapath)
    df$RMD <- round((round(df$mz) - df$mz) / df$mz * 10 ^ 6)
    df$OMD <- round((round(df$mz) - df$mz) * 10 ^ 3)
    # high order mass defect computation
    
    mdh1 <- getmdh(df$mz,cus = input$cus1)
    mdh2 <- getmdh(df$mz,cus = input$cus2)
    # change colname
    name1 <- paste0(colnames(mdh1),'_p1')
    name2 <- paste0(colnames(mdh2),'_p2')
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
    if (input$single == "Single") {
      plotly::plotlyOutput("DTPlot1")
    } else{
      shiny::fluidRow(shiny::column(6, plotly::plotlyOutput("DTPlot1")),
                      shiny::column(6, plotly::plotlyOutput("DTPlot2")))
    }
  })
  output$plotctr <- shiny::renderUI({
    if (input$single == "Single") {
      shiny::fluidRow(
        shiny::h4("Plot controls"),
        shiny::tags$br(),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'xvar1',
            label = 'X variable for plot',
            choices = names(MD_data())
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'yvar1',
            label = 'Y variable for plot',
            choices = names(MD_data())
          )
        ),
        shiny::column(
          12,
          shiny::selectInput(
            inputId = 'zvar1',
            label = 'Symbol variable for plot',
            choices = list(`NULL` = 'NA',`Variable` = names(MD_data()))
            ,
            selected = 'NULL'
          )
        )
      )
      
    } else{
      shiny::fluidRow(
        shiny::h4("Plot controls"),
        shiny::tags$br(),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'xvar1',
            label = 'X variable for Plot 1',
            choices = names(MD_data()),
            selected = names(MD_data())[1]
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'yvar1',
            label = 'Y variable for Plot 1',
            choices = names(MD_data()),
            selected = names(MD_data())[4]
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'xvar2',
            label = 'X variable for Plot 2',
            choices = names(MD_data()),
            selected = names(MD_data())[1]
          )
        ),
        shiny::column(
          6,
          shiny::selectInput(
            inputId = 'yvar2',
            label = 'Y variable for Plot 2',
            choices = names(MD_data()),
            selected = names(MD_data())[4]
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
    }
  })
  output$plotctr2 <- shiny::renderUI({
    if (input$single == "Single") {
      shiny::fluidRow(
        shiny::tags$br(),
        shiny::textInput('x1', 'x axis label', input$xvar1),
        shiny::textInput('y1', 'y axis label', input$yvar1)
      )
    } else{
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
    }
    
  })
  #### For MD Plot Panel ####
  
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
    
    MDplot_y1 <-
      m[, input$yvar1]
    
    MDplot_x1 <-
      m[, input$xvar1]
    
    if(input$zvar1 == 'NA'){
      MDplot_z1 <- 1
    }else{
      MDplot_z1 <- m[, input$zvar1]
    }
    
    # Checkbox option for size of markers by intensity
    if (input$ins) {
      intensity <- m$intensity
    } else{
      intensity <- NULL
    }
    
    if (input$single == "Double") {
      MDplot_x2 <-
        m[, input$xvar2]
      
      MDplot_y2 <-
        m[, input$yvar2]
      
      if(input$zvar2 == 'NA'){
        MDplot_z2 <- 1
      }else{
        MDplot_z2 <- m[, input$zvar2]
      }
      
    }
    
    # highlight selected rows in the scatterplot
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
    
    # Plot 2
    if (input$single == "Double") {
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
    }
    # highlight selected rows in the table
    output$x1 <- DT::renderDT({
      T_out1 <- m[d$selection(), ]
      dt <-
        DT::datatable(
          m,
          editable = TRUE,
          rownames = FALSE,
          filter = "top"
        )
      if (NROW(T_out1) == 0) {
        dt
      } else {
        T_out1
      }
    })
    
    # exporting the annotated data
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
