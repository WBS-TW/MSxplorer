#' Extract most intense MS1 ions from MS2 EIC 
#'
#' @param filepath string path to mzML file
#' @param flagfragments data.frame with variables: name, mz, ms_level ("ms1, "ms2) 
#' @param numTopIons num number of most intensive MS1 ions to extract
#' @param diff num range of mass (Da) to extract from the mz specified in flagfragments (mz-diff, mz+diff)
#' @param mserr num mass error (Da) of the range boundaries specified in diff (mz-diff+-mserr, mz+dif+-mserr)
#' @param rtrange num the retention time range to extract MS1 EIC (in sec). Error occurs if set too low.
#'
#' @return overlaid EIC from MS2 and MS1
#' @export
#'
#' @examples
#' fl <- "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/B6 batch std_1_F,2_1.mzML"
#' PFSA_frags <- data.frame(name = c("FSO3", "SO3"), mz = c(98.9552, 79.9558), ms_level = c("ms2", "ms2"))
#' plotTopMS1Peaks(filepath = fl, flagfragments = PFSA_frags, numTopIons = 3)
plotTopMS1Peaks <- function(filepath, flagfragments, numTopIons = 10, diff = 0.01, mserr = 0.01, rtrange = 0.3) {
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Click on the peaks ions to select retention time and then click Get MS1"),
    miniUI::miniContentPanel(
      shiny::fillRow(flex = c(NA,1),
              shiny::fillCol(width = "100px",
                      shiny::textOutput("rtselect")),
              shiny::fillCol(flex = c(1),
                      plotly::plotlyOutput("plot1"), 
                      plotly::plotlyOutput("plot2")
                      )
              )
      ),
    miniUI::miniButtonBlock(
      shiny::actionButton("getMS1", "Get MS1"),
      shiny::actionButton("exp_excel", "Export selection to Excel"),
      shiny::actionButton("exp_metfrag", "Copy selection in Metfrag format")
      )
    )
  
  server <- function(input, output, session) {
    
    data_prof <- MSnbase::readMSData(filepath, mode = "onDisk", centroided = TRUE)
    
    output$plot1 <- plotly::renderPlotly({
      MSXploreR::plotEIC(filepath = filepath, featlist = flagfragments, diff = diff, mserr = mserr)
    })  
    
    rtr <- NULL

   # Output the clicked point to be input in next plot
    shiny::observeEvent(plotly::event_data("plotly_click"), {
      rtr <- as.data.frame(plotly::event_data("plotly_click"))[[3]]
      rtr <- c(round((rtr - rtrange), 1), round((rtr+ rtrange), 1))
      rtr <<- rtr
      
      output$rtselect <- renderText(paste0("RT:", rtr[1], "-", rtr[2]))
      
    })
    

    # Plot new EIC of MS1 from the clicked point, when pressing "Get MS1"
    shiny::observeEvent(input$getMS1, {
      MS1 <- data_prof %>%
        MSnbase::filterRt(rtr) %>%
        MSnbase::filterMsLevel(1)
      
    
      
      MS1_spec <- data.frame(mz = MS1[[1]]@mz, intensity = MS1[[1]]@intensity) %>%
        dplyr::arrange(desc(intensity)) %>%
        dplyr::slice_head(n = numTopIons)

      featlist <- MS1_spec %>%
        dplyr::mutate(name = paste0("S", 1:numTopIons), mz = MS1_spec$mz, ms_level = "ms1") %>%
        dplyr::select(name, mz, ms_level, -intensity)


      featlist <- featlist %>% dplyr::bind_rows(flagfragments)
      
      output$plot2 <- plotly::renderPlotly({
        MSXploreR::plotEIC(filepath = fl, featlist = featlist, diff = diff, mserr = diff)
      })
      
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
    }
  
  #runGadget(ui, server, viewer = dialogViewer("MS2 EIC", width = 700, height = 700))
  shiny::runGadget(ui, server, viewer = shiny::browserViewer())
}

