#' Extract most intense MS1 ions from MS2 EIC using Shiny 
#'
#' @param filepath string. Path to mzML file
#' @param featlist data.frame. Variable names should be: name, mz, ms_level (has to be either "ms1, "ms2) 
#' @param numTopIons num number of most intensive MS1 ions to extract
#' @param diff num range of mass (Da) to extract from the mz specified in featlist (mz-diff, mz+diff)
#' @param mserr num mass error (Da) of the range boundaries specified in diff (mz-diff+-mserr, mz+dif+-mserr)
#' @param rtWindow num the retention time range to extract MS1 EIC (in sec). Error occurs if set too low.
#'
#' @return overlaid EIC from MS2 and MS1 in an interactive Shiny environment
#' @export
#'
#' @examples
#' fl <- "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/B6 batch std_1_F,2_1.mzML"
#' 
#' PFSA_frags <- data.frame(name = c("FSO3", "SO3"), mz = c(98.9552, 79.9558), ms_level = c("ms2", "ms2"))
#' 
#' plotTopMS1Peaks(filepath = fl, featlist = PFSA_frags, numTopIons = 3)
#' 
#' fl <- "D:/Raw_data/Dust_Florian/LC/LC_neg_mzMLQC_rep1.mzML"
#' Br_frags <- data.frame(name = c("Br79", "Br81"), mz = c(78.9183, 80.9163), ms_level = c("ms2", "ms2"))
#' plotTopMS1Peaks(filepath = fl, featlist = Br_frags, numTopIons = 3)


plotTopMS1Peaks <- function(filepath, featlist, numTopIons = 10, diff = 0.01, mserr = 0.01, rtWindow = 0.3) {
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Click on the peaks ions to select retention time and then click Get MS1"),
    miniUI::miniContentPanel(
      shiny::fillRow(flex = c(NA,1),
                     shiny::fillCol(width = "110px",
                                    shiny::textOutput("rtselect"),
                                    shiny::textOutput("rtWindow"),
                                    shiny::actionButton("sync1", "Sync 1<-2")),
                     shiny::fillCol(flex = c(1, 1),
                                    plotly::plotlyOutput("plot1" ,height = "100%"), 
                                    plotly::plotlyOutput("plot2", height = "100%")
                     )
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("getMS1", "Get MS1"),
      shiny::actionButton("exp_excel", "Export selection to Excel"),
      shiny::actionButton("exp_metfrag", "Copy selection in Metfrag format")
    )
  )
  
#----------SERVER--------------------------------------------

  server <- function(input, output, session) {
    
    data_prof <- MSnbase::readMSData(filepath, mode = "onDisk", centroided = TRUE)
    hd <- MSnbase::fData(data_prof)
    ms1 <- which(hd$msLevel == 1)
    ms2 <- which(hd$msLevel == 2)
    
    rtselms1 <- hd$retentionTime[ms1] > hd$retentionTime[ms1][1] & hd$retentionTime[ms1] < tail(hd$retentionTime[ms1], n=1)
    rtselms2 <- hd$retentionTime[ms2] > hd$retentionTime[ms2][1] & hd$retentionTime[ms2] < tail(hd$retentionTime[ms2], n=1)
    
    M_plot <- vector("list", nrow(featlist))
    rtr <- NULL
    rtselection1 <- NULL

    p1 <- plotly::plot_ly() 
    p2 <- plotly::plot_ly() 
    
    output$plot1 <- plotly::renderPlotly({
      
      
      for (i in 1:nrow(featlist)) {
        M <- MSnbase::MSmap(data_prof, 
                            scans = if(featlist$ms_level[i] == "ms1") {ms1[rtselms1]} 
                            else if(featlist$ms_level[i] == "ms2") {ms2[rtselms2]}, 
                            lowMz = featlist$mz[i]-diff, 
                            highMz = featlist$mz[i]+diff, 
                            resMz = mserr, 
                            hd = hd, 
                            zeroIsNA = FALSE)
        M_plot[[i]] <- data.frame(rt = M@rt, int = M@map) %>%
          replace(., is.na(.), 0) %>%
          dplyr::mutate(sum_int = rowSums(.[grep("int", names(.))], na.rm = TRUE))
        p1 <- plotly::add_lines(p1, 
                                x = M_plot[[i]][["rt"]], 
                                y = M_plot[[i]][["sum_int"]], 
                                name = paste0(round(featlist$mz[i],4), "_", featlist$ms_level[i])
        ) # do not use $ dollar sign for subsetting list
      }
      p1 %>% plotly::config(showTips = FALSE) %>% plotly::layout(showlegend = TRUE)
    })  
    
    
    
#### Show the clicked point to be input in next plot
    shiny::observeEvent(plotly::event_data("plotly_click"), {
      rtr <- as.data.frame(plotly::event_data("plotly_click"))[[3]]
      rtr <- c(round((rtr - rtWindow), 1), round((rtr+ rtWindow), 1))
      rtr <<- rtr
      
      output$rtselect <- renderText(paste0("RT:", rtr[1], "-", rtr[2]))
      
    })
    
#### Plot new EIC of MS1 from the clicked point, when pressing "Get MS1"
    shiny::observeEvent(input$getMS1, {
      MS1 <- data_prof %>%
        MSnbase::filterRt(rtr) %>% #rtr is from taken from the global environment rtr <<-
        MSnbase::filterMsLevel(1)
      
      
      
      MS1_spec <- data.frame(mz = MS1[[1]]@mz, intensity = MS1[[1]]@intensity) %>%
        dplyr::arrange(desc(intensity)) %>%
        dplyr::slice_head(n = numTopIons)
      
      featlist2 <- MS1_spec %>%
        dplyr::mutate(name = paste0("S", 1:numTopIons), mz = MS1_spec$mz, ms_level = "ms1") %>%
        dplyr::select(name, mz, ms_level, -intensity)
      
      
      featlist2 <- featlist2 %>% dplyr::bind_rows(featlist)
      
      output$plot2 <- plotly::renderPlotly({
        for (i in 1:nrow(featlist2)) {
          M <- MSnbase::MSmap(data_prof, 
                              scans = if(featlist2$ms_level[i] == "ms1") {ms1[rtselms1]} 
                              else if(featlist2$ms_level[i] == "ms2") {ms2[rtselms2]}, 
                              lowMz = featlist2$mz[i]-diff, 
                              highMz = featlist2$mz[i]+diff, 
                              resMz = mserr, 
                              hd = hd, 
                              zeroIsNA = FALSE)
          M_plot[[i]] <- data.frame(rt = M@rt, int = M@map) %>%
            replace(., is.na(.), 0) %>%
            dplyr::mutate(sum_int = rowSums(.[grep("int", names(.))], na.rm = TRUE))
          p2 <- plotly::add_lines(p2, 
                                  x = M_plot[[i]][["rt"]], 
                                  y = M_plot[[i]][["sum_int"]], 
                                  name = paste0(round(featlist2$mz[i],4), "_", featlist2$ms_level[i])
          ) # do not use $ dollar sign for subsetting list
        }
        p2 %>% plotly::config(showTips = FALSE) %>% plotly::layout(showlegend = TRUE)
    
      })  
      
    })

    
#### Show selected rt range
    #### Sync ranges of plotly1 <- plotly2 
    shiny::observeEvent(plotly::event_data("plotly_relayout"), {
      rtselection1 <- plotly::event_data("plotly_relayout")
      rtselection1 <- c(rtselection1[1], rtselection1[2])
      output$rtWindow <- shiny::renderText(paste0("RT range:", "\n", rtselection1))
    })
    
    
#### Sync ranges of plotly1 <- plotly2 
    shiny::observeEvent(input$sync1, {
      rtselection1 <- plotly::event_data("plotly_relayout")
      rtselection1 <- c(rtselection1[1], rtselection1[2])
      print(rtselection1)
      output$rtWindow <- shiny::renderText(paste0(rtselection1))
      
      output$plot1 <- plotly::renderPlotly({
        
        
        for (i in 1:nrow(featlist)) {
          M <- MSnbase::MSmap(data_prof, 
                              scans = if(featlist$ms_level[i] == "ms1") {ms1[rtselms1]} 
                              else if(featlist$ms_level[i] == "ms2") {ms2[rtselms2]}, 
                              lowMz = featlist$mz[i]-diff, 
                              highMz = featlist$mz[i]+diff, 
                              resMz = mserr, 
                              hd = hd, 
                              zeroIsNA = FALSE)
          M_plot[[i]] <- data.frame(rt = M@rt, int = M@map) %>%
            replace(., is.na(.), 0) %>%
            dplyr::mutate(sum_int = rowSums(.[grep("int", names(.))], na.rm = TRUE))
          p1 <- plotly::add_lines(p1, 
                                  x = M_plot[[i]][["rt"]], 
                                  y = M_plot[[i]][["sum_int"]], 
                                  name = paste0(round(featlist$mz[i],4), "_", featlist$ms_level[i])
          ) # do not use $ dollar sign for subsetting list
        }
        p1 %>% 
          plotly::config(showTips = FALSE) %>%
          plotly::layout(xaxis = list(range = c(rtselection1[[1]], rtselection1[[2]])),
                         showlegend = TRUE)
          
      })  
      
    })

    
    
    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
  }
  
  #runGadget(ui, server, viewer = dialogViewer("MS2 EIC", width = 700, height = 700))
  shiny::runGadget(ui, server, viewer = shiny::browserViewer())
}

