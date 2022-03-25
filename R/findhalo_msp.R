library(shiny)
library(shinydashboard)
library(rhandsontable)
library(tidyverse)
library(plotly)


#####functions#######

#----findhalo_rt------
findhalo_rt <- function(mz,
                        intensity,
                        sf = 79/78.917789,
                        step = 0.001,
                        stepsd1=0.003,
                        stepsd2=0.005,
                        mzc=700,
                        cutoffint = 1000,
                        cutoffr=0.4,
                        rt,
                        clustercf =5){
  mzr <- round(mz)
  sm <- mz*sf
  sd <- ceiling(sm)-sm
  smsd <- ifelse(mz<=mzc,stepsd1,stepsd2)
  smstep <- seq(0,1,step)
  rt <- rt
  
  data <- cbind.data.frame(mz=mz,
                           mzr = mzr,
                           sm = sm,
                           sd =sd,
                           intensity=intensity,
                           rt = rt)
  data2 <<- data
  
  result <- NULL
  for(i in 1:length(smstep)){
    maxi = smstep[i]+smsd
    mini = smstep[i]-smsd
    index = sd<maxi & sd>mini
    
    li <- data[index & intensity > cutoffint,]
    mzt <- mzr[index & intensity > cutoffint]
    rtt <- rt[index & intensity > cutoffint]
    
    if(length(mzt)>=2){
      c <- cutree(hclust(dist(mzt)),h=clustercf)
      t <- cutree(hclust(dist(rtt)), h = clustercf)
      u <- paste0(c, t)
      cn <- length(unique(u))
      lit <- cbind.data.frame(li,u,i)
      for (j in 1:cn){
        li2 <- lit[lit[,7]==j,]
        mzt2 <- lit$mzr[lit[,7]==j]
        
        if(length(mzt2)>=2){
          if(length(unique(li2$intensity))>1){
            ratio <- max(li2$intensity[li2$intensity != max(li2$intensity)]) / max(li2$intensity)
            diff <- abs(li2$mzr[round(li2$intensity) == round(max(li2$intensity[li2$intensity != max(li2$intensity)]))] - li2$mzr[which.max(li2$intensity)])
          }else{
            ratio <- 1
            diff <- abs(li2$mzr[1]-li2$mzr[2])
          }
          
          if(ratio>cutoffr&round(diff)==2){
            result <- rbind.data.frame(result,li2)
          }
        }
      }
    }
  }
  return(result[!duplicated(result$mz), ])
}

#----findhalo_no_rt---------
findhalo_no_rt <- function(mz,
                           intensity,
                           sf = 79/78.917789,
                           step = 0.001,
                           stepsd1=0.003,
                           stepsd2=0.005,
                           mzc=700,cutoffint = 1000,
                           cutoffr=0.4,
                           clustercf = 5){
  mzr <- round(mz)
  sm <- mz*sf
  sd <- ceiling(sm)-sm
  smsd <- ifelse(mz<=mzc,stepsd1,stepsd2)
  smstep <- seq(0,1,step)
  
  data <- cbind.data.frame(mz=mz, mzr = mzr, sm = sm,sd =sd, intensity=intensity)
  data2 <<- data
  
  result <- NULL
  for(i in 1:length(smstep)){
    maxi = smstep[i]+smsd
    mini = smstep[i]-smsd
    index = sd<maxi & sd>mini
    
    li <- data[index&intensity>cutoffint,]
    mzt <- mzr[index&intensity>cutoffint]
    
    if(length(mzt)>=2){
      c <- cutree(hclust(dist(mzt)),h= clustercf)
      cn <- length(unique(c))
      lit <- cbind.data.frame(li,c,i)
      for (j in 1:cn){
        li2 <- lit[lit[,6]==j,]
        mzt2 <- lit$mzr[lit[,6]==j]
        
        if(length(mzt2)>=2){
          if(length(unique(li2$intensity))>1){
            ratio <- max(li2$intensity[li2$intensity != max(li2$intensity)]) / max(li2$intensity)
            diff <- abs(li2$mzr[round(li2$intensity) == round(max(li2$intensity[li2$intensity != max(li2$intensity)]))] - li2$mzr[which.max(li2$intensity)])
          }else{
            ratio <- 1
            diff <- abs(li2$mzr[1]-li2$mzr[2])
          }
          
          if(ratio>cutoffr&round(diff)==2){
            result <- rbind.data.frame(result,li2)
          }
        }
      }
    }
  }
  return(result[!duplicated(result$mz), ])
}

#-----read_msp-------

read_msp <- function(msp_file) {
  
  # read msp file
  msp <- readLines(msp_file, warn = FALSE)
  msp <- stringr::str_replace_all(msp, "Comments:.*(Name:)", "Comment removed due to incorrect name format")
  
  n <- grep("Name:", msp, ignore.case = TRUE)
  
  msp <- split(msp, cumsum(seq_along(msp) %in% n))
  
  peak_tbl_list <- list()
  
  for (i in seq_along(msp)) {
    
    m_ind <- msp[[i]]
    
    numpeaks <- grep("Num Peaks", m_ind, ignore.case = TRUE)
    peaks <- m_ind[numpeaks+1:length(m_ind)]
    peaks <- na.omit(peaks)
    peak_list <- as.list(peaks)
    
    peak_tbl <- data.frame(mz = NA, intensity = NA, annotation = NA)
    
    maxlength <- 3L
    
    # check that some peak tables are separated by tab or space
    if(grepl("\t", peak_list[[1]]) == TRUE){
      for (i in seq_along(peak_list)) {
        i <- unlist(strsplit(peak_list[[i]], "\t"))
        i <- c(i, rep(NA, maxlength-length(i)))
        peak_tbl <- rbind(peak_tbl, i)
      }
    }else if (grepl("\t", peak_list[[1]]) == FALSE) {
      for (i in seq_along(peak_list)) {
        i <- unlist(strsplit(peak_list[[i]], " "))
        i <- c(i, rep(NA, maxlength-length(i)))
        peak_tbl <- rbind(peak_tbl, i)
      }
    }
    
    find_allNA <- which(apply(peak_tbl, 1, function(x)all(is.na(x))))
    
    peak_tbl <- peak_tbl[-find_allNA,]
    peak_tbl$mz <- as.numeric(peak_tbl$mz)
    peak_tbl$intensity <- as.numeric(peak_tbl$intensity)
    
    peak_tbl_list <- append(peak_tbl_list, list(peak_tbl))
    
    
  }
  
  
  return(peak_tbl_list)
  
}



#------UI-------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = "Findhalo: msp file"),
  dashboardSidebar(width = 300,
                   sidebarMenu(
                     fileInput("file1", "Select msp file"),
                     checkboxInput("check_rt", "Include retention times"),
                     numericInput("msp_id", "Which msp id to read", value = 1, min = 1),
                     actionButton("click_file", "Select"),
                     fluidRow(
                       column(4, checkboxInput("check_int", "Show intensity")),
                       column(4, checkboxInput("check_overlay", "Overlay data"))
                     ),
                     numericInput("md_nom", "Nomin mass", value = 79, width = "60%"),
                     numericInput("md_exact", "Exact mass", value = 78.917789),
                     sliderInput("Y_step", "Y_step", min = 0.001, max = 0.1, value = 0.001),
                     numericInput("X_change_sd", "X_change_sd", value = 700),
                     sliderInput("Y1_sd", "Y1_sd", min = 0.001, max = 0.1, value = 0.003),
                     sliderInput("Y2_sd", "Y2_sd", min = 0.001, max = 0.1, value = 0.005),
                     numericInput("threshold_diff", "Threshold difference", value = 100/25),
                     numericInput("threshold_int_min", "Minimum intensity threshold", value = 1000),
                     sliderInput("cluster_cf", "Cluster height", min = 1, max = 20, step = 1, value = 5)
                   )
  ),
  dashboardBody(
    fluidRow(
      box(plotlyOutput("plot1"),  width = 12)),
    fluidRow(
      box(rHandsontableOutput("table1"),  width = 12))
  )
)


#-------------SERVER--------------------------------------------------------------

server <- function(input, output, session){
  raw_data <- reactive({
    req(input$file1)
    df <- read_msp(input$file1$datapath)[[input$msp_id]]
    return(df)
  })
  
  
  result <- reactive({
    if(input$check_rt){
      raw_data <- raw_data()
      t <- findhalo_rt(mz = raw_data$mz, 
                       intensity = raw_data$intensity, 
                       sf = input$md_nom/input$md_exact, 
                       step = input$Y_step,stepsd1=input$Y1_sd, 
                       stepsd2=input$Y2_sd,
                       mzc=input$X_change_sd,
                       cutoffint = input$threshold_int_min, 
                       cutoffr=1/input$threshold_diff, 
                       rt = raw_data$rt, 
                       clustercf = input$cluster_cf
      )
      return(t)
      
    }else{
      raw_data <- raw_data()
      t <- findhalo_no_rt(mz = raw_data$mz,
                          intensity = raw_data$intensity,
                          sf = input$md_nom/input$md_exact,
                          step = input$Y_step,stepsd1=input$Y1_sd, 
                          stepsd2=input$Y2_sd,
                          mzc=input$X_change_sd,
                          cutoffint = input$threshold_int_min,
                          cutoffr=1/input$threshold_diff,
                          clustercf = input$cluster_cf)
      return(t)
    }
  })
  
  
  # OBSERVEEVENT #
  observeEvent(input$click_file, {
    
    
    
    
    ## PLOTTING OUTPUT##
    
    ###Datatable###
    output$table1 <- renderRHandsontable({
      rhandsontable(result(), width = 700, height = 400) %>%
        hot_cols(fixedColumnsLeft = 1) %>%
        hot_rows(fixedRowsTop = 1)
      
    })
    
    ###Plot###
    output$plot1 <- renderPlotly({
      
      
      if(input$check_rt){
        df <- result()
        
        if (input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i,
              hovertext = paste("RT: ", df$rt, "\n",
                                "mz: ", df$mz)
            )
        } else {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              color = df$i,
              hovertext = paste("RT: ", df$rt, "\n",
                                "mz: ", df$mz))
        }
        
        if (input$check_overlay) {
          p <- p %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, hovertext = paste("RT: ", df$rt, "\n",
                                                                            "mz: ", df$mz))
        } else {
          p <- p
        }
        
        if (input$check_overlay & input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i,
              hovertext = paste("RT: ", df$rt, "\n",
                                "mz: ", df$mz)
            ) %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, size = ~data2$intensity, hovertext = paste("RT: ", df$rt, "\n",
                                                                                                     "mz: ", df$mz))
        } else {
          p <- p
        }
        
        
      }else{
        
        df <- result()
        
        if (input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i)
        } else {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              color = df$i)
        }
        
        if (input$check_overlay) {
          p <- p %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2)
        } else {
          p <- p
        }
        
        if (input$check_overlay & input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i
            ) %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, size = ~data2$intensity)
        } else {
          p <- p
        }
        
        
      }
      
      
      
      
    })
    
  }) # END OF OBSERVEEVENT #
  
  
}

#------------RUN APP-----------------------------------------------

shinyApp(ui, server)

# testing performance
# profvis::profvis(
#   runApp(shinyApp(ui, server))
# )