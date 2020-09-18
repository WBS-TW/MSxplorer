
# plot function to extract EIC from MS1 or MS2
plotEIC <- function(filepath, featlist, err = 0.005, mserr = 0.005) {
  data_prof <- readMSData(filepath, mode = "onDisk", centroided = TRUE)
  hd <- fData(data_prof)
  ms1 <- which(hd$msLevel == 1)
  ms2 <- which(hd$msLevel == 2)
  err <- err
  mserr <- mserr
  
  rtselms1 <- hd$retentionTime[ms1]
  rtselms2 <- hd$retentionTime[ms2]
  
  M_plot <- vector("list", nrow(featlist))
  
  p <- plot_ly() 
  
  for (i in 1:nrow(featlist)) {
    M <- MSmap(data_prof, 
               scans = if(featlist$ms_level[i] == "ms1") {ms1[rtselms1]} 
               else if(featlist$ms_level[i] == "ms2") {ms2[rtselms2]}, 
               lowMz = featlist$mz[i]-err, 
               highMz = featlist$mz[i]+err, 
               resMz = mserr, 
               hd = hd, 
               zeroIsNA = FALSE)
    M_plot[[i]] <- data.frame(rt = M@rt, int = M@map) %>%
      replace(., is.na(.), 0) %>%
      mutate(sum_int = rowSums(.[grep("int", names(.))], na.rm = TRUE))
    p <- add_lines(p, x = M_plot[[i]][["rt"]], y = M_plot[[i]][["sum_int"]], name = paste0(featlist$mz[i], "_", featlist$ms_level[i])
    ) # do not use $ dollarsign for subsetting list
  }
  p %>% config(showTips = FALSE)
}
