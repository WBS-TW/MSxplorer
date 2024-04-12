 
#' Plot extracted ion chromatograms
#'
#'Plot the EIC from MS1 or MS2 from raw data using MSnbase
#' @param filepath Path to the raw data file(s)
#' @param featlist A data frame containing the columns: name, mz, ms_level
#' @param diff A single value specifying the mass range from the specified m/z values in featlist (mz-diff, mz+diff)
#' @param mserr The mass error in Dalton
#'
#' @return Overlayed extracted chromatograms
#' @export
#'
#' @examples
#' fl <- "D:\\TEST\\28Rterr spiked 1000pg_1_D_4_1.mzML"

#' plotEIC(filepath = filepath, featlist = featlist)

plotHilbertCurve <- function(filepath, featlist, diff = 0.005, mserr = 0.005) {
  data_prof <- MSnbase::readMSData(filepath, mode = "onDisk", centroided = TRUE)
  hd <- MSnbase::fData(data_prof)
  ms1 <- which(hd$msLevel == 1)
  ms2 <- which(hd$msLevel == 2)

  rtselms1 <- hd$retentionTime[ms1] > hd$retentionTime[ms1][1] & hd$retentionTime[ms1] < tail(hd$retentionTime[ms1], n=1)
  rtselms2 <- hd$retentionTime[ms2] > hd$retentionTime[ms2][1] & hd$retentionTime[ms2] < tail(hd$retentionTime[ms2], n=1)
  
  M_plot <- vector("list", nrow(featlist))
  
  p <- plotly::plot_ly() 
  
  for (i in 1:nrow(featlist)) {
    M <- MSnbase::MSmap(data_prof, 
               scans = if(featlist$ms_level[i] == "ms1") {ms1[rtselms1]} 
               else if(featlist$ms_level[i] == "ms2") {ms2[rtselms2]}, 
               lowMz = featlist$mz[i]-diff, 
               highMz = featlist$mz[i]+diff, 
               resMz = mserr, 
               hd = hd, 
               zeroIsNA = FALSE)
    M_plot[[i]] <- data.frame(rt = M@rt, int = M@map) |>
      replace(., is.na(.), 0) |>
      dplyr::mutate(sum_int = rowSums(.[grep("int", names(.))], na.rm = TRUE))
    p <- plotly::add_lines(p, x = M_plot[[i]][["rt"]], y = M_plot[[i]][["sum_int"]], name = paste0(round(featlist$mz[i],4), "_", featlist$ms_level[i])
    ) # do not use $ dollar sign for subsetting list
  }
  p |> plotly::config(showTips = FALSE) 
  #|> plotly::layout(title = grepl(filepath, ))
}
