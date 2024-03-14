
# NOT WORKING YET

#' Export to NIST msp format
#' Code modified from Stettin et al. doi: 10.3390/metabo10040143
#'
#' @param data an XSAnnotate object generated from CAMERA
#' @param min.clustersize numeric how many fragments need to be in a group for it to be included in the results. 5 - trace compounds; 20 - mostly high quality spectra
#' @param num.files numeric (length = 1) specify how many files are included in the analysis
#' @param out.file character name of the output msp file (.msp will be automatically added)
#' @param alkaneSeries 
#'
#' @return
#' @export
#'
#' @examples
#' alkaneSeries <- data.frame(Cnum = c(11, 13, 15, 17, 19, 21, 23, 25), rt = c(4.90, 8.00, 11.13, 14.05, 16.70, 19.37, 22.45, 25.77))
#' 
exportMSP <- function(data, out.file, min.clustersize = 15, num.files, alkaneSeries = NULL) {
  
  #the following script identifies all "compounds" considered too small according to min.clustersize
  peaktable <- getPeaklist(data, intval="maxo")
  pcgroups <- sort(unique(as.numeric(peaktable$pcgroup)))
  lspectra <- NULL
  extra <- NULL
  small.clusters <- NULL
  big.clusters <- NULL
  n <- 1
  for(x in pcgroups) {
    clustersize <- length(peaktable[peaktable$pcgroup==x,ncol(peaktable)])
    if(clustersize < as.numeric(min.clustersize))
    {
      small.clusters <- c(small.clusters, x)
    }
    else
    {
      big.clusters <- c(big.clusters, x)
      group1 <- peaktable[peaktable$pcgroup==x,]
      group1 <- group1[, -(2:(ncol(group1)-3-num.files))]
      group1 <- group1[, -(ncol(group1):(ncol(group1)-2))]
      decider <- NULL
      for(i in 2:ncol(group1))
      {
        decider <- c(decider, sum(group1[, i]))
      }
      highest <- which(decider==max(decider), arr.ind = TRUE)
      group1 <- data.frame(group1[, 1], group1[, highest[1]+1])
      colnames(group1) <- c("mz", "int")
      lspectra[[n]] <- group1
      n <- n+1
    }
  }
  
  # remove peaks with zero intensities
  lspectra <- lapply(lspectra, function(x) {x[x$int != 0,]})
  
  
  reduced.peaktable <- getReducedPeaklist(xs8, method = "sum", intval = "into", default.adduct.info = "maxint", mzrt.range = FALSE, npeaks.sum = FALSE, cleanup = FALSE)
  if(is.null(small.clusters)==FALSE) {
    for(z in small.clusters)
    {
      reduced.peaktable <- reduced.peaktable[-(which(reduced.peaktable[, ncol(reduced.peaktable)]==z)), ]
    }
  }
  
  #write.csv(reduced.peaktable, file = "peaktable_MatchedFilter.csv") #This creates a peaktable .csv file with only "compounds" instead of m/z features
  extra <- data.frame(Name = paste("Unknown", big.clusters, "RT =",round(reduced.peaktable$rt/60, 2) ), Class = "Unknown", RT = round(reduced.peaktable$rt/60, 2),  stringsAsFactors = FALSE)
  
  ####KOVATS#####
  #Adding Kovats retention index to the extra object to write to msp
  RTI <- data.frame(rt = extra$RT)

  
  RI <- vector(length = nrow(RTI))
  for (i in seq_len(nrow(RTI))) {
    m <- dplyr::case_when(
      RTI$rt[i] >= alkaneSeries$rt[1] & RTI$rt[i] < alkaneSeries$rt[2] ~ alkaneSeries$Cnum[1],
      RTI$rt[i] >= alkaneSeries$rt[2] & RTI$rt[i] < alkaneSeries$rt[3] ~ alkaneSeries$Cnum[2],
      RTI$rt[i] >= alkaneSeries$rt[3] & RTI$rt[i] < alkaneSeries$rt[4] ~ alkaneSeries$Cnum[3],
      RTI$rt[i] >= alkaneSeries$rt[4] & RTI$rt[i] < alkaneSeries$rt[5] ~ alkaneSeries$Cnum[4],
      RTI$rt[i] >= alkaneSeries$rt[5] & RTI$rt[i] < alkaneSeries$rt[6] ~ alkaneSeries$Cnum[5],
      RTI$rt[i] >= alkaneSeries$rt[6] & RTI$rt[i] < alkaneSeries$rt[7] ~ alkaneSeries$Cnum[6],
      RTI$rt[i] >= alkaneSeries$rt[7] & RTI$rt[i] <= alkaneSeries$rt[8] ~ alkaneSeries$Cnum[7]
    )
    
    n <- dplyr::case_when(
      RTI$rt[i] >= alkaneSeries$rt[1] & RTI$rt[i] < alkaneSeries$rt[2] ~ alkaneSeries$Cnum[2],
      RTI$rt[i] >= alkaneSeries$rt[2] & RTI$rt[i] < alkaneSeries$rt[3] ~ alkaneSeries$Cnum[3],
      RTI$rt[i] >= alkaneSeries$rt[3] & RTI$rt[i] < alkaneSeries$rt[4] ~ alkaneSeries$Cnum[4],
      RTI$rt[i] >= alkaneSeries$rt[4] & RTI$rt[i] < alkaneSeries$rt[5] ~ alkaneSeries$Cnum[5],
      RTI$rt[i] >= alkaneSeries$rt[5] & RTI$rt[i] < alkaneSeries$rt[6] ~ alkaneSeries$Cnum[6],
      RTI$rt[i] >= alkaneSeries$rt[6] & RTI$rt[i] < alkaneSeries$rt[7] ~ alkaneSeries$Cnum[7],
      RTI$rt[i] >= alkaneSeries$rt[7] & RTI$rt[i] <= alkaneSeries$rt[8] ~ alkaneSeries$Cnum[8]
    )
    
    
    RI[i] <- round(100*n + 100*(m-n) * (RTI$rt[i] - alkaneSeries[alkaneSeries$Cnum == n,]$rt)/(alkaneSeries[alkaneSeries$Cnum == m,]$rt - alkaneSeries[alkaneSeries$Cnum == n,]$rt), 0)
  }
  
  extra <- cbind(extra, RI)
  
  ######END KOVATS######
  
  # create the msp object with spectra and all meta information in the extra object (NAME, RETENTIONTIME, RETENTIONINDEX)
  export.msp <- construct.msp(lspectra, extra)
  
  #This creates a NIST MS Search compatible .msp file with all compound pseudospectra
  write.msp(export.msp, file = paste0(out.file, ".msp")) 
  
  
}