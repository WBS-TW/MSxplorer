#' getmass - helper function to MDPlotR
#' Get the exact mass from chemical formula using rcdk
#' @param data 
#'
#' @return vector
#' @export 
getmass <- function(data) {
  data("isotopes")
  if (grepl('-', data)) {
    name <- unlist(strsplit(data, '-'))
    iso1 <- as.double(enviPat::isopattern(chemforms = name[1], isotopes = isotopes)[[1]][[1,1]])
    iso2 <- as.double(enviPat::isopattern(chemforms = name[2], isotopes = isotopes)[[1]][[1,1]])
    cus <- iso1 - iso2
    
  } else if (grepl("/", data)) {
    name <- unlist(strsplit(data, "/"))
    frac <- as.double(name[2])
    iso <- as.double(enviPat::isopattern(chemforms = name[1], isotopes = isotopes)[[1]][[1,1]])
    cus <- iso / frac
    
  }else{
    cus <- as.double(enviPat::isopattern(chemforms = data, isotopes = isotopes)[[1]][[1,1]])
  }
  return(cus)
}
