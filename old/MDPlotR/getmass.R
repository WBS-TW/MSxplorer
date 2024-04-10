


#' getmass
#' Get the exact mass from chemical formula using rcdk
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
#' 
getmass <- function(data) {
  if (grepl('-', data)) {
    name <- unlist(strsplit(data, '-'))
    iso1 <- rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
    iso2 <- rcdk::get.isotopes.pattern(rcdk::get.formula(name[2]))
    cus <- as.numeric(iso1[max(iso1[, 2]), 1]) - as.numeric(iso2[max(iso2[, 2]), 1])
    
    } else if (grepl("/", data)) {
      name <- unlist(strsplit(data, "/"))
      frac <- as.numeric(name[2])
      iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
      cus <- as.numeric(iso[max(iso[, 2]), 1]) / frac
      
    }else{
    iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(data))
    cus <- as.numeric(iso[max(iso[, 2]), 1])
  }
  return(cus)
}