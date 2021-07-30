

#' Mass defect calculations
#'
#' @param mz 
#' @param cus 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
#' 
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
    omd <- mz * round(cus[1]) / cus[1] #rename omd so not to be confused with OMD in MDPlotR.R
    sumd <- cus[2] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      # Here, first order using "round" instead of ceiling as Roach et al. did
      MD1 <- round(round(omd) - omd, digits = 6)
      md2 <- round(round(sumd) - sumd, digits = 6)
      smd <-  MD1 / md2
      MD2 <- round(round(smd) - smd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2)
      
      
    } else if (method == 'floor') {
      MD1 <- round(floor(omd) - omd, digits = 6)
      md2 <- round(floor(sumd) - sumd, digits = 6)
      smd <-  MD1 / md2
      MD2 <- round(floor(smd) - smd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2)
      
    } else {
      MD1 <- round(ceiling(omd) - omd, digits = 6)
      md2 <- round(ceiling(sumd) - sumd, digits = 6)
      smd <-  MD1 / md2
      MD2 <- round(ceiling(smd) - smd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2)
    }
  } else if (length(cus) == 3) {
    omd <- mz * round(cus[1]) / cus[1]
    sumd <- cus[2] * round(cus[1]) / cus[1]
    tumd <- cus[3] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      MD1 <- round(round(omd) - omd, digits = 6)
      md2 <- round(round(sumd) - sumd, digits = 6)
      md3 <- round(round(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(round(smd) - smd, digits = 6)
      md3 <- round(round(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD3 <- round(round(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
      
    } else if (method == 'floor') {
      MD1 <- round(floor(omd) - omd, digits = 6)
      md2 <- round(floor(sumd) - sumd, digits = 6)
      md3 <- round(floor(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(floor(smd) - smd, digits = 6)
      md3 <- round(floor(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD3 <- round(floor(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
      
    } else{
      MD1 <- round(ceiling(omd) - omd, digits = 6)
      md2 <- round(ceiling(sumd) - sumd, digits = 6)
      md3 <- round(ceiling(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(ceiling(smd) - smd, digits = 6)
      md3 <- round(ceiling(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD3 <- round(ceiling(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    }
    
  } else if (length(cus) > 3) {
    message("Sorry, only three MD base units are allowed!")
    omd <- mz * round(cus[1]) / cus[1]
    sumd <- cus[2] * round(cus[1]) / cus[1]
    tumd <- cus[3] * round(cus[1]) / cus[1]
    
    if (method == 'round') {
      MD1 <- round(round(omd) - omd, digits = 6)
      md2 <- round(round(sumd) - sumd, digits = 6)
      md3 <- round(round(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(round(smd) - smd, digits = 6)
      md3 <- round(round(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD3 <- round(round(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
      
    } else if (method == 'floor') {
      MD1 <- round(floor(omd) - omd, digits = 6)
      md2 <- round(floor(sumd) - sumd, digits = 6)
      md3 <- round(floor(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(floor(smd) - smd, digits = 6)
      md3 <- round(floor(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD1_3 <- round(floor(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
      
    } else{
      MD1 <- round(ceiling(omd) - omd, digits = 6)
      md2 <- round(ceiling(sumd) - sumd, digits = 6)
      md3 <- round(ceiling(tumd) - tumd, digits = 6)
      smd <-  MD1 / md2
      tsmd <- md3 / md2
      MD2 <- round(ceiling(smd) - smd, digits = 6)
      md3 <- round(ceiling(tsmd) - tsmd, digits = 6)
      tmd <- MD2 / md3
      MD3 <- round(ceiling(tmd) - tmd, digits = 6)
      re <- cbind.data.frame(mz,MD1,MD2,MD3)
    }
    
  } else{
    
    if (method == 'round') {
      omd <- mz * round(cus) / cus
      MD1 <- round(round(omd) - omd, digits = 6)
      re <- cbind.data.frame(mz,MD1)
      
    } else if (method == 'floor') {
      omd <- mz * floor(cus) / cus
      MD1 <- round(floor(omd) - omd, digits = 6)
      re <- cbind.data.frame(mz,MD1)
      
    } else{
      omd <- mz * ceiling(cus) / cus
      MD1 <- round(ceiling(omd) - omd, digits = 6)
      re <- cbind.data.frame(mz,MD1)
    }
  }
  return(re)
}