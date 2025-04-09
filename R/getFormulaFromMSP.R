# Parse the chemical formula from msp files as input to HRMF

#' getFormulaFromMSO
#'
#' @param file 
#'
#' @return character vector of chemical formula
#' 
#' 


getFormulaFromMSP <- function(file) {
msp <- readLines(file, warn = FALSE)

# Removing comments that include Name: to avoid incorrect indexing

msp <- stringr::str_replace_all(msp, "Comments:.*(Name:)", "Comment removed due to incorrect name format")

Names <- grep("Name:", msp, ignore.case = TRUE)

msp <- split(msp, cumsum(seq_along(msp) %in% Names))

msp_formula <- character()

for (i in seq_along(msp)) {
  
  m_ind <- msp[[i]]
  
  Formula_index <- stringr::str_detect(m_ind, "FORMULA:")
  Formula <- m_ind[Formula_index]
  Chemical_Formula <- str_extract(Formula, "(?<=:)(.+)") |> str_squish()
  
  msp_formula[i] <- Chemical_Formula
  
}
return(msp_formula)
}
