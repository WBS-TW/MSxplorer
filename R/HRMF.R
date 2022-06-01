# Actually, using Sirius can help predict and annotate all peaks.

# TODO
# 1. first Generate the monoisotopic mass -> HRMF monoiso score
# 2. then add isotope patter from isopattern -> HRMF full iso score
# 3. add score on isotopic match RMS -> % iso match score
# 4. Combined score 

# 5. loop through top n candidates? (how to get this? from NIST? or some external file generated from NIST search)
# 6. loop through all names in msp file. How to get annotated formula?


# use rcdk and isopattern
# https://pubs.acs.org/doi/full/10.1021/acs.analchem.5b01503
# https://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/

# example 2-Methoxy-5-methylaniline, C8H11NO

# read_msp
library(devtools)
library(dplyr)
library(enviPat)
library(stringr)
source("./R/read_msp.R")
data("isotopes")

# state the mass accuracy in ppm
mass_accuracy <- 10

compound <- read_msp("./data/2-Methoxy-5-methylaniline.msp")[[1]] %>%
  select(mz, intensity)
rownames(compound) <- NULL

formula <- "C8H11NO"

atoms <- as.data.frame(rcdk::get.formula(formula)@isotopes) %>%
  mutate(number = as.numeric(number)) %>%
  mutate(mass = as.numeric(mass))

element_range <- list()
for (i in 1:nrow(atoms)) {
  element_range[[i]] <- c(atoms[i,][[1]], 0, atoms[i,][[2]])
  }

# the windows, in Da, state the +- window from the accurate mass based on the setting on mass accuracy
windows <- mass_accuracy/10^6*compound$mz

# initiate an empty list
match_comp <- list()


for (i in seq_along(compound$mz)) {
  
match <- rcdk::generate.formula(
  compound$mz[[i]],
  window = windows[[i]],
  elements = element_range,
  validation = FALSE,
  charge = 1
)

match_list <- list(found = logical(), annodf = data.frame(), isopat = data.frame())

if(length(match) > 0) {
match_list[["found"]] <- TRUE

annodf <- data.frame(Mass = match[[1]]@mass,
                     Formula = match[[1]]@string,
                     MassError_ppm = round((compound$mz[i]-match[[1]]@mass)/match[[1]]@mass*10^6, 1))
match_list[["annodf"]] <- annodf

chemforms <- str_extract(match_list[["annodf"]][["Formula"]], "(?<=\\[).+?(?=\\])")

isopat <- as.data.frame(isopattern(isotopes = isotopes, chemforms = chemforms,
                                   threshold = 0.001, 
                                   charge = 1, 
                                   emass = 0.00054858, 
                                   plotit = FALSE, 
                                   algo=1, 
                                   rel_to = 0,
                                   verbose = TRUE)) %>%
  rename(Mass = 1,
         Abundance = 2)

match_list[["isopat"]] <- isopat


} else {
  match_list[["found"]] <- FALSE
  match_list[["annodf"]] <- NULL
  match_list[["isopat"]] <- NULL
}

match_comp[[i]] <- match_list
}






use#' High resolution mass filtering
#'
#' @param mass 
#'
#' @return
#' @export
#'
#' @examples
HRMF <- function(mass) {
  

  
}
