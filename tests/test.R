compound[1,1] <- 215.085526
compound[1,2] <- 100
compound[2,1] <- 216.0888
compound[2,2] <- 19
compound <- compound[1:2, 1:2]

formula <- "C18H14"



########
compound <- data.frame(mz = c(263.77794, 264.78134, 265.77589, 266.77929, 267.77384, 268.77724, 269.77180, 270.77520), 
                       intensity = c(34.26, 0.75, 100, 2.2, 97.29, 2.14, 31.55, 0.69))
formula <- "C2H3Br3"


#####Compare with HRF from Tracefinder for PAHs
compound <- read.csv("./data/acenaphthylene_deconvolutedTraceFinder.csv") %>%
  select(mz, intensity)
formula <- "C12H8"


###



# the read_msp was sourced. Change this when finalizing the function
# Need to manually input formula or extract from msp file? Since unknowns are not annotated with chemical formula and not shown in msp

#msp <- "./data/2-Methoxy-5-methylaniline.msp"
#formula <- "C8H11NO"
# 
# msp <- "./data/MTM00088_GC-EI-FT_POSITIVE_YATIGPZCMOYEGE-UHFFFAOYSA-N.msp"
# formula <- "C14H8Br6O2"

#msp <- "./data/MTM00126_GC-EI-FT_POSITIVE_WYEHFWKAOXOVJD-UHFFFAOYSA-N.msp"
#formula <- "C19H11F5N2O2"

library(dplyr)
library(enviPat)
library(stringr)
library(tidyr)
library(rcdk)
#from Florian dust sample MSDIAL
msp <- "./data/benzothiazole_dust.msp"
formula <- "C7H5NS" 

HRMF_output <- HRMF(msp = msp, formula = formula)
HRMF_scores <- HRMF_output$HRMF_total

# msp <- "./data/CI Pigment yellow 151_QCxMS.msp"
# formula <-"C18H15N5O5" 


####

msp <- "./data/Spectrum_0_202111221149_NIST_cut.msp"
formula <- "C7H5NS"

HRMF_output <- HRMF(msp = msp, formula = formula)
