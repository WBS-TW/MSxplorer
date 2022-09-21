

# test using GC orbitrap 2,6.TDI
file <- "D:\\R_projects\\MTM_HRMS_LIB\\3_MSP\\MTM00013_GC-EI-FT_POSITIVE_RUELTTOHQODFPA-UHFFFAOYSA-N.msp"

formula <- as.character(readxl::read_xlsx("D:\\R_projects\\MTM_HRMS_LIB\\0_Main_Lib\\MTM_HRMS_MAINLIB.xlsx")[10,5])

HRMF_output <- HRMF(file = file, formula = formula, charge = 1, mass_accuracy = 5, intensity_cutoff = 50000, IR_RelAb_cutoff = 1)
HRMF_scores <- HRMF_output$HRMF_total
HRMF_compound <- HRMF_output$compound
HRMF_allions <- HRMF_output$all_ions

# Theoretical 2,6-TDI by QCxMS

file <- "D:\\Projects\\PredictEIspectra\\QCxMS_Runs\\2_6-Diisocyanatotoluene\\2_6-Diisocyanatotoluene.msp"

formula <- as.character(readxl::read_xlsx("D:\\R_projects\\MTM_HRMS_LIB\\0_Main_Lib\\MTM_HRMS_MAINLIB.xlsx")[10,5])

HRMF_output <- HRMF(file = file, formula = formula, charge = 1, mass_accuracy = 5, intensity_cutoff = 1, IR_RelAb_cutoff = 1)
HRMF_scores <- HRMF_output$HRMF_total
HRMF_compound <- HRMF_output$compound
HRMF_allions <- HRMF_output$all_ions


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
file <- "./data/benzothiazole_dust.msp"
formula <- "C7H5NS" 

HRMF_output <- HRMF(file = file, formula = formula)
HRMF_scores <- HRMF_output$HRMF_total
HRMF_compound <- HRMF_output$compound
HRMF_allions <- HRMF_output$all_ions


# msp <- "./data/CI Pigment yellow 151_QCxMS.msp"
# formula <-"C18H15N5O5" 


####

file <- "./data/Spectrum_0_202111221149_NIST_cut.msp"
formula <- "C7H5NS"

HRMF_output <- HRMF(file = file, formula = formula)

# octocrylene 
file <- "D:\\TEST\\Spectrum_0_202111221149_NIST_ID630.msp"
formula <- "C4H10O6S2"
formula <- "C4H11O2P"
HRMF_output <- HRMF(file = file, formula = formula)



library(rcdk)
mass <- 69.03353
window <- 0.000345
me <- 0.00054858

element_range <- list(c('C', 0,4), c('H', 0,11), c('O', 0,2), c('P', 0,1))
match <- rcdk::generate.formula(
  mass+me,
  window = window,
  elements = element_range,
  validation = FALSE,
  charge = 1
)


mass - window
mass + window
chemical <- get.formula("C4H5O", charge = 1)
chemical@mass




