#' HRMF: High resolution mass filtering for GC HRMS data
#'
#' @param file 
#' @param formula character The chemical formula
#' @param mass_accuracy integer The mass accuracy in ppm 
#' @param intensity_cutoff integer  The absolute intensity cutoff stated in the msp
#' @param IR_RelAb_cutoff integer The relative abundance cutoff for the isotopologues relative to the monoisotopic mass
#'
#' @return list   list of outputs consisting of dataframes 
#' @export
#'
#' @examples hit <- HRMF(file = "D://R_projects//MSxplorer//data//2-Methoxy-5-methylaniline_C8H11NO.msp", formula = "C8H11NO", IR_RelAb_cutoff = 1, intensity_cutoff = 100000)


#' # TO CHECK
# # Might try out isopwrap instead to get more accurate resolution calculations. 
# check if charge mass calculations are  accurate (in isopattern use 0 or 1 charge)? Can compare with Tracefinder
# rename variables and names to clarify their functions, rename "reversed"..
# Inputs: msp file with multiple msp, a csv with different chemical formula to query for each individual msp (nrow should be same as number of msps)
# add negative and positive charge and adduct options in case of CI? But then this should only be done for the first rcdk calc and not isopattern?
# Verify correct mass for odd vs even electron ions when matching rcdk formulae 

# Check: how is the windows parameter calculated in the rcdk::get.formula function. Verify the calculations. 

# Might actually be easier to implement an API to chemcalc online "from monoisotopic mass" tool. But m/z not same as with isopattern, why?

# TO FIX
# rcdk::generate.formula does not take into account the charge (e.g. loss of electron for +1) when finding formula
# Calculations m/z of fragments does not take into account the loss of electron. Double check with chemcalc.org
# However, some fragments are also even electron ions which then do not need to account for loss of electrons?
# Isotopic masses differs between envipat and ChemCalc MF. Check why!
# apply seven golden rules, add DBE, remove negative DBE

# TODO
# 1. first Generate the monoisotopic mass -> HRMF monoiso score (charge but not added to formula by rcdk)
# 2. then add isotope patter from isopattern -> HRMF full iso score (OK)
# 2.1 Add rel ab cutoff to remove low theoretical isotopologues, (default 1% ?). Do this in isopattern (OK)
# 2.2. Filter absolute intensities, e.g >50000 (OK, added intensity_cutoff)
# 3. add score on isotopic match RMS -> % iso match score
# 3.1. Figure of merit (FoM): https://www.youtube.com/watch?v=8akAr3foa1o
# 4. Combined score 
# 4.1. scores https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8779335/

# 5. loop through top n candidates? (how to get this? from NIST? or some external file generated from NIST search)
# 6. loop through all names in msp file. How to get annotated formula?
# 7. Use the theoretical mz of all isotopes as input to xcms to perform new peak picking?



#' High resolution mass filtering for GC HRMS data
# 
# https://pubs.acs.org/doi/full/10.1021/acs.analchem.5b01503
# https://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/
# subformula graph: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00776-y

# Reverse HRMF: 
# One example is applying a high-resolution mass filter (HRMF) (Kwiecien et al., 2015). HRMF calculates the percentage of fragment ions with 
# formulae which can be predicted when setting atom constraints for formula matching to only those contained in the proposed molecular formula. 
#Reverse HRMF can also be used, but this approach limits scoring to only peaks found in the library, ignoring other peaks in the experimental 
#spectra for scoring purposes. Reverse HRMF reduces influences of artifacts during deconvolution on scoring, hence reducing false negatives, 
#but may also increase false positives when experimental peaks are real and not found in the library. 



# read_msp
library(dplyr)
library(rcdk) #v3.6.0
library(enviPat)
library(stringr)
library(tidyr)
 
 
HRMF <- function(file, formula, charge = 1, mass_accuracy = 5, intensity_cutoff = 1, IR_RelAb_cutoff = 1) {

  source("./R/read_msp.R") # this one should be omitted when the package can be loaded
  data(list = "isotopes", package = "enviPat") # this is needed by isopattern to calculate the isotopic patterns
  
  compounds <- read_msp(file)
  compound <- compounds[[1]] %>% select(mz, intensity) # NOTE: need to be able to iterate through whole msp later? But need formula for all ind hits in the msp
  rownames(compound) <- NULL
  
  
  # the windows, in Da, defining the +- Da window from which the accurate mass is based on for setting the mass accuracy
  # THIS STILL CAUSES ERROR WITH generate.formula with too narrow window. Used a workaround for this in below comment on generate.formula()
  
  
  compound <- compound %>%
    mutate(mz_min = mz - mass_accuracy/10^6*compound$mz,
           mz_max = mz + mass_accuracy/10^6*compound$mz) %>%
    filter(intensity > {{intensity_cutoff}}) %>%
    select(mz, mz_min, mz_max, intensity)
  
  # using rcdk to generate all elements in the chemical formula
  atoms <- as.data.frame(rcdk::get.formula(formula)@isotopes) %>%
    mutate(number = as.numeric(number)) %>%
    mutate(mass = as.numeric(mass))
  
  # define the minimum (zero) and maximum elements for rcdk to match theoretical formula with measured masses 
  element_range <- list()
  for (i in 1:nrow(atoms)) {
    element_range[[i]] <- c(atoms[i,][[1]], 0, atoms[i,][[2]])
  }
  
  # Calculating the number of matched peaks
  
  # initiate an empty list
  match_comp <- vector(mode = "list", length = length(compound$mz))
  
  # This windows vector needs to be executed after the intensity filtering, otherwise the length is incorrect
  windows <- mass_accuracy/10^6*compound$mz
  
  ############# Trick generate.formula into generating the correct mass for charged species by adding/removing electron mass (me) ###########
  if(charge > 0){ #for positive charge
    me <- 0.00054858
    }else if (charge == 0){ #for neutral
      me <- 0
      }else { #for negative charge
        me <- -(0.00054858)
        }
   ####################################################################################################
  
  # TO CHECK: some ions are even electron ions, will the electron mass be taken into account for even vs odd ions?
  for (i in seq_along(compound$mz)) {
    match <- rcdk::generate.formula(
      compound$mz[[i]] + me, ## Workaround for generate.formula bug. Remove if bug is fixed ##
      #compound$mz[[i]],
      window = windows[[i]],
      elements = element_range,
      validation = FALSE,
      charge = charge
    )
    
    match_list <- list(annotated = logical(), annodf = data.frame(), isopat = data.frame())
    
    if(length(match) > 0) {
      match_list[["annotated"]] <- TRUE
      
      annodf <- data.frame(MonoIso_mz = match[[1]]@mass,
                           MonoIonFormula = match[[1]]@string,
                           MassError_ppm = round((compound$mz[i]-match[[1]]@mass)/match[[1]]@mass*10^6, 2))
      match_list[["annodf"]] <- annodf
      
      # extract the monoisotopic ion formula for input into isopattern
      # FIX ERROR: loss of electron not calculated for fragments!
      chemforms <- str_extract(match_list[["annodf"]][["MonoIonFormula"]], "(?<=\\[).+?(?=\\])")
      
      isopat <- as.data.frame(enviPat::isopattern(isotopes = isotopes, chemforms = chemforms,
                                         threshold = IR_RelAb_cutoff, # this one needs to be used as input variable in function
                                         charge = charge, # should be 1 charge. rcdk already already uses +1 but chemical formula should be the same, CHECK and VERIFY!
                                         emass = 0.00054858, 
                                         plotit = FALSE, 
                                         algo=1, 
                                         rel_to = 0,
                                         verbose = TRUE)) %>%
        rename(Iso_mz = 1, Abundance = 2) %>%
        mutate(Abundance = round(Abundance, 1)) %>%
        mutate(MonoIsoFormula = chemforms) %>%
        select(MonoIsoFormula, everything()) %>%
        rename_with( ~str_remove(., paste0(chemforms, ".")), .cols = -(1:3))
      
      match_list[["isopat"]] <- isopat
      
      
    } else {
      match_list[["annotated"]] <- FALSE
      match_list[["annodf"]] <- NULL
      match_list[["isopat"]] <- NULL
    }
    
    match_comp[[i]] <- match_list
  }
  
  # binding all objects in the list
  all_ions <- bind_rows(match_comp)
  
  # This provides all theoretical ions
  all_ions <- tidyr::unnest(all_ions, cols = c(annotated, annodf, isopat)) %>%
    ungroup() %>%
    filter(annotated == TRUE) %>%
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))
  
  all_elems <- all_ions %>%
    select(-c(1:7)) %>%
    mutate(IsoFormula = as.character("")) %>%
    select(`12C`, `13C`, everything())
  
  
  for (i in 1:nrow(all_elems)) {
    elem <- NULL
    for (j in 1:ncol(all_elems)) {
      if(all_elems[i,j] > 0) {
        elem_iso <- paste0("[", colnames(all_elems[i,j]), "]", all_elems[i,j])
        elem <- paste0(elem, elem_iso)
      }
    }
    all_elems$IsoFormula[i] <- elem
  }
  
  all_elems <- all_elems %>%
    select(IsoFormula)
  
  all_ions <- all_ions %>%
    mutate(Isoformula = all_elems$IsoFormula)
  
  all_ions <- all_ions %>%
    mutate(Detected_mz = 0) %>%
    mutate(Detected_int = 0) %>%
    rename(Annotated = annotated)
  
  for (i in seq_along(all_ions$Iso_mz)) {
    for (j in seq_along(compound$mz)) {
      
      if(all_ions$Iso_mz[i] > compound$mz_min[j] & all_ions$Iso_mz[i] < compound$mz_max[j]){
        all_ions$Detected_mz[i] <- compound$mz[j]
        all_ions$Detected_int[i] <- compound$intensity[j]}
    }
  }
  
  # this mz window is for the theoretical mz to use with "reversed" HRF
  windows2 <- mass_accuracy/10^6*all_ions$Iso_mz
  
  all_ions <- all_ions %>%
    mutate(MassError_ppm = case_when(Detected_mz > 1 ~ round((Detected_mz-Iso_mz)/Iso_mz*10^6,1))) %>%
    group_by(MonoIsoFormula) %>%
    mutate(Expected_int = round(Abundance/100*max(Detected_int),0)) %>% #adding the expected abundance based on the theor isotope ratio
    mutate(Detected_RelAb = round(Detected_int/max(Detected_int)*100, 1)) %>%
    ungroup() %>%
    mutate(Iso_mz_min = Iso_mz-windows2,
           Iso_mz_max = Iso_mz+windows2) %>%
    rename(Theor_RelAb = Abundance) %>%
    filter(Expected_int > {{intensity_cutoff}}) %>%
    select(MonoIsoFormula, Iso_mz, Iso_mz_min, Iso_mz_max, Theor_RelAb, Isoformula, Detected_mz, MassError_ppm, Detected_int, Detected_RelAb, Expected_int)
  #filter(!is.na(MassError_ppm))
  
  
  compound <- compound %>% mutate(Theor_mz = 0, MonoIsoFormula = "", IsoFormula = "", Theor_RelAb = 0)
  
  for (i in seq_along(compound$mz)) {
    for (j in seq_along(all_ions$Iso_mz)) {
      
      if(compound$mz[i] > all_ions$Iso_mz_min[j] & compound$mz[i] < all_ions$Iso_mz_max[j]){
        compound$Theor_mz[i] <- all_ions$Iso_mz[j]
        compound$IsoFormula[i] <- all_ions$Isoformula[j]
        compound$MonoIsoFormula[i] <- all_ions$MonoIsoFormula[j]
        compound$Theor_RelAb[i] <- all_ions$Theor_RelAb[j]}
    }
  }
  
  compound <- compound %>%
    mutate(MassError_ppm = case_when(Theor_mz > 1 ~ round((mz-Theor_mz)/Theor_mz*10^6,1))) %>%
    group_by(MonoIsoFormula) %>%
    mutate(Detected_RelAb = round(intensity/max(intensity)*100, 1)) %>%
    ungroup()
  
  
  all_iso <- length(all_ions$MassError_ppm)
  
  #this generates a table with all theoretical ions and isotopologues for the formula and match with msp
  HRMF_theor_formula <- all_ions %>%
    mutate(Detected_RelAb = replace_na(Detected_RelAb, 0)) %>% #remove NaN in Detected_RelAb column. Check why NaN is generated..
    #drop_na(MassError_ppm) %>%
    summarise(peak_count_forw = all_iso,
              df_theortomsp = round(sum(!is.na(MassError_ppm))/all_iso*100, 0),
              HRMF_theor_score = round(sum(Detected_mz*Detected_RelAb)/sum(Iso_mz*Theor_RelAb),2) # check if this reversed ratio correct?
    ) 
  
  
  all_iso_cmp <- length(compound$MassError_ppm)
  
  # #this generates a table with all match from msp to theoretical ions and isotopologues from formula
  
  HRMF_msp_formula <- compound %>%
    #drop_na(MassError_ppm) %>%
    summarise(peak_count_rev = all_iso_cmp,
              df_msptotheor = round(sum(!is.na(MassError_ppm))/all_iso_cmp*100, 0),
              HRMF_msp_score = round(sum(Theor_mz*Theor_RelAb)/sum(mz*Detected_RelAb),2)
    )
  
  # calculate Figure of Merit (FoM) (https://www.youtube.com/watch?v=8akAr3foa1o)
  sumint <- 0L
  for (i in seq_along(all_ions$Theor_RelAb)) {
    absint <- abs(all_ions$Theor_RelAb[i] - all_ions$Detected_RelAb[i]) / max(all_ions$Theor_RelAb[i],all_ions$Detected_RelAb[i])
    sumint <- sumint + absint
  }
  
  sumint <- sumint/length(all_ions$Theor_RelAb)
  FoM <- round(1-sumint, 2)
  
  # combine all scores in one table
  HRMF_scores <- cbind(HRMF_theor_formula, HRMF_msp_formula, FoM)
  
  return(list(all_ions = all_ions, compound = compound, HRMF_scores = HRMF_scores))


  
}
