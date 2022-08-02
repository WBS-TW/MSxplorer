#' HRMF: High resolution mass filtering for GC HRMS data
#'
#' @param msp 
#' @param formula 
#' @param mass_accuracy 
#' @param intensity_cutoff 
#' @param IR_RelAb_cutoff 
#'
#' @return
#' @export
#'
#' @examples
#' # TO CHECK

# check if charge mass calculations are  accurate (in isopattern use 0 or 1 charge)? Can compare with Tracefinder
# rename variables and names to clarify their functions, rename "reversed"..
# Inputs: msp with multiple msp, a csv with different chemical formula to query for each individual msp (nrow should be same as number of msps)
# add negative and positive charge and adduct options in case of CI? But then this should only be done for the first rcdk calc and not isopattern?


# TODO
# 1. first Generate the monoisotopic mass -> HRMF monoiso score (OK)
# 2. then add isotope patter from isopattern -> HRMF full iso score (OK)
# 2.1 Add rel ab cutoff to remove low theoretical isotopologues, (default 1% ?). Do this in isopattern (OK)
# 2.2. Filter absolute intensities >50000 (OK)
# 3. add score on isotopic match RMS -> % iso match score
# 4. Combined score 
# scores https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8779335/

# 5. loop through top n candidates? (how to get this? from NIST? or some external file generated from NIST search)
# 6. loop through all names in msp file. How to get annotated formula?
# 7. Use the theoretical mz of all isotopes as input to xcms to perform new peak picking?



#' High resolution mass filtering for GC HRMS data

#' # use rcdk and isopattern
# https://pubs.acs.org/doi/full/10.1021/acs.analchem.5b01503
# https://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/

# read_msp
# library(dplyr)
# library(enviPat)
# library(stringr)
# library(tidyr)
# library(rcdk)

HRMF <- function(msp, formula, mass_accuracy = 10, intensity_cutoff = 50000, IR_RelAb_cutoff = 1) {

  source("./R/read_msp.R")
  data("isotopes") # this is needed by isopattern to calculate the isotopic patterns
  
  compounds <- read_msp(msp)
  compound <- compounds[[1]] %>% select(mz, intensity)
  rownames(compound) <- NULL
  
  # state the intensity cutoff
  
  #intensity_cutoff <- intensity_cutoff
  
  # state the mass accuracy in ppm
  #mass_accuracy <- mass_accuracy
  
  # the windows, in Da, defining the +- Da window from which the accurate mass is based on for setting the mass accuracy
  windows <- mass_accuracy/10^6*compound$mz
  
  compound <- compound %>%
    mutate(mz_min = mz-windows,
           mz_max = mz+windows) %>%
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
  match_comp <- list()
  
  
  for (i in seq_along(compound$mz)) {
    
    match <- rcdk::generate.formula(
      compound$mz[[i]],
      window = windows[[i]],
      elements = element_range,
      validation = FALSE,
      charge = 1
    )
    
    match_list <- list(annotated = logical(), annodf = data.frame(), isopat = data.frame())
    
    if(length(match) > 0) {
      match_list[["annotated"]] <- TRUE
      
      annodf <- data.frame(MonoIso_mz = match[[1]]@mass,
                           MonoIonFormula = match[[1]]@string,
                           MassError_ppm = round((compound$mz[i]-match[[1]]@mass)/match[[1]]@mass*10^6, 2))
      match_list[["annodf"]] <- annodf
      
      # extract the ion monoisotopic ion formula for input into isopattern
      chemforms <- str_extract(match_list[["annodf"]][["MonoIonFormula"]], "(?<=\\[).+?(?=\\])")
      
      isopat <- as.data.frame(isopattern(isotopes = isotopes, chemforms = chemforms,
                                         threshold = IR_RelAb_cutoff, # this one needs to be used as input variable in function
                                         charge = 1, # should be 1 charge but rcdk already is +1, check
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
  
  # this mz window is for the theoretical mz to use for reversed HRF
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
  
  HRMF_forward <- all_ions %>%
    mutate(Detected_RelAb = replace_na(Detected_RelAb, 0)) %>% #remove NaN in Detected_RelAb column. Check why NaN is generated..
    #drop_na(MassError_ppm) %>%
    summarise(peak_count_forw = all_iso,
              df_forw = round(sum(!is.na(MassError_ppm))/all_iso*100, 0),
              HRF_score = round(sum(Detected_mz*Detected_RelAb)/sum(Iso_mz*Theor_RelAb),2) # is this reversed ratio correct?
    ) 
  
  
  all_iso_cmp <- length(compound$MassError_ppm)
  
  # rename as this is not reversed HRMF
  HRMF_reverse <- compound %>%
    #drop_na(MassError_ppm) %>%
    summarise(peak_count_rev = all_iso_cmp,
              df_rev = round(sum(!is.na(MassError_ppm))/all_iso_cmp*100, 0),
              RHRF_score = round(sum(Theor_mz*Theor_RelAb)/sum(mz*Detected_RelAb),2)
    )
  
  
  
  HRMF_total <- cbind(HRMF_forward, HRMF_reverse)
  
  return(list(all_ions = all_ions, compound = compound, HRMF_total = HRMF_total))


  
}
