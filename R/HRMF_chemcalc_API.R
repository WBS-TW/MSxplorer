#' HRMF: High resolution mass filtering for GC HRMS data
#'
#' @param file 
#' @param formula character The chemical formula
#' @param mass_accuracy integer The mass accuracy in ppm 
#' @param intensity_cutoff integer  The absolute intensity cutoff stated in the msp
#' @param IR_RelAb_cutoff integer The relative abundance cutoff for the isotopologues relative to the monoisotopic mass
#'


library(rjson)

mf="C10H20Cl2"

query <- paste("https://www.chemcalc.org/chemcalc/mf?mf=",mf,sep="")
resJSON <- fromJSON(paste(readLines(query), collapse=""))

resJSON$em
resJSON$mw
