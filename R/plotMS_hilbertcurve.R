#' Plot Hilbert curve for MS data using the new RforMassSpectrometry metapackage
#'
#'

#'
#' @return Hilbert curve of XXX
#' @export
#'
#' @examples
#' fl <- "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/B6 batch std_1_F,2_1.mzML"
#' featlist <- readxl::read_xlsx("D:\\R_projects\\MSXploreR\\tests\\featlist.xlsx", sheet = "PFSAs")
#' plotEIC(filepath = filepath, featlist = featlist)

library(tidyverse)
library(HilbertCurve)
library(RforMassSpectrometry)

