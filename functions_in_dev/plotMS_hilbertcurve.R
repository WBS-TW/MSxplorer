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
#' vignette("HilbertCurve")

BiocManager::install("RforMassSpectrometry/SpectraVis")

library(tidyverse)
library(HilbertCurve)
library(RforMassSpectrometry)
library(mzR)
library(Spectra)
library(SpectraVis)
library(IRanges)

#testing Spectra package for LC positive MSe data
fl <- "D:\\OneDrive - Linköpings universitet\\Raw_data\\Dust_Florian\\LC\\LC_pos_mzML\\IS_RS_rep3.mzML"

sp <- Spectra(fl)
sp

# check the ms1 and ms2 spectrum index
ms_1 <- filterPrecursorScan(sp, 100)
ms_2

SpectraVis::browseSpectra(sp)
SpectraVis::plotlySpectra(sp[100])
  

# testing pipe to filter ms1 and intensity

sp_filt_ms1_int <- sp |> 
  Spectra::filterMsLevel(1) |> 
  Spectra::filterIntensity(1000) |> 
  Spectra::filterRt(c(50, 60))


rt <- Spectra::rtime(sp_filt_ms1_int)
mz <- Spectra::mz(sp_filt_ms1_int[1])@listData[[1]]
int <- Spectra::intensity(sp_filt_ms1_int[1])@listData[[1]]

# for loop

rt <- numeric()
mz <- numeric()
int <- numeric()

for (i in seq_along(sp_filt_ms1_int)) {
  int0 <- Spectra::intensity(sp_filt_ms1_int[i])@listData[[1]]
  int <- append(int, int0)
}

ir <- IRanges(start, end)

hc <- HilbertCurve(1, length(int), level = 10, mode = "pixel")

hc_layer(hc, ir)


