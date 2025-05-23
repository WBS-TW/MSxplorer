% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRMF.R
\name{HRMF}
\alias{HRMF}
\title{HRMF: High resolution mass filtering for GC HRMS data}
\usage{
HRMF(
  file,
  formula,
  select_compound = 1,
  charge = 1,
  mass_accuracy = 5,
  intensity_cutoff = 1,
  IR_RelAb_cutoff = 1,
  return_detailed_list = FALSE
)
}
\arguments{
\item{file}{}

\item{formula}{character The chemical formula}

\item{select_compound}{number or a character}

\item{mass_accuracy}{integer The mass accuracy in ppm}

\item{intensity_cutoff}{integer  The absolute intensity cutoff stated in the msp}

\item{IR_RelAb_cutoff}{integer The relative abundance cutoff for the isotopologues relative to the monoisotopic mass}
}
\value{
list   list of outputs consisting of dataframes
}
\description{
HRMF: High resolution mass filtering for GC HRMS data
}
\examples{
hit <- HRMF(file = "D://R_projects//MSxplorer//data//2-Methoxy-5-methylaniline_C8H11NO.msp", formula = "C8H11NO", IR_RelAb_cutoff = 1, intensity_cutoff = 100000)
# TO CHECK
High resolution mass filtering for GC HRMS data
}
