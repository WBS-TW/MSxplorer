test_that("plotting works", {
  
  library(readxl)
  filepath <-  "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/QA_QC_1_F,3_1.mzML"
  filepath <-  "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/B6 batch std_1_F,2_1.mzML"
  filepath <-  "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/10Rwater A_1_C_8_1.mzML"
  filepath <-  "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/19Rlarvae_1_C,3_1.mzML"
  featlist <- readxl::read_xlsx("D:\\R_projects\\MSXploreR\\tests\\featlist.xlsx", sheet = "PFSAs")
  
  featlist <- readxl::read_xlsx("D:\\R_projects\\MSXploreR\\tests\\featlist.xlsx", sheet = "novel")
  
  # featlist <- data.frame(name = c("susp1", "susp2", "susp3", "susp4", "susp5", "susp6"), 
  #                        mz = c(247.9618, 298.9429, 505.0516, 348.9398, 79.9630, 79.9630), 
  #                        ms_level = c("ms1", "ms1", "ms1", "ms1", "ms1", "ms2"))
  
  
  
  plotEIC(filepath = filepath, featlist = featlist)
  plotEIC(filepath = filepath, featlist = featlist, diff = 0.01, mserr = 0.01)
  
})
