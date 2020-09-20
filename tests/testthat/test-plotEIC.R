test_that("plotting works", {
  
  
  featlist <- data.frame(name = c("susp1", "susp2", "susp3", "susp4", "susp5", "susp6"), 
                         mz = c(247.9618, 298.9429, 505.0516, 348.9398, 79.9630, 79.9630), 
                         ms_level = c("ms1", "ms1", "ms1", "ms1", "ms1", "ms2"))
  
  filepath <-  "D:/Raw_data/Kallinge/New_analysis_20200414/centroid/QA_QC_1_F,3_1.mzML"
  
  plotEIC(filepath = filepath, featlist = featlist)

  
})
