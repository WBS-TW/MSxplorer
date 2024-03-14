


convert_msp_to_csv <- function(file) {

file <- "D:/Raw_data/Dust_Florian/GC/Raw_data/mzML/MSDIAL_export/Spectrum_2_2022382144_plus_MTMLibs.msp"
file <- "D:/TEST/2022-03-17_MTM_HRMS_LIB_RECETEOX_THERMO_combinedQCxMS.msp"
file <- "D:/TEST/2022-03-17_MTM_HRMS_LIB_RECETEOX_THERMO_combinedQCxMS_NISTformat.msp"
file <- "D:/TEST/NIST_msp/Spectrum_2_2022382144_combinedQCxMS_NIST.MSP"
msp <- readLines(file, warn = FALSE)

# Removing comments that include Name: to avoid incorrect indexing

msp <- stringr::str_replace_all(msp, "Comments:.*(Name:)", "Comment removed due to incorrect name format")

msp <- stringr::str_replace_all(msp, c("NAME:" = "Name:", "RETENTIONINDEX:" = "RI:"))


n <- grep("Name:", msp, ignore.case = TRUE)
RI <- grep("RI:", msp, ignore.case = TRUE)
TI_tab <- grepl("RI:\t", msp, ignore.case = TRUE)

msp2 <- msp[n]
msp2 <- stringr::str_replace_all(msp2, "Name:", "")
msp2 <- stringr::str_replace_all(msp2, ",", "_")
msp2 <- stringr::str_replace_all(msp2, "\\|", "_")
msp2 <- stringr::str_replace_all(msp2, "\\=", "_")
msp2 <- stringr::str_squish(msp2)

write(msp2, file = stringr::str_replace(file, "msp", "csv"), append = FALSE)

}