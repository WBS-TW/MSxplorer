

#' TEST spectral library search using MSPepSearch
#'
#' @param msp 
#' @param arguments 
#'
#' @return
#' @export
#'
#' @examples
msp_libsearch <- function(msp, arguments) {

path <- "D:/Program/MSPepSearch/"
file <- "MSPepSearch64.exe"
lib <- "D:\\Program\\NIST14\\MTM_HRMS_RECETOX_THERMO_20220518"
msp <- "D:\\Raw_data\\Dust_Florian\\GC\\Raw_data\\mzML\\Spectrum_0_202111221149_NIST.msp"
# msp <- "\\MTM_HRMS_RECETOX_THERMO_20220518"

arguments <- c("G", "u", "s", "v", "i", "q", "h", 
               "/ZI", "0.1", 
               "/ZIPPM", "10", 
               "/MPPM", "30", 
               "/MzLimits", "50", "-1", 
               #"/RI", "t60r20", # include retention index, not working yet!
               "/MinMF", "500",
               "/OnlyFound", 
               "/HITS", "10", 
               "/OutChemForm",
               "/OutSrchComment",
               #"/PATH", "D:\\Program\\NIST14",
               "/LIB", lib, 
               #"/MAIN", "mainlib",
                "/INP", msp
               )



msp_in <- sys::exec_internal(cm = paste0(path, file), args = arguments)

output_err <- unlist(stringr::str_split(rawToChar(msp_in$stderr), "\\r\\n"))
output_list <- rawToChar(msp_in$stdout)
output_list <- unlist(stringr::str_split(output_list, ">"))[[4]]
output_table <- data.table::fread(output_list)

return(output_table)

}

