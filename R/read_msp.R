#' read_msp
#'
#' @param file, a string with the path of the msp file in either NIST or RIKEN format
#'
#' @return a list of dataframes from individual compounds in the msp file containing mz, intensity and annotation (if available)
#' @export
#'
#' @examples 
#' file <- "./data/GC_orbitrapSTD_RIKEN.msp"
#' peak_table_list <- read_msp(file)

read_msp <- function(file) {
  
# reading msp file
msp <- readLines(file, warn = FALSE)

# Removing comments that include Name: to avoid incorrect indexing

msp <- stringr::str_replace_all(msp, "Comments:.*(Name:)", "Comment removed due to incorrect name format")

Names <- grep("Name:", msp, ignore.case = TRUE)

msp <- split(msp, cumsum(seq_along(msp) %in% Names))

peak_tbl_list <- list()

for (i in seq_along(msp)) {
  
m_ind <- msp[[i]]

Name <- stringr::str_extract(m_ind[1], "(?m)(?<=\\b(?i)Name:).*$")
Name <- stringr::str_replace(Name, " ", "")
numpeaks <- grep("Num Peaks", m_ind, ignore.case = TRUE)
peaks <- m_ind[numpeaks+1:length(m_ind)]
peaks <- na.omit(peaks)
peak_list <- as.list(peaks)

peak_tbl <- data.frame(mz = NA, intensity = NA, annotation = NA)

maxlength <- 3L

# check that some peak tables are separated by tab or space
if(grepl("\t", peak_list[[1]]) == TRUE){
  for (j in seq_along(peak_list)) {
    peakrow <- unlist(strsplit(peak_list[[j]], "\t"))
    peakrow <- c(peakrow, rep(NA, maxlength-length(peakrow)))
    peak_tbl <- rbind(peak_tbl, peakrow)
    }
  }else if (grepl("\t", peak_list[[1]]) == FALSE) {
    for (j in seq_along(peak_list)) {
      peakrow <- unlist(stringr::str_split(peak_list[[j]], " ", n = 3))
      peakrow <- c(peakrow, rep(NA, maxlength-length(peakrow)))
      peak_tbl <- rbind(peak_tbl, peakrow)
    }
    }

find_allNA <- which(apply(peak_tbl, 1, function(x)all(is.na(x))))

peak_tbl <- peak_tbl[-find_allNA,]

peak_tbl$mz <- as.numeric(peak_tbl$mz)
peak_tbl$intensity <- as.numeric(peak_tbl$intensity)


peak_tbl_list <- append(peak_tbl_list, list(peak_tbl))
names(peak_tbl_list)[i] <- Name


}

# rownames(peak_tbl_list) <- NULL

return(peak_tbl_list)

}
