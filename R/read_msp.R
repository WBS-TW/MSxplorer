



#' read_msp
#'
#' @param file, a string with the path of the msp file 
#'
#' @return
#' @export
#'
#' @examples 
#' file <- "./data/GC_orbitrapSTD.msp"
#' peak_table_list <- read_msp(file)

read_msp <- function(file) {
  
# reading msp file
msp <- readLines(file, warn = FALSE)

# Removing comments that include Name: to avoid incorrect indexing

msp <- stringr::str_replace_all(msp, "Comments:.*(Name:)", "Comment removed due to incorrect name format")

n <- grep("Name:", msp, ignore.case = TRUE)

msp <- split(msp, cumsum(seq_along(msp) %in% n))

peak_tbl_list <- list()

for (i in seq_along(msp)) {
  
m_ind <- msp[[i]]

numpeaks <- grep("Num Peaks", m_ind, ignore.case = TRUE)
peaks <- m_ind[numpeaks+1:length(m_ind)]
peaks <- na.omit(peaks)
peak_list <- as.list(peaks)

peak_tbl <- data.frame(mz = NA, intensity = NA, annotation = NA)

maxlength <- 3L

# check that some peak tables are separated by tab or space
if(grepl("\t", peak_list[[1]]) == TRUE){
  for (i in seq_along(peak_list)) {
    i <- unlist(strsplit(peak_list[[i]], "\t"))
    i <- c(i, rep(NA, maxlength-length(i)))
    peak_tbl <- rbind(peak_tbl, i)
    }
  }else if (grepl("\t", peak_list[[1]]) == FALSE) {
    for (i in seq_along(peak_list)) {
      i <- unlist(strsplit(peak_list[[i]], " "))
      i <- c(i, rep(NA, maxlength-length(i)))
      peak_tbl <- rbind(peak_tbl, i)
    }
    }

find_allNA <- which(apply(peak_tbl, 1, function(x)all(is.na(x))))

peak_tbl <- peak_tbl[-find_allNA,]
peak_tbl$mz <- as.numeric(peak_tbl$mz)
peak_tbl$intensity <- as.numeric(peak_tbl$intensity)

peak_tbl_list <- append(peak_tbl_list, list(peak_tbl))


}


return(peak_tbl_list)

}
