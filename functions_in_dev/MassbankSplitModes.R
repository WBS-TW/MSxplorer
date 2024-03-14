

msp_file <- "./data/MassBankEU_NIST_cropped.msp"

# read msp file
m <- readLines(msp_file, warn = FALSE)

# change comment field to avoid problem with "Names:"
# comment <- grep("Comments:", m)
# m <- m[-(comment)]

m <- stringr::str_replace_all(m, "Comments:.*(Name:)", "Comment removed due to incorrect name format")

n <- grep("Name:", m)

# get all ion modes
im <- grep("Ion_mode:", m)
#ion_modes <- unlist(lapply(im, function(x) substring(m[x], 11, nchar(m[x]))))


# Writing positive mode
for (i in seq_along(n)) {
  if (i < length(n)) {
    ind <- m[seq(n[i],n[i+1]-1)]
    im <- grep("Ion_mode:", ind)
    ion_mode <- unlist(lapply(im, function(x) substring(ind[x], 11, nchar(ind[x]))))
  } else if (i == length(n)) {
    ind <- m[seq(n[i], length(m))]
    im <- grep("Ion_mode:", ind)
    ion_mode <- unlist(lapply(im, function(x) substring(ind[x], 11, nchar(ind[x]))))
  }
  
  if(ion_mode=="POSITIVE" || ion_mode == "P" || ion_mode == "POS") {
    sink("MassBankEU_NIST_POSITIVE.msp", append = TRUE)
    for (j in seq_along(ind)) {
      cat(paste0(ind[j], "\n"))
    }
    sink()
  } else NULL
}

# Writing negative mode
for (i in seq_along(n)) {
  if (i < length(n)) {
    ind <- m[seq(n[i],n[i+1]-1)]
    im <- grep("Ion_mode:", ind)
    ion_mode <- unlist(lapply(im, function(x) substring(ind[x], 11, nchar(ind[x]))))
  } else if (i == length(n)) {
    ind <- m[seq(n[i], length(m))]
    im <- grep("Ion_mode:", ind)
    ion_mode <- unlist(lapply(im, function(x) substring(ind[x], 11, nchar(ind[x]))))
  }
  
  if(ion_mode=="NEGATIVE" || ion_mode == "N" || ion_mode == "NEG") {
    sink("MassBankEU_NIST_NEGATIVE.msp", append = TRUE)
    for (j in seq_along(ind)) {
      cat(paste0(ind[j], "\n"))
    }
    sink()
  } else NULL
}

