library(tidyverse)
library(glue)

# Define parameters
params <- tibble::tibble(
  exe      = "MSPepSearch64.exe",
  query    = "Gusviqh",
  zi       = 0.01,
  zipppm   = 50,
  mppm     = 50,
  mzlimits = "50 -1",
  minmf    = 500,
  hits     = 10,
  lib      = "D:/Program/NIST20/MSSEARCH/MTM_HRMS_RECETOX_THERMO_20220518",
  inp      = "D:/Tempfolder/NIST_EI/multiNIST_EI.msp",
  outmgf   = "D:/Tempfolder/test.mgf",
  outtab   = "D:/Tempfolder/test.tsv"
)

# Build command
cmd <- glue(
  "{params$exe} {params$query} /ZI {params$zi} /ZIPPM {params$zipppm} /MPPM {params$mppm} ",
  "/MzLimits {params$mzlimits} /MinMF {params$minmf} /OnlyFound /HITS {params$hits} ",
  "/LIB {params$lib} /INP {params$inp} /OutChemForm /OUTMGF {params$outmgf} /OUTTAB {params$outtab}"
)

# Run command
system(cmd)


file_path <- "D:/Tempfolder/test.tsv"

# Read all lines
lines <- read_lines(file_path)

# Find start and end of the table
start_idx <- which(str_detect(lines, "^Unknown\\s+Rank"))  # header line
end_idx <- which(str_detect(lines, "^>")) %>% keep(~ . > start_idx) %>% first() - 1

# Extract only table lines
table_lines <- lines[start_idx:end_idx]

# Read into tibble
results <- read_tsv(
  I(table_lines),
  col_names = TRUE,
  na = ""
)

print(results)