
##Use self organizing map to group MS spectra

#library(kohonen)
#library(dplyr)

# dat <- readxl::read_xlsx("D:/Raw_data/Dust_Florian/LC/LC_neg_mzML/alignment_result_area_TEST.xlsx")
# dat_norm <- dat %>%
#   select(-c(1:3)) %>%
#   select(1,4, 8)
# dat_norm <- scale(dat_norm)
#   
# 
# g <- somgrid(xdim = 5, ydim = 5, topo = "hexagonal")
# map <- som(dat_norm, grid = g, alpha = c(0.05, 0.01), radius = 1)
# plot(map, type = "changes")
# plot(map, type = "count")
# plot(map, type = "mapping")
