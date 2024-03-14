
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

####------##### 
# From http://sombrero.nathalievialaneix.eu/articles/e-doc-relationalSOM.html#basic-package-description
# library(tidyverse)
# library(SOMbrero)
# 
# dat <- readxl::read_xlsx("D:/Raw_data/Dust_Florian/LC/LC_neg_mzML/alignment_result_area_TEST.xlsx")
# dat_norm <- dat %>%
#   select(-c(1:3)) %>%
#   select(1,4, 8)
# dat_norm <- scale(dat_norm)
# 
# data("lesmis")
# plot(lesmis, vertex.size = 0)
# 
# # Training the SOM
# set.seed(622)
# mis.som <- trainSOM(x.data=dissim.lesmis, type = "relational", nb.save = 10,
#                     init.proto = "random", radius.type = "letremy")
# plot(mis.som, what="energy")
# 
# # Resulting clustering
# mis.som$clustering
# table(mis.som$clustering)
# 
# plot(mis.som, what = "obs", type = "names")
# plot(mis.som, what = "add", type = "graph", var = lesmis)
# plot(mis.som, what = "prototypes", type = "lines")  +
#   guides(color = guide_legend(keyheight = 0.5, ncol = 2, label.theme = element_text(size = 6))) + 
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# plot(mis.som, what = "prototypes", type = "barplot")  +
#   guides(fill = guide_legend(keyheight = 0.5, ncol = 2, label.theme = element_text(size = 6))) + 
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# plot(mis.som, what = "prototypes", type = "poly.dist")
# plot(mis.som, what = "prototypes", type = "smooth.dist")
# plot(mis.som, what = "prototypes", type = "umatrix")
# plot(mis.som, what = "prototypes", type = "mds")
# plot(mis.som, what = "prototypes", type = "grid.dist")
# 
# plot(lesmis, vertex.label.color = rainbow(25)[mis.som$clustering], 
#      vertex.size = 0)
# legend(x = "left", legend = 1:25, col = rainbow(25), pch = 19)
# 
# # Compute super clusters
# plot(superClass(mis.som))
# sc.mis <- superClass(mis.som, k = 5)
# summary(sc.mis)
# 
# table(sc.mis$cluster)
# plot(sc.mis)
# 
# plot(sc.mis, what = "prototypes", type = "grid")
# plot(sc.mis, what = "prototypes", type = "lines")
# plot(sc.mis, what = "prototypes", type = "mds")
# plot(sc.mis, type = "dendro3d")
# 
# plot(lesmis, vertex.size = 0, 
#      vertex.label.color = rainbow(5)[sc.mis$cluster[mis.som$clustering]])
# legend(x = "left", legend = paste("SC", 1:5), col = rainbow(5), pch = 19)
# 
# projectIGraph(sc.mis, lesmis)
# par(mar = rep(0,4))
# plot(sc.mis, what = "add", type = "projgraph", variable = lesmis, s.radius = 2)
