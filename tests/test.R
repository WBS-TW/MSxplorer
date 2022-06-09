compound[1,1] <- 215.085526
compound[1,2] <- 100
compound[2,1] <- 216.0888
compound[2,2] <- 19
compound <- compound[1:2, 1:2]

formula <- "C18H14"



########
compound <- data.frame(mz = c(263.77794, 264.78134, 265.77589, 266.77929, 267.77384, 268.77724, 269.77180, 270.77520), 
                       intensity = c(34.26, 0.75, 100, 2.2, 97.29, 2.14, 31.55, 0.69))
formula <- "C2H3Br3"


#####Compare with HRF from Tracefinder for PAHs
compound <- read.csv("./data/acenaphthylene_deconvolutedTraceFinder.csv") %>%
  select(mz, intensity)
formula <- "C12H8"
