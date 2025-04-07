# fitting a function to the Wilson loop

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

r1 <- 1
r2 <- 2
r3 <- 3

spatialExtent <- spatialExtentArray[1]

for(beta in betaArray){

  cat("results for L T beta",spatialExtent, temporalExtent, beta, "\n" )

  fitr1 <- readRDS(paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_", spatialExtent, "_", temporalExtent,"_", beta, "_0/wrtFitResults_r1.rds"))
  fitr2 <- readRDS(paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_", spatialExtent, "_", temporalExtent,"_", beta, "_0/wrtFitResults_r2.rds"))
  fitr3 <- readRDS(paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_", spatialExtent, "_", temporalExtent,"_", beta, "_0/wrtFitResults_r3.rds"))

  Vr1 <- fitr1$t0[[2]]
  Vr2 <- fitr2$t0[[2]]
  Vr3 <- fitr3$t0[[2]]

  bsVr1 <- fitr1$t[,2]
  bsVr2 <- fitr2$t[,2]
  bsVr3 <- fitr3$t[,2]

  g1 <- r1^2*(Vr2-Vr1)
  g2 <- r2^2*(Vr3-Vr2)

  bsg1 <- r1^2*(bsVr2-bsVr1)
  bsg2 <- r2^2*(bsVr3-bsVr2)

  cat(g1, sd(bsg1), "\n")
  cat(g2, sd(bsg2), "\n")
  cat("\n")
}
