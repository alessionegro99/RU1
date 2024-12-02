# Extracting effective mass fitting the effective mass plateau

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

t1Mat <- matrix(data = NA, nrow = 3, ncol = length(betaArray))
t2Mat <- matrix(data = NA, nrow = 3, ncol = length(betaArray))

t1Mat[1, ] <- c(4)
t2Mat[1, ] <- c(12)

t1Mat[2, ] <- c(4)
t2Mat[2, ] <- c(12)

t1Mat[3, ] <- c(5)
t2Mat[3, ] <- c(9)

distances <- c(1, 2, 3)

for(spatialExtent in spatialExtentArray){
  j <- 0
  for(beta in betaArray){
    j <- j+1
    for(r in distances){
      cat("Now L = ", spatialExtent, " beta = ", beta, " r = ", r, "\n")

      t1 <- t1Mat[r, j]
      t2 <- t2Mat[r, j]

      tmp <- extractEffectiveMass(spatialExtent
                                  , temporalExtent
                                  , beta
                                  , sizeWLoops
                                  , thermSkip
                                  , r = r
                                  , t1 = t1
                                  , t2 = t2)
      summary(tmp)

      chi2red <- tmp$effmassfit$t0[[2]]/tmp$effmassfit$dof

      saveRDS(tmp, paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_"
                          , spatialExtent
                          , "_"
                          , temporalExtent
                          ,"_"
                          , beta
                          ,"_0/plateauData_r"
                          , r,".rds"))

      pdf(paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/heatbath_"
                 , spatialExtent
                 , "_"
                 , temporalExtent
                 ,"_"
                 , beta
                 ,"_0/mEffPlateau_r_"
                 , r
                 , ".pdf"))

      ylimopt = c(0.95*tmp$effmassfit$t0[[1]],(1 + 0.15*r)*tmp$effmassfit$t0[[1]])

      plotMEffPlatWithLegend(tmp
                             , beta
                             , t1
                             , t2
                             , spatialExtent
                             , temporalExtent
                             , chi2 = chi2red
                             , ylim = ylimopt)

      dev.off()
    }
  }
}
