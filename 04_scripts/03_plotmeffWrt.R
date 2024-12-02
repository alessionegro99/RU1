source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")


# plotting mEff for every r and t
for(spatialExtent in spatialExtentArray)
{
  cat("Currently at L ", spatialExtent, "\n")
  for(beta in betaArray){
    cat("Plotting effective mass for beta ", beta, "\n")

    ylim <- c(list(c(0.1, 0.3), c(0.2, 0.4), c(0.3, 0.5)), rep(list(c(0,1)), 32))
    xlim <- c(list(c(0, 24), c(0, 16), c(0, 16)), rep(list(c(0,20)), 29))

    plotEffectiveMass(spatialExtent
                      , temporalExtent
                      , beta
                      , sizeWLoops
                      , thermSkip
                      , ylim = ylim
                      , xlim = xlim)
  }
}
