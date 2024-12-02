# plotting the MC history for every Wilson loop up to nsmax*ntmax extent

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

for(spatialExtent in spatialExtentArray){
  cat("Currently L ", spatialExtent, "\n")
  for(beta in betaArray)
  {
    cat("Plotting MC history for the plaquette of beta ", beta, "\n")
    plotMCHistory(spatialExtent
                  , temporalExtent
                  , sizeWLoops
                  , beta
                  , smallestOnly = TRUE
                  , extra = 0)
    cat("Bootstrap analysis plaquette for ", beta, "\n")
    bootstrapAnalysis(spatialExtent
                      , temporalExtent
                      , sizeWLoops
                      , beta
                      , thermSkip
                      , smallestOnly = TRUE)
  }
}




