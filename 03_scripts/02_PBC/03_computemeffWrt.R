source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

# saving mEff for every r and t
for(spatialExtent in spatialExtentArray)
{
  cat("Currently at L ", spatialExtent, "\n")
  for(beta in betaArray){
    cat("Computing effective mass for beta ", beta, "\n")
    saveEffectiveMass(spatialExtent
                      , temporalExtent
                      , beta
                      , sizeWLoops
                      , thermSkip)
  }
}



