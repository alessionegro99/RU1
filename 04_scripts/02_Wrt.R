# plotting W_r(t) for every r and t

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

for(beta in betaArray){
  plotWrt(spatialExtent
          , temporalExtent
          , beta
          , sizeWLoops
          , thermSkip
          , performFit = FALSE)
}
