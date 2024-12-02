# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

sMax <- 8

for(beta in 1.55)
{
  for(tMax in 8:16)
  {
    plotPlateauFit(spatialExtent
                   , temporalExtent
                   , beta
                   , sizeWLoops
                   , thermSkip
                   , spatialExtentMax = sMax
                   , temporalExtentMax = tMax)
  }
}
