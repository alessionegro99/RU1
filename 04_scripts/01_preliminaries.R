# plotting the MC history for every Wilson loop up to nsmax*ntmax extent

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

for(beta in betaArray)
{
  plotMCHistory(spatialExtent
                , temporalExtent
                , sizeWLoops
                , beta
                , smallestOnly = TRUE)

  bootstrapAnalysis(spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , beta
                    , thermSkip
                    , smallestOnly = TRUE)
}


