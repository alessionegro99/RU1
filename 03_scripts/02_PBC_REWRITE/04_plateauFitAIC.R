# Extracting effective masses using AIC

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

sMax <- 8

for(beta in betaArray)
{
  print(paste("Currently beta =", beta))
  for(tMax in 12:16)
  {
    mEff <- computeEffectiveMassAIC(spatialExtent
                            , temporalExtent
                            , beta
                            , sizeWLoops
                            , thermSkip
                            , spatialExtentMax = sMax
                            , temporalExtentMax = tMax)

    saveRDSToData(fileToSave = mEff
                  , saveNameNoExt = paste0("mEfftMax", tMax, "sMax", sMax)
                  , spatialExtent = spatialExtent
                  , temporalExtent = temporalExtent
                  , invCoupling = beta
                  , sizeWLoops = sizeWLoops)
  }
}


