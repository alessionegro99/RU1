# Extracting effective masses using AIC

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

sMax <- 8

for(beta in betaArray)
{
  for(tMax in 8:16)
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


