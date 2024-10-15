# plotting W_r(t) for every r and t

for(beta in betaArray){
  source("~/projects/stepscaling/RU1/03_functions/header.R")
  source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

  tmp <- inputFileName(spatialExtent, temporalExtent, beta, sizeWLoops)
  writePath <- writePath(tmp)

  if (!dir.exists(writePath)) {
    dir.create(writePath, recursive = TRUE)
    cat("Directory created:", writePath, "\n")
  }

  wr<- plotWrt(spatialExtent
          , temporalExtent
          , beta
          , sizeWLoops
          , thermSkip
          , performFit = FALSE)

  saveRDS(wr, paste0(writePath, "wrt.rds"))
}

# plotting mirrored W(r,t) and W(R-r,T)

for(beta in betaArray)
{
  source("~/projects/stepscaling/RU1/03_functions/header.R")
  source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

  plotWrtMirrored(spatialExtent
                  , temporalExtent
                  , beta
                  , sizeWLoops
                  , labStep = 3)
}

