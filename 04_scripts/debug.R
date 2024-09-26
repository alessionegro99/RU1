# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

tMax <- 16
rMax <- 1

invCoupling <- 1.65

betaPrecision <- sprintf("%.6f", invCoupling)

inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)
dataFile <- dataFile(spatialExtent, temporalExtent, betaPrecision)

dataPath <- dataPath(inputFileName)
plotPath <- plotPath(inputFileName)
writePath <- writePath(inputFileName)

if (!dir.exists(plotPath)) {
  dir.create(plotPath, recursive = TRUE)
  cat("Directory created:", plotPath, "\n")
}
if (!dir.exists(writePath)) {
    dir.create(writePath, recursive = TRUE)
    cat("Directory created:", writePath, "\n")
}

for(i in 1:rMax)
{
  Wt <- wLoopToCf(dataPath
                    , dataFile
                    , spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , i
                    , thermSkip)
  Wt <- bootstrap.cf(Wt, boot.R = bootSamples, boot.l = blockSize)

  effectiveMass <- bootstrap.effectivemass(Wt, type = "log")

  # due to noise some NaNs could pop up
  finiteEffectiveMasses <- sum(!is.na(effectiveMass[[2]]))

  # AIC procedure for the original data
  wAIC <- list()
  normZ <- 0

  # performing constant fits for all possible continuous ranges within
  # t = 1 and t = tMax with a minimum of at least three support points
  for (j in 0:(finiteEffectiveMasses-3))
  {
    for (k in (j+2):(finiteEffectiveMasses-1))
    {
      # performing the fit only if both points are finite!
      if(is.na(effectiveMass$t0[j+1]) == FALSE
          && is.na(effectiveMass$t0[k+1]) == FALSE)
      {
        # ### TO DO ### check if the replacing
        # of nans has a significant effect on the fit
        fitEffectiveMass <- fit.effectivemass(effectiveMass
                                                , t1 = j
                                                , t2 = k
                                                , useCov = FALSE
                                                , replace.na = FALSE)

        mt1t2 <- fitEffectiveMass$effmassfit$t0[[1]]
        dmt1t2 <- fitEffectiveMass$effmassfit$se
        chi2 <- fitEffectiveMass$chisqr[[1,1]]

        if(mt1t2>0 && mt1t2/dmt1t2>=1)
        {
          # ### TO DO: ### check 2 * (k-j) or (k-j) ??
          wAICt1t2 <- exp(-0.5 * (chi2 + 2 + 2 * (k - j)))
          # check if summing to 1
          normZ <- normZ + wAICt1t2

          wAIC[[length(wAIC) + 1]] <- wAICt1t2 * mt1t2
        }
      }
    }
  }

  wAIC <- unlist(wAIC)/normZ
  mEff0 <- sum(wAIC)

  # bootstrapping the AIC procedure
  mEffBoot <- list()

  for(counterBootSamples in 1:bootSamples)
  {
    wAIC <- list()
    normZ <- 0

    # performing constant fits for all possible continuous ranges within
    # t = 1 and t = tMax with a minimum of at least three support points
    for (j in 0:(finiteEffectiveMasses-3))
    {
      for (k in (j+2):(finiteEffectiveMasses-1))
      {
        # performing the fit only if both points are finite!
        if(is.na(effectiveMass$t0[j+1]) == FALSE
            && is.na(effectiveMass$t0[k+1]) == FALSE)
        {
          # ### TO DO ### check if the replacing
          # of nans has a significant effect on the fit
          fitEffectiveMass <- fit.effectivemass(effectiveMass
                                                  , t1 = j
                                                  , t2 = k
                                                  , useCov = FALSE
                                                  , replace.na = FALSE)

          mt1t2 <- fitEffectiveMass$effmassfit$t[[counterBootSamples,1]]
          dmt1t2 <- fitEffectiveMass$effmassfit$se
          chi2 <- fitEffectiveMass$effmassfit$t[counterBootSamples,2]

          if(mt1t2>0 && mt1t2/dmt1t2>=1)
          {
            # ### TO DO: ### check 2 * (k-j) or (k-j) ??
            wAICt1t2 <- exp(-0.5 * (chi2 + 2 + 2 * (k - j)))

            normZ <- normZ + wAICt1t2

            wAIC[[length(wAIC) +1]] <- wAICt1t2 * mt1t2
          }
        }
      }
    }
    wAIC <- unlist(wAIC)/normZ
    mEffBoot[[length(mEffBoot)+1]] <- sum(wAIC)
  }

  print(mEff0)
  print(unlist(mEffBoot)-mEff0)
}
