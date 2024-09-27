# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

tMax <- 16
rMax <- 3

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

mEffAIC0 <- c(rep(0, rMax))
mEffAICBoot <- matrix(0, nrow = bootSamples, ncol = rMax)

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

  mEfft <- bootstrap.effectivemass(Wt, type = "log")

  m0Lst <- list()
  w0Lst <- list()
  p0Lst <- list()
  mBootLst <- list()
  wBootLst <- list()
  pBootLst <- list()

  for (j in 1:(tMax - 1 - 2))
  {
    for (k in (j + 2):(tMax - 1))
    { ##ADD NUMBER OF NANS IN THE BSSAMPLES EQUAL TO ZERO!!!
      if(is.na(mEfft$effMass[[j]]) == FALSE
         && is.na(mEfft$effMass[[k]]) == FALSE
         && is.na(mEfft$deffMass[[j]]) == FALSE
         && is.na(mEfft$deffMass[[k]]) == FALSE
         && mEfft$effMass[[j]] > 0
         && mEfft$effMass[[k]] > 0
         && abs(mEfft$effMass[[j]]/mEfft$deffMass[[j]]) > 1
         && abs(mEfft$effMass[[k]]/mEfft$deffMass[[k]]) > 1
         )
      {
        fitmt1t2 <- fit.effectivemass(mEfft
                                      , t1 = j
                                      , t2 = k
                                      , useCov = FALSE
                                      , replace.na = FALSE)
        ## ADD CONDITION ON FITTED PLATEAUX
        mt1t20 <- fitmt1t2$effmassfit$t0[[1]]
        chi20 <- fitmt1t2$effmassfit$t0[[2]]
        mt1t2Boot <- fitmt1t2$effmassfit$t[, 1]
        chi2Boot <- fitmt1t2$effmassfit$t[, 2]

        wt1t20 <- exp(-0.5 * (chi20 + 2 + 2 * (k - j)))
        pt1t20 <- mt1t20 * wt1t20
        wt1t2Boot <- exp(-0.5 * (chi2Boot + 2 + 2 * (k - j)))
        pt1t2Boot <- mt1t2Boot * wt1t2Boot

        w0Lst[[length(w0Lst) + 1]] <- wt1t20
        p0Lst[[length(p0Lst) + 1]] <- pt1t20
        wBootLst[[length(wBootLst) + 1]] <- wt1t2Boot
        pBootLst[[length(pBootLst) + 1]] <- pt1t2Boot
      }
    }
  }

  p0 <- Reduce(`+`, p0Lst)
  normZ0 <- Reduce(`+`, w0Lst)
  pBoot <- Reduce(`+`, pBootLst)
  normZBoot <- Reduce(`+`, wBootLst)

  mEffAIC0[i] <- p0/normZ0
  mEffAICBoot[, i] <- pBoot/normZBoot
}

mEffAIC <- list(t0 = mEffAIC0, t = mEffAICBoot)
