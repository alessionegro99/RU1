# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

tMax = 14
rMax = 3

invCoupling = 1.55

temporalExtentMax = tMax
spatialExtentMax = rMax

############################################

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

mEffAICWeight0 <- c(rep(0, spatialExtentMax))
errAICWeight0 <- c(rep(0, spatialExtentMax))
mEffAICWeightBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)
errAICWeightBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)

mEffAICQuantile0 <- c(rep(0, spatialExtentMax))
errAICQuantile0 <- c(rep(0, spatialExtentMax))
mEffAICQuantileBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)
errAICQuantileBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)

pdf(paste0(plotPath, "/mEffDistributions.pdf"))
for(i in 1:spatialExtentMax)
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
  mBootLst <- list()

  mSq0Lst <- list()
  mSqBootLst <- list()

  dm0Lst <- list()

  w0Lst <- list()
  wBootLst <- list()

  for (j in 1:(temporalExtentMax - 1 - 2))
  {
    for (k in (j + 2):(temporalExtentMax - 1))
    {
      if(is.na(mEfft$effMass[[j]]) == FALSE
         && is.na(mEfft$effMass[[k]]) == FALSE
         && is.na(mEfft$deffMass[[j]]) == FALSE
         && is.na(mEfft$deffMass[[k]]) == FALSE
         && mEfft$effMass[[j]] > 0
         && mEfft$effMass[[k]] > 0
         && abs(mEfft$effMass[[j]]/mEfft$deffMass[[j]]) > 1
         && abs(mEfft$effMass[[k]]/mEfft$deffMass[[k]]) > 1
         && length(which(is.na(mEfft$t[, j]))) == 0
         && length(which(is.na(mEfft$t[, k]))) == 0
      )
      {
        fitmt1t2 <- fit.effectivemass(mEfft
                                      , t1 = j
                                      , t2 = k
                                      , useCov = FALSE
                                      , replace.na = TRUE)

        m0Lst[[length(m0Lst) + 1]] <- fitmt1t2$effmassfit$t0[[1]]
        mBootLst[[length(mBootLst) + 1]] <- fitmt1t2$effmassfit$t[, 1]

        dm0Lst[[length(dm0Lst) + 1]] <- fitmt1t2$effmassfit$se

        mSq0Lst[[length(mSq0Lst) + 1]] <- fitmt1t2$effmassfit$t0[[1]] ** 2
        mSqBootLst[[length(mSqBootLst) + 1]] <- fitmt1t2$effmassfit$t[, 1] ** 2

        # computing AIC weights

        w0Lst[[length(w0Lst) + 1]] <- exp(-0.5 * (fitmt1t2$effmassfit$t0[[2]] + 2 + 2 * (k - j)))
        wBootLst[[length(wBootLst) + 1]] <- exp(-0.5 * (fitmt1t2$effmassfit$t[, 2] + 2 + 2 * (k - j)))
      }
    }
  }

  if (length(m0Lst) == 0)
  {
    print(paste0("No meaningful fits available for r = "
                 , i
                 , " interrupting procedure ..."))
    break
  }

  # a global shift does not

  lapply(w0Lst, function(x, subtAmnt) x - subtAmnt, subtAmnt = min(unlist(w0Lst)))
  #wBootLst <- subtractColMins(wBootLst)

  normZ0 <- Reduce(`+`, w0Lst)
  normZBoot <- Reduce(`+`, wBootLst)

  # weighting procedure to compute the effective mass
  # first on the original data, then on the bootstrap samples

  wm0Lst <- Map(`*`, m0Lst, w0Lst)
  wmBootLst <- Map(`*`, mBootLst, wBootLst)

  wmSq0Lst <- Map(`*`, mSq0Lst, w0Lst)
  wmSqBootLst <- Map(`*`, mSqBootLst, wBootLst)

  m0 <- Reduce(`+`, wm0Lst)
  mBoot <- Reduce(`+`, wmBootLst)

  mSq0 <- Reduce(`+`, wmSq0Lst)
  mSqBoot <- Reduce(`+`, wmSqBootLst)

  mEffAICWeight0[i] <- m0/normZ0
  errAICWeight0[i] <- mSq0/normZ0 - (m0/normZ0) ** 2
  mEffAICWeightBoot[, i] <- mBoot/normZBoot
  errAICWeightBoot[, i] <- 1/normZBoot * mSqBoot - (mBoot/normZBoot) ** 2

  # quantile procedure to compute the effective mass
  # first on the original data, then on the bootstrap samples

  intervalBoot <- c(0.5*min(unlist(mBootLst)), 1.5*max(unlist(mBootLst)))
  mEff0Tmp <- uniroot(f = cdfShiftQuantile
                      , quantile = 0.5
                      , means = m0Lst
                      , sds = dm0Lst
                      , weights = w0Lst
                      , norm = normZ0
                      , interval = intervalBoot
                      , tol = 1e-12)
  mEffBootTmp <- Map(findRootCdfShiftQuantileBoot
                     , mBoot = transpose(mBootLst)
                     , wBoot = transpose(wBootLst)
                     , normBoot = normZBoot
                     , MoreArgs = list(quantile = 0.5
                                       , intervalUniRoot = intervalBoot
                                       , stDevBoot = dm0Lst))

  q160Tmp <- uniroot(f = cdfShiftQuantile
                     , quantile = 0.16
                     , means = m0Lst
                     , sds = dm0Lst
                     , weights = w0Lst
                     , norm = normZ0
                     , interval = intervalBoot
                     , tol = 1e-12)
  q16BootTmp <- Map(findRootCdfShiftQuantileBoot
                    , mBoot = transpose(mBootLst)
                    , wBoot = transpose(wBootLst)
                    , normBoot = normZBoot
                    , MoreArgs = list(quantile = 0.16
                                      , intervalUniRoot = intervalBoot
                                      , stDevBoot = dm0Lst))


  q840Tmp <- uniroot(f = cdfShiftQuantile
                     , quantile = 0.84
                     , means = m0Lst
                     , sds = dm0Lst
                     , weights = w0Lst
                     , norm = normZ0
                     , interval = intervalBoot
                     , tol = 1e-12)
  q84BootTmp <- Map(findRootCdfShiftQuantileBoot
                    , mBoot = transpose(mBootLst)
                    , wBoot = transpose(wBootLst)
                    , normBoot = normZBoot
                    , MoreArgs = list(quantile = 0.84
                                      , intervalUniRoot = intervalBoot
                                      , stDevBoot = dm0Lst))

  err0Tmp <- abs(0.5 * (q160Tmp$root - q840Tmp$root))
  errBootTmp <- abs(0.5 * (unlist(q16BootTmp) - unlist(q84BootTmp)))

  mEffAICQuantile0[i] <- mEff0Tmp$root
  mEffAICQuantileBoot[, i] <- unlist(mEffBootTmp)

  errAICQuantile0[i] <- err0Tmp
  errAICQuantileBoot[, i] <- errBootTmp

  # plotting the sum of pnorms and dnorms for visualization of the results
  # assuming gaussian distributions of the fit results

  xPlot <- seq(mEffAICQuantile0[i] - 3 * errAICQuantile0[i]
               , mEffAICQuantile0[i] + 3 * errAICQuantile0[i]
               , length.out = 1000)

  # first plotting the pdf
  yPlot <- sapply(xPlot, function(x) dnormWeightedSum(x, m0Lst, dm0Lst, w0Lst))
  plot(xPlot
       , yPlot
       , type = "l"
       , main = "Weighted Sum of Normal Distributions"
       , xlab = "x"
       , ylab = "Weighted Density"
       , col = "black"
       , lwd = 2)

  # then plotting the cdf
  yPlot <- sapply(xPlot
                  , function(x) pnormWeightedSumNormalized(x
                                                           , m0Lst
                                                           , dm0Lst
                                                           , w0Lst
                                                           , normZ0))
  plot(xPlot
       , yPlot
       , type = "l"
       , main = "Weighted Sum of Cumulative Normal Distributions"
       , xlab = "x"
       , ylab = "Weighted CDF"
       , col = "black"
       , lwd = 2)

  q50x <- mEffAICQuantile0[i]
  q50y <- pnormWeightedSumNormalized(x = q50x
                                     , m0Lst
                                     , dm0Lst
                                     , w0Lst
                                     , normZ0)
  q16x <- q160Tmp$root
  q16y <- pnormWeightedSumNormalized(x = q16x
                                     , m0Lst
                                     , dm0Lst
                                     , w0Lst
                                     , normZ0)
  q84x <- q840Tmp$root
  q84y <- pnormWeightedSumNormalized(x = q84x
                                     , m0Lst
                                     , dm0Lst
                                     , w0Lst
                                     , normZ0)

  abline(v = q50x, col = "red")
  abline(v = q16x, col = "red")
  abline(v = q84x, col = "red")
  abline(h = q50y, col = "red")
  abline(h = q16y, col = "red")
  abline(h = q84y, col = "red")
}
dev.off()

weightingProcedure <- list(res0 = mEffAICWeight0, err0 = sqrt(errAICWeight0), rest = mEffAICWeightBoot, errt = sqrt(errAICWeightBoot))
quantileProcedure <- list(res0 = mEffAICQuantile0, err0 = errAICQuantile0, rest = mEffAICQuantileBoot, errt = errAICQuantileBoot)
mEffAIC <- list(wP = weightingProcedure, qP = quantileProcedure)

labStep = 3

labelli <- c(seq(1, temporalExtent - labStep, by = labStep), paste0("T = ", temporalExtent), seq(temporalExtent - labStep, 1, by = labStep))

