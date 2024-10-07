source("~/projects/stepscaling/RU1/03_functions/header.R")

# fit functions

saveRDSToData <- function(fileToSave
                          , saveNameNoExt
                          , spatialExtent
                          , temporalExtent
                          , invCoupling
                          , sizeWLoops)
{
  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)

  writePath <- writePath(inputFileName)

  saveRDS(fileToSave, paste0(writePath, saveNameNoExt, ".rds"))
}

exponential <- function (par, x, boot.R, ...)
{
  par[1] * exp( -par[2] * x )
}

# various helpful functions

dnormWeightedSum <- function(x, means, sds, weights)
{
  if (length(means) != length(sds) || length(sds) != length(weights)) {
    stop("All input lists (means, sds, weights) must have the same length")
  }

  weightedSum <- sum(mapply(function(mu, sigma, w) {
    w * dnorm(x, mean = mu, sd = sigma)
  }, means, sds, weights))

  return(weightedSum)
}

pnormWeightedSumNormalized <- function(x, means, sds, weights, norm)
{
  # Check that all lists are the same length
  if (length(means) != length(sds) || length(sds) != length(weights)) {
    stop("All input lists (means, sds, weights) must have the same length")
  }

  # Use mapply to compute the weighted sum of normal CDFs
  weightedSum <- sum(mapply(function(mu, sigma, w) {
    w * pnorm(x, mean = mu, sd = sigma)
  }, means, sds, weights))

  return(weightedSum/norm)
}

cdfShiftQuantile <- function(x, means, sds, weights, quantile, norm){
  pnormWeightedSumNormalized(x, means, sds, weights, norm) - quantile
}

findRootCdfShiftQuantileBoot <- function(mBoot
                                           , stDevBoot
                                           , wBoot
                                           , quantile
                                           , normBoot
                                           , intervalUniRoot)
{

  findRoot <- function(x){
    cdfShiftQuantile(x
                     , means = mBoot
                     , sds = stDevBoot
                     , weights = wBoot
                     , quantile = quantile
                     , norm = normBoot)
  }

  root <- uniroot(findRoot, interval = intervalUniRoot, tol = 1e-12)$root
  return(root)
}

# extracts the Wilson loop(s) timeSeries and puts it into a cf-type object
wLoopToCf <- function (dataPath
                       , dataFile
                       , spatialExtent
                       , temporalExtent
                       , sizeWLoops
                       , r
                       , skipRows)
{
  # set thermalization time to 0 if not given
  if (missing(skipRows)) skipRows <- 0

  # reading the dataFile into a matrix type object and extracting its dimension
  tmp <- as.matrix(read.table(paste(dataPath, dataFile, sep = "")))
  tmp <- tmp[-(1:skipRows), ]
  matrixDim <- dim(tmp)

  # number of different loops measured W(r,t)  = number of columns
  nconfs <- matrixDim[2] - 1

  # maximum of the extents of the Wilson loops measured
  nsMax <- sizeWLoops * spatialExtent
  ntMax <- sizeWLoops * temporalExtent

  # extracting <W(r,t)> for fixed r
  colIndex <- seq(from = r, to = (ntMax-1)*nsMax + r, by = nsMax)

  newcf <- cf_orig(cf = tmp[, colIndex])
  newcf <- cf_meta(newcf, nrObs = 1, Time = ntMax)

  return(invisible(newcf))
}

# plots the Monte-Carlo history of the Wilson loop observable
plotMCHistory <- function(spatialExtent
                          , temporalExtent
                          , sizeWLoops
                          , invCoupling
                          , smallestOnly)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)
  dataFile <- dataFile(spatialExtent
                       , temporalExtent
                       , betaPrecision)

  dataPath <- dataPath(inputFileName)
  plotPath <- plotPath(inputFileName)

  if (!dir.exists(plotPath))
  {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  if(smallestOnly == TRUE)
  {
    nsMax <- 1
    ntMax <- 1
  }
  else
  {
    nsMax <- sizeWLoops * temporalExtent
    ntMax <- sizeWLoops * temporalExtent
  }

  pdf(paste0(plotPath, "thermalization.pdf"))
  for(i in 1:nsMax)
  {
    wt <- wLoopToCf(dataPath
                    , dataFile
                    , spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , i)
    for(j in 1:ntMax)
    {
      timeSeries <- wt$cf[, j]
      plot(timeSeries
           , xlab = TeX(r"($n_{configs}$)")
           , ylab = TeX(r"($\langle\textit{W}(r,t)\rangle$)")
           , main = TeX(sprintf(paste(r"(MC history of W( r =)"
                                      , i
                                      , r"(, t =)"
                                      , j
                                      , r"(), L =)"
                                      , spatialExtent
                                      , r"(, T =)"
                                      , temporalExtent
                                      , r"(, $\invCoupling$ =)"
                                      , invCoupling))))
    }
  }
  dev.off()
}

# computes and plots the bootstrap error for the Wilson loop observable
bootstrapAnalysis <- function(spatialExtent
                              , temporalExtent
                              , sizeWLoops
                              , invCoupling
                              , thermSkip
                              , smallestOnly)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)
  dataFile <- dataFile(spatialExtent
                       , temporalExtent
                       , betaPrecision)

  dataPath <- dataPath(inputFileName)
  plotPath <- plotPath(inputFileName)

  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  if(smallestOnly == TRUE)
  {
    nsMax <- 1
    ntMax <- 1
  }
  else
  {
    nsMax <- sizeWLoops * temporalExtent
    ntMax <- sizeWLoops * temporalExtent
  }

  pdf(paste0(plotPath, "bootstrap.pdf"))
  for(i in 1:nsMax){
    Wt <- wLoopToCf(dataPath
                    , dataFile
                    , spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , i)
    for(j in 1:ntMax){
      timeSeries <- Wt$cf[, j]
      bootstrap.analysis(timeSeries
                         , skip = thermSkip
                         , boot.R = bootSamples
                         , boot.l = blockSizeAnalysis
                         , pl = TRUE)
      title(main = TeX(sprintf(paste(r"(Bootstrap error for W( r =)"
                                     , i
                                     , r"(, t =)"
                                     , j
                                     , r"(), r = )", spatialExtent, r"(, T =)"
                                     , temporalExtent, r"(, $\invCoupling$ =)"
                                     , invCoupling))))
    }
  }
  dev.off()
}

# plots the Wilson loop as a function of t for every value of r
plotWrt <- function(spatialExtent
                    , temporalExtent
                    , invCoupling
                    , sizeWLoops
                    , thermSkip
                    , performFit)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)
  dataFile <- dataFile(spatialExtent
                       , temporalExtent
                       , betaPrecision)

  dataPath <- dataPath(inputFileName)
  plotPath <- plotPath(inputFileName)

  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  nsMax <- spatialExtent * sizeWLoops
  ntMax <- temporalExtent * sizeWLoops

  pdf(paste0(plotPath, "/wilsonloop.pdf"))
  for(i in 1:nsMax){

    Wt <- wLoopToCf(dataPath
                    , dataFile
                    , spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , i
                    , thermSkip)
    Wt <- bootstrap.cf(Wt, boot.R = bootSamples, boot.l = blockSize)

    plot(Wt,
         xlab = TeX(r"($t$)"),
         ylab = TeX(r"($W_r(t)$)"),
         main = TeX(sprintf(paste(r"($W_r$(t))"
                                  , r"(, r =)"
                                  , i
                                  , r"(, L =)"
                                  , spatialExtent
                                  , r"(, T =)"
                                  , temporalExtent
                                  , r"(, $\invCoupling$ =)"
                                  , invCoupling))))

    # expecting exponential decrease, plotting in logscale
    plot(Wt,
         log = "y",
         xlab = TeX(r"($t$)"),
         ylab = TeX(r"($W_r(t)$)"),
         main = TeX(sprintf(paste(r"(logscale, $W_r$(t))"
                                  , r"(, r =)"
                                  , i
                                  , r"(, L =)"
                                  , spatialExtent
                                  , r"(, T =)"
                                  , temporalExtent, r"(, $\invCoupling$=)"
                                  , invCoupling))))
    if(performFit == TRUE)
    {
      fit.result <- bootstrap.nlsfit(fn = exponential,
                                       par.guess = rep(1, 2),
                                       bsamples = Wt$cf.tsboot$t,
                                       x = 1:ntMax,
                                       y = Wt$cf.tsboot$t0,
                                       CovMatrix = NULL,
                                       na.rm = TRUE)
      summary(fit.result)
      plot(fit.result,
             xlab = TeX(r"($t$)"),
             ylab = TeX(r"($W_r(t)$)"),
             main = TeX(sprintf(paste(r"($W_r$(t))"
                                      , r"(, r =)"
                                      , i
                                      , r"(, L =)"
                                      , spatialExtent
                                      , r"(, T =)"
                                      , temporalExtent
                                      , r"(, $\invCoupling$=)"
                                      , invCoupling))))
      }
    }
  dev.off()
}


# computes the effective mass from the value of Wrt
plotEffectiveMass <- function(spatialExtent
                             , temporalExtent
                             , invCoupling
                             , sizeWLoops
                             , thermSkip)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops)
  dataFile <- dataFile(spatialExtent, temporalExtent, betaPrecision)

  dataPath <- dataPath(inputFileName)
  plotPath <- plotPath(inputFileName)

  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  nsMax <- spatialExtent * sizeWLoops

  pdf(paste0(plotPath, "/meff.pdf"))
  for(i in 1:nsMax){
    Wt <- wLoopToCf(dataPath
                    , dataFile
                    , spatialExtent
                    , temporalExtent
                    , sizeWLoops
                    , i
                    , thermSkip)
    Wt <- bootstrap.cf(Wt, boot.R = bootSamples, boot.l = blockSize)

    mEfft <- bootstrap.effectivemass(Wt, type = "log")

    plot(mEfft,
         xlab = TeX(r"($m_{eff}$)"),
         ylab = TeX(r"(t)"),
         main = TeX(sprintf(paste(r"($m_eff$(t), type = log, r =)"
                           , i
                           , r"(, L =)"
                           , spatialExtent
                           , r"(, T =)"
                           , temporalExtent
                           , r"(, $\beta$ =)"
                           , invCoupling))))
  }
  dev.off()
}

# plateau extraction based on the Akaike information criterion
# ### ADD REFERENCE HERE ###

computeEffectiveMassAIC <- function(spatialExtent
                                    , temporalExtent
                                    , invCoupling
                                    , sizeWLoops
                                    , thermSkip
                                    , spatialExtentMax
                                    , temporalExtentMax)
{
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
  mEffAICWeightBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)

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
    dm0Lst <- list()
    w0Lst <- list()
    p0Lst <- list()
    mBootLst <- list()
    wBootLst <- list()
    pBootLst <- list()

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

          mt1t20 <- fitmt1t2$effmassfit$t0[[1]]
          dmt1t20 <- fitmt1t2$effmassfit$se
          chi20 <- fitmt1t2$effmassfit$t0[[2]]
          mt1t2Boot <- fitmt1t2$effmassfit$t[, 1]
          chi2Boot <- fitmt1t2$effmassfit$t[, 2]

          ## check definition 2 * (k - j)

          wt1t20 <- exp(-0.5 * (chi20 + 2 + 2 * (k - j)))
          pt1t20 <- mt1t20 * wt1t20
          wt1t2Boot <- exp(-0.5 * (chi2Boot + 2 + 2 * (k - j)))
          pt1t2Boot <- mt1t2Boot * wt1t2Boot

          m0Lst[[length(m0Lst) + 1]] <- mt1t20
          dm0Lst[[length(dm0Lst) + 1]] <- dmt1t20
          w0Lst[[length(w0Lst) + 1]] <- wt1t20
          p0Lst[[length(p0Lst) + 1]] <- pt1t20
          mBootLst[[length(mBootLst) + 1]] <- mt1t2Boot
          wBootLst[[length(wBootLst) + 1]] <- wt1t2Boot
          pBootLst[[length(pBootLst) + 1]] <- pt1t2Boot
        }
      }
    }

    normZ0 <- Reduce(`+`, w0Lst)
    normZBoot <- Reduce(`+`, wBootLst)

    # weighting procedure to compute the effective mass
    # first on the original data, then on the bootstrap samples

    p0 <- Reduce(`+`, p0Lst)
    pBoot <- Reduce(`+`, pBootLst)

    mEffAICWeight0[i] <- p0/normZ0
    mEffAICWeightBoot[, i] <- pBoot/normZBoot

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

  weightingProcedure <- list(res0 = mEffAICWeight0, rest = mEffAICWeightBoot)
  quantileProcedure <- list(res0 = mEffAICQuantile0, err0 = errAICQuantile0, rest = mEffAICQuantileBoot, errt = errAICQuantileBoot)
  mEffAIC <- list(wP = weightingProcedure, qP = quantileProcedure)

  return(mEffAIC)
}
