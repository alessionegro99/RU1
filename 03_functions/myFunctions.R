source("~/projects/stepscaling/RU1/03_functions/header.R")

# fit functions

exponential <- function (par, x, boot.R, ...) par[1] * exp( -par[2] * x )

# various helpful functions

dnormWeightedSum <- function(x, means, sds, weights) {
  if (length(means) != length(sds) || length(sds) != length(weights)) {
    stop("All input lists (means, sds, weights) must have the same length")
  }

  weightedSum <- sum(mapply(function(mu, sigma, w) {
    w * dnorm(x, mean = mu, sd = sigma)
  }, means, sds, weights))

  return(weightedSum)
}

pnormWeightedSumNormalized <- function(x, means, sds, weights, norm) {
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

cdfShiftQuantile <- function(x, means, sds, weights, quantile, norm)
  pnormWeightedSumNormalized(x, means, sds, weights, norm) - quantile

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

  mEffAIC0 <- c(rep(0, spatialExtentMax))
  mEffAICBoot <- matrix(0, nrow = bootSamples, ncol = spatialExtentMax)

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
    #countNA <- length(which(is.na(mEfft$t)))

    m0Lst <- list()
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
                                        , replace.na = FALSE)

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
}

