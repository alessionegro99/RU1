source("~/projects/stepscaling/RU1/03_functions/header.R")

# fit functions

exponential <- function (par, x, boot.R, ...) par[1] * exp( -par[2] * x )

# various helpful functions

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
  nsMax = sizeWLoops * spatialExtent
  ntMax = sizeWLoops * temporalExtent

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

    effectiveMass <- bootstrap.effectivemass(Wt, type = "log")

    plot(effectiveMass,
         ylim = c(0,1),
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
# (https://arxiv.org/abs/2008.01069)
computeEffectiveMassAIC <- function()
{

}
