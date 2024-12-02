
count_sum_of_squares_with_duplicates <- function(m) {
  count <- 0

  # Iterate over possible values of a and b
  for (a in 0:floor(sqrt(m))) {
    b_squared <- m - a^2
    if (b_squared >= 0) {
      b <- sqrt(b_squared)
      if (b == floor(b)) { # Check if b is an integer
        count <- count + 1
      }
    }
  }
  return(count)
}

# rsq = x*x + y*y
wLoopToCfNP <- function (dataPath
                         , dataFile
                         , spatialExtent
                         , temporalExtent
                         , rsq
                         , skip = 0)
  {
    header <- read.table(paste0(dataPath, dataFile), comment.char = "", nrows = 1)
    header <- header[2:(length(header)-1)]

    x_vals <- as.numeric(sub(".*x=([0-9]+).*", "\\1", header))
    t_vals <- as.numeric(sub(".*t=([0-9]+).*", "\\1", header))
    y_vals <- as.numeric(sub(".*y=([0-9]+).*", "\\1", header))

    u <- lapply(1:length(x_vals), function(i) c(x_vals[i], t_vals[i], y_vals[i]))

    v <- unlist(lapply(1:length(u), function(j) u[[j]][1]^2 + u[[j]][3]^2))

    ii <- which(v == rsq)

    tmp <- as.matrix(read.table(paste0(dataPath, dataFile), skip = skip))
    tmp <- tmp[, ii]

    newcf <- cf_orig(cf = tmp)
    newcf <- cf_meta(newcf, nrObs = 1, Time = temporalExtent)

    newcf$cf <- matrix(as.numeric(newcf$cf), ncol = dim(tmp)[2], nrow = dim(tmp)[1])

    return(invisible(newcf))
  }

# filenames heatbath

inputFileName <- function(spatialExtent, temporalExtent, invCoupling, extra = 0)
{
  paste0("heatbath", "_", spatialExtent, "_", temporalExtent, "_", invCoupling,"_", extra)
}

dataFile <- function(spatialExtent, temporalExtent, invCoupling)
{
  paste0("result2p1d.u1potential.rotated.Nt", temporalExtent, ".Ns", spatialExtent, ".b", invCoupling, ".xi1.000000.nape0.alpha1.000000finedistance")
}

dataPath <- function(x)
{
  paste0("/home/negro/projects/stepscaling/RU1/01_rawdata/heatbath/", x,"/omeas/")
}

plotPath <- function(x)
{
  paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/", x, "/")
}

writePath <- function(x)
{
  paste0("/home/negro/projects/stepscaling/RU1/02_output/data/",x, "/")
}

dataFileNP <- function(spatialExtent, temporalExtent, invCoupling)
{
  paste0("result2p1d.u1potential.Nt", temporalExtent, ".Ns", spatialExtent, ".b", invCoupling, ".xi1.000000.nape0.alpha1.000000nonplanar")
}

plotMEffPlatWithLegend <- function(fitEffMass
                                   , beta
                                   , t1
                                   , t2
                                   , spatialExtent
                                   , temporalExtent
                                   , chi2 = 0
                                   , ylim = NULL){

  mEff <- fitEffMass$effmassfit$t0[[1]]
  dmEff <- fitEffMass$effmassfit$se

  plot(fitEffMass
       , xlim = c(0, t2 + 5)
       , ylim = ylim
       , xlab = "t"
       , ylab = "mEff"
       , main = paste0("beta = "
                       , beta
                       , " L = "
                       , spatialExtent
                       , " T = "
                       , temporalExtent
                       , " r = "
                       , r
                       , " t1 = "
                       , t1 + 1
                       , " t2 = "
                       , t2 + 1))

  decimal_places <- -floor(log10(dmEff))  # Number of decimal places to keep in `var`

  formatted_mEff <- sprintf(paste0("%.", decimal_places, "f"), mEff)
  formatted_dmEff <- sprintf("%.0f", dmEff * 10^(decimal_places))

  legend_text <- c(paste0("mEff == ", formatted_mEff, "(", formatted_dmEff, ")")
                   , paste0("chi[red]^{2} == ", chi2))



  # Add legend with formatted text
  legend("topright", legend = parse(text = legend_text))
}

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
                                 , extra = 0)

  writePath <- writePath(inputFileName)

  saveRDS(fileToSave, paste0(writePath, saveNameNoExt, ".rds"))
}

exponential <- function (par, x, boot.R, ...)
{
  par[1] * exp( -par[2] * x )
}

# various helpful functions

plotWrtMirrored <- function(spatialExtent
                            , temporalExtent
                            , invCoupling
                            , sizeWLoops
                            , labStep)
{
  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , extra = 0)

  writePath <- writePath(inputFileName)
  plotPath <- plotPath(inputFileName)

  wrt <- readRDS(paste0(writePath, "wrt.rds"))

  x1 <- c(1:temporalExtent)
  x2 <- c((2 * temporalExtent - 1):temporalExtent)

  pdf(paste0(plotPath, "wrtMirrored.pdf"))
  for(r in c(1:(floor(spatialExtent/2))))
  {
    y1 <- wrt[[r]]$cf.tsboot$t0
    y2 <- wrt[[spatialExtent - r]]$cf.tsboot$t0

    dy1 <- apply(wrt[[r]]$cf.tsboot$t, 2, sd)
    dy2 <- apply(wrt[[spatialExtent-r]]$cf.tsboot$t, 2, sd)

    labelli <- c(seq(1, temporalExtent - labStep, labStep), paste0("T = ", temporalExtent), seq(temporalExtent- labStep, 1,  - labStep))

    plotwitherror(x = x1
                  , y = y1
                  , dy = dy1
                  , xlim = c(1, 2 * temporalExtent - 1)
                  , xaxt = 'n'
                  , col = "red"
                  , ylab = "<W(r,t)>"
                  , xlab = "t"
                  , main = paste0("<W> for r1 = ", r, ", r2 = ", spatialExtent - r, ", L = ", spatialExtent, " T = ", temporalExtent, ", beta = ", beta))
    axis(1, at = seq(1, 2 * temporalExtent - 1, 3), labels = labelli)
    plotwitherror(x = x2
                  , y = y2
                  , dy = dy2
                  , col = "blue"
                  , xlim = c(1, 2 * temporalExtent - 1)
                  , rep = TRUE)

    legend(x = "topright",
           legend = c(paste0("W(r = ", r, ")"), paste0("W(r = ",  spatialExtent - r , ")")),
           col = c("red", "blue"),
           lwd = 2)
  }
  dev.off()
}

r2F <- function(potential, r1, r2)
{
  g0 <- r1 ** 2 * (potential$res0[r2]-potential$res0[r1])/(r2-r1)
  gBoot <- r1 ** 2 * (potential$rest[, r2]-potential$rest[, r1])/(r2-r1)

  res <- list(g0 = g0, gBoot = gBoot)

  return(res)
}

subtractColMins <- function(vecLst)
{
  n <- length(vecLst[[1]])

  mins <- sapply(1:n, function(i) min(sapply(vecLst, function(v) v[i])))

  newLst <- lapply(vecLst, function(v){
    v - mins
  })

  return(newLst)
}

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

cdfShiftQuantile <- function(x, means, sds, weights, quantile, norm)
{
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

# extracts the Polyakov loop(s) timeSeries and puts it into a cf-type object
pLoopToCf <- function (dataPath
                       , dataFile
                       , spatialExtent
                       , temporalExtent
                       , sizeWLoops
                       , r
                       , skipRows = 0)
{
  # reading the dataFile into a matrix type object and extracting its dimension
  tmp <- as.matrix(read.table(paste0(dataPath, dataFile)))

  return(invisible(newcf))
}

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
                          , smallestOnly = TRUE
                          , extra)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , extra = 0)
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

      myacf <- computeacf(timeSeries, 300)
      plot(myacf)
      summary(myacf)
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
                                 , extra = 0)
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
                                 , extra = 0)
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

  Wr <- list()

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
    Wr[[length(Wr) + 1]] <- Wt
    }
  dev.off()

  return(Wr)
}

# extracts effective mass from the plateau fit

extractEffectiveMass <- function(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , sizeWLoops
                                 , thermSkip
                                 , r
                                 , t1
                                 , t2)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , extra = 0)
  dataFile <- dataFile(spatialExtent, temporalExtent, betaPrecision)

  dataPath <- dataPath(inputFileName)
  plotPath <- plotPath(inputFileName)
  writePath <- writePath(inputFileName)

  if (!dir.exists(writePath)) {
    dir.create(writePath, recursive = TRUE)
    cat("Directory created:", writePath, "\n")
  }
  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  Wt <- wLoopToCf(dataPath
                  , dataFile
                  , spatialExtent
                  , temporalExtent
                  , sizeWLoops
                  , r
                  , thermSkip)

  Wt <- bootstrap.cf(Wt, boot.R = bootSamples, boot.l = blockSize, seed = 1442556)

  mEfft <- bootstrap.effectivemass(Wt, type = "log")

  fitmt1t2 <- fit.effectivemass(mEfft
                                , t1 = t1
                                , t2 = t2
                                , useCov = FALSE
                                , replace.na = TRUE)

  return(fitmt1t2)
}

saveEffectiveMass <- function(spatialExtent
                              , temporalExtent
                              , invCoupling
                              , sizeWLoops
                              , thermSkip)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , extra = 0)
  dataFile <- dataFile(spatialExtent, temporalExtent, betaPrecision)

  dataPath <- dataPath(inputFileName)
  writePath <- writePath(inputFileName)

  if (!dir.exists(writePath)) {
    dir.create(writePath, recursive = TRUE)
    cat("Directory created:", writePath, "\n")
  }

  nsMax <- spatialExtent * sizeWLoops

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
    saveRDS(mEfft, paste0(writePath, "/mEff_r_", i, ".rds"))
  }
}

# plots mEff from existing rds file
plotEffectiveMass <- function(spatialExtent
                             , temporalExtent
                             , invCoupling
                             , sizeWLoops
                             , thermSkip
                             , ylim = NULL
                             , xlim = NULL)
{
  betaPrecision <- sprintf("%.6f", invCoupling)

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling
                                 , extra = 0)
  dataFile <- dataFile(spatialExtent, temporalExtent, betaPrecision)

  writePath <- writePath(inputFileName)
  plotPath <- plotPath(inputFileName)

  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
    cat("Directory created:", plotPath, "\n")
  }

  nsMax <- spatialExtent * sizeWLoops

  pdf(paste0(plotPath, "/mEff.pdf"))
  for(i in 1:nsMax){

    mEfft <- readRDS(paste0(writePath, "/mEff_r_", i, ".rds"))

    if (is.null(ylim) == 0)
    {
      tmpy <- ylim[[i]]
    }
    else
      tmpy <- NULL

    if (is.null(xlim) == 0)
    {
      tmpx <- xlim[[i]]
    }
    else
      tmpx <- NULL

    plot(mEfft,
         ylim = tmpy,
         xlim = tmpx,
         ylab = TeX(r"($m_{eff}$)"),
         xlab = TeX(r"(t)"),
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

# plateau extraction through consequential fits

plotPlateauFit <- function(spatialExtent
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
                                 , extra = 0)
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

  pdf(paste0(plotPath, "/mEffConstantFits.pdf"))
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

    # plots every fit for a minum of 3 support points up to temporalExtentMax

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

          plot(fitmt1t2, xlab = "t", ylab = "mEff", main = paste0("Plateau fit from", j, " to ", k))
        }
      }
    }
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
                                 , extra = 0)
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

  return(mEffAIC)
}
