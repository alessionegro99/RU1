# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

# ntmax <- T*frac
# nsmax <- L*frac
#
# for(beta in betaseries){
#   betaprecision <- sprintf("%.6f", beta)
#
#   inputfilename <- inputfilename(L, T, beta, frac)
#   datafile <- datafile(L, T, betaprecision)
#
#   datapath <- datapath(inputfilename)
#   plotpath <- plotpath(inputfilename)
#   writepath <- writepath(inputfilename)
#
#   if (!dir.exists(plotpath)) {
#     dir.create(plotpath, recursive = TRUE)
#     cat("Directory created:", plotpath, "\n")
#   }
#   if (!dir.exists(writepath)) {
#     dir.create(writepath, recursive = TRUE)
#     cat("Directory created:", writepath, "\n")
#   }
#
#   pdf(paste0(plotpath, "/fitmeff.pdf"))
#   for(i in 1:nsmax){
#     write("r,t1,t2,meff,dmeff,chi2,dof,chi2red", file = paste0(writepath, "effectivemassfit_r", i,".csv"))
#
#     Wt <- wlooptocf(datapath, datafile, L, T, frac, i, skiprows = therm_skip)
#     Wt <- bootstrap.cf(Wt, boot.R = bssamples, boot.l = blocksize)
#
#     meffvec <- bootstrap.effectivemass(Wt, type = "log")
#
#     my.bssamplesfit.list <- list()
#
#     num.meffvalues <- sum(!is.na(meffvec[[2]]))
#
#     write("", file = paste0(writepath, "effectivemassfit_r", i,".csv"), append = TRUE)
#
#     for (j in 1:(num.meffvalues-3)){
#       for (k in (j+2):(num.meffvalues-1)){
#         if(is.na(meffvec$t0[j+1]) == FALSE && is.na(meffvec$t0[k+1]) == FALSE){
#           # check if the replacing of nans has a significant effect on the fit
#           effmass <- fit.effectivemass(meffvec, t1 = j, t2 = k, useCov = FALSE, replace.na = FALSE)
#
#           t1 <- effmass$effmassfit$t1
#           t2 <- effmass$effmassfit$t2
#           meff <- effmass$effmassfit$t0[[1]]
#           dmeff <- effmass$effmassfit$se
#           chi2 <- effmass$effmassfit$chisqr
#           dof <- effmass$dof
#           chi2red <- chi2/dof
#
#           line = paste(i, t1, t2, meff, dmeff, chi2red, sep = ",")
#           write(line, file = paste0(writepath, "effectivemassfit_r", i, ".csv"), append = TRUE)
#
#           my.bssamplesfit.list[[length(my.bssamplesfit.list)+1]] <- effmass$effmassfit$t[, 1]
#
#           plot(effmass,
#                 xlab = "t",
#                 ylab = TeX(r"($m_eff$)"),
#                 main = TeX(sprintf(paste0(r"(const fit to $m_eff$(t))", r"(, r =)", i, r"(, L =)", L, r"(, T =)", T, r"(, $\beta$ =)", beta))))
#         }
#       }
#     }
#     saveRDS(my.bssamplesfit.list, paste0(writepath, "bssamplesfit_r", i, ".rds"))
#   }
#   dev.off()
# }

for(beta in betaArray)
{
computeEffectiveMassAIC(spatialExtent
                        , temporalExtent
                        , beta
                        , sizeWLoops
                        , thermSkip
                        , tMax = 10
                        , rMax = 1)
}


# starting from the second point, for |t2-t1|>=2 produces every possible fit

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

tMax <- 16
rMax <- 14

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

mEffAICt0 <- c(rep(0, rMax))
mEffAICt <- matrix(0, nrow = bootSamples, ncol = rMax)

########################################

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

  # bootstrapping the AIC procedure
  mEffBoot <- list()

  for(counterBootSamples in 1:(bootSamples+1))
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

          if(counterBootSamples == (bootSamples + 1))
          {
            mt1t2 <- fitEffectiveMass$effmassfit$t0[[1]]
            chi2 <- fitEffectiveMass$chisqr[[1,1]]
          }

          else
          {
            mt1t2 <- fitEffectiveMass$effmassfit$t[[counterBootSamples,1]]
            chi2 <- fitEffectiveMass$effmassfit$t[counterBootSamples,2]
          }

          dmt1t2 <- fitEffectiveMass$effmassfit$se

          if(i==1)
            print(paste(i, j, k, counterBootSamples, mt1t2, dmt1t2))

          if(mt1t2 > 0 && mt1t2/dmt1t2 >= 1 && is.na(mt1t2) == FALSE && is.na(dmt1t2) == FALSE)
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

    if(counterBootSamples <= bootSamples)
    {
      mEffBoot[[length(mEffBoot)+1]] <- sum(wAIC)
    }
  }

  mEffAICt0[i] <- sum(wAIC)
  mEffAICt[, i] <- unlist(mEffBoot)

}

mEffLst <- list(t0 = mEffAICt0, t = mEffAICt)
