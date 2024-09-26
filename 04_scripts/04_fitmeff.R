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
