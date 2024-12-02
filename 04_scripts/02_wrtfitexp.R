# fitting a function to the Wilson loop

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

spatialExtent <- spatialExtentArray[1]

for(beta in betaArray){
  wrt <- readRDS(paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_",spatialExtent,"_", temporalExtent,"_", beta, "_0/wrt.rds"))
  fn <- function(par, x, boot.r, ...) par[1] * exp(-par[2] * x) + par[3] * exp(-par[2]*(temporalExtent - x))
  fncosh <- function(par, x, boot.r, ...) par[1] * exp(-par[2] * x)

  rWloop <- 4
  usecosh <- FALSE

  wt <- wrt[[rWloop]]

  iLim <- c(5)
  fLim <- c(16)

  pdf(paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/heatbath_",spatialExtent,"_", temporalExtent,"_", beta, "_0/fits_r_", rWloop,".pdf"))

  for (i in seq_along(iLim)){
    plotwitherror(x = c(1:wt$Time)
                  , y = wt$cf.tsboot$t0
                  , dy = wt$tsboot.se
                  , log = "y"
                  , xlab = "t"
                  , ylab = paste0("W(r = ", rWloop ,", t)")
                  , main = paste0("W(r = ", rWloop ,", t), L = 3 , T = 16, invCoupling = ", beta))

    START <- iLim[i]
    END <- fLim[i]

    x <- c(START:END)
    value <- wt$cf.tsboot$t0[START:END]
    bsamples <- wt$cf.tsboot$t[,START:END]

    if(usecosh == FALSE)
    fit.result <- bootstrap.nlsfit(fn, c(0.7, 0.5,0.01), value, x, bsamples, use.minpack.lm = FALSE)
    if(usecosh == TRUE)
    fit.result <- bootstrap.nlsfit(fncosh, c(1, 0.5), value, x, bsamples, use.minpack.lm = FALSE)

    plot(fit.result, rep = TRUE, col.line = "black", col.band = "grey")

    A0 <- fit.result$t0[1]
    dA0 <- fit.result$se[1]
    if(useSym == FALSE){
      B0 <- fit.result$t0[3]
      dB0 <- fit.result$se[3]
    }
    mEff <- fit.result$t0[2]
    dmEff <- fit.result$se[2]
    chi2 <- fit.result$chisqr

    decimal_places <- function(x)
      return (-floor(log10(x)))

    formatted_mEff <- sprintf(paste0("%.", decimal_places(dmEff), "f"), mEff)
    formatted_dmEff <- sprintf("%.0f", dmEff * 10^(decimal_places(dmEff)))

    formatted_A0 <- sprintf(paste0("%.", decimal_places(dA0), "f"), A0)
    formatted_dA0 <- sprintf("%.0f", dA0 * 10^(decimal_places(dA0)))

    if(usecosh == FALSE){
    formatted_B0 <- sprintf(paste0("%.", decimal_places(dB0), "f"), B0)
    formatted_dB0 <- sprintf("%.0f", dB0 * 10^(decimal_places(dB0)))
    }
    if (usecosh == FALSE)
    legend_text <- c( "f(t) == A[0] * e^{-E[0] * t} + B[0] * e^{ - E[0] * (T-t)}"
                     , paste0("mEff == ", formatted_mEff, "(", formatted_dmEff, ")")
                     , paste0("chi[red]^{2} == ", chi2/fit.result$dof)
                     , paste0("A[0] == ", formatted_A0, "(", formatted_dA0, ")")
                     , paste0("B[0] == ", formatted_B0, "(", formatted_dB0, ")")
                     )
    if(usecosh == TRUE)
    legend_text <- c( "f(t) == A[0] * cosh(E[0]*(T/2-t))"
                      , paste0("mEff == ", formatted_mEff, "(", formatted_dmEff, ")")
                      , paste0("chi[red]^{2} == ", chi2/fit.result$dof)
                      , paste0("A[0] == ", formatted_A0, "(", formatted_dA0, ")")
    )


    legend("topright", legend = parse(text = legend_text))

  }

  dev.off()

  saveRDS(fit.result, paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_", spatialExtent, "_", temporalExtent, "_", beta, "_0/wrtFitResults_r", rWloop,".rds"))

  summary(fit.result)
}
