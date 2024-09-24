# plotting W_r(t) for every r and t

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/wloop.R")

ntmax <- T*frac
nsmax <- L*frac

for(beta in betaseries){
  betaprecision <- sprintf("%.6f", beta)

  inputfilename <- inputfilename(L, T, beta, frac)
  datafile <- datafile(L, T, betaprecision)

  datapath <- paste0("/home/negro/projects/stepscaling/RU1/01_rawdata/heatbath/", inputfilename,"/omeas/")
  plotpath <- paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/", inputfilename, "/")

  if (!dir.exists(plotpath)) {
    dir.create(plotpath, recursive = TRUE)
    cat("Directory created:", plotpath, "\n")
  }

  pdf(paste0(plotpath, "/wilsonloop.pdf"))
  for(i in 1:nsmax){
    Wt <- wlooptocf(datapath, datafile, L, T, frac, i, skiprows = therm_skip)
    Wt <- bootstrap.cf(Wt, boot.R = bssamples, boot.l = blocksize)

    plot(Wt,
         xlab = TeX(r"($t$)"),
         ylab = TeX(r"($W_r(t)$)"),
         main = TeX(sprintf(paste(r"($W_r$(t))", r"(, r=)", j, r"(, L=)", L, r"(, T=)", T, r"(, $\beta$=)", beta))))

    # expecting exponential decrease, plotting in logscale
    plot(Wt,
         log = "y", ,
         xlab = TeX(r"($t$)"),
         ylab = TeX(r"($W_r(t)$)"),
         main = TeX(sprintf(paste(r"(logscale, $W_r$(t))", r"(, r=)", j, r"(, L=)", L, r"(, T=)", T, r"(, $\beta$=)", beta))))

    fit.result <- bootstrap.nlsfit(fn = exponential,
                                     par.guess = rep(1, 2),
                                     bsamples = Wt$cf.tsboot$t,
                                     x = 1:ntmax,
                                     y = Wt$cf.tsboot$t0,
                                     CovMatrix = NULL,
                                     na.rm = TRUE)
    summary(fit.result)
    plot(fit.result,
           xlab = TeX(r"($t$)"),
           ylab = TeX(r"($W_r(t)$)"),
           main = TeX(sprintf(paste(r"($W_r$(t))", r"(, r=)", j, r"(, L=)", L, r"(, T=)", T, r"(, $\beta$=)", beta))))
  }
  dev.off()
}
