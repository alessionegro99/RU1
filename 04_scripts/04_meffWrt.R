# plotting m_eff for every r and t

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

  pdf(paste0(plotpath, "/meff.pdf"))
  for(i in 1:nsmax){

    Wt <- wlooptocf(datapath, datafile, L, T, sizewloops = frac, i, skiprows = therm_skip)
    Wt <- bootstrap.cf(Wt, boot.R = bssamples, boot.l = blocksize)

    meff <- bootstrap.effectivemass(Wt, type = "log")

    plot(meff,
         xlab = TeX(r"($m_{eff}$)"),
         ylab = TeX(r"(t)"),
         main = TeX(sprintf(paste(r'($m_{eff}, type = "log", r =$)', i))))
  }
  dev.off()
}
