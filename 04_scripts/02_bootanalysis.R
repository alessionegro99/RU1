# plotting the bootstrap error for every Wilson loop up to nsmax*ntmax extent

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

  pdf(paste0(plotpath, "bootstrap.pdf"))
  for(i in 1:nsmax){
    Wt <- wlooptocf(datapath, datafile, L, T, frac, i)
    for(j in 1:ntmax){
      timeseries <- Wt$cf[, j]
      bootstrap.analysis(timeseries, skip = therm_skip, boot.R = bssamples, boot.l = blocksize_analysis, pl = TRUE)
      title(main = TeX(sprintf(paste(r"(Bootstrap error for W( r=)", i, r"(,t=)", j, r"(), L=)", L, r"(, T=)", T, r"(, $\beta$=)", beta))))
    }
  }
  dev.off()
}
