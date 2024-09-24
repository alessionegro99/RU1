# plotting the MC history for every Wilson loop up to nsmax*ntmax extent

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/wloop.R")

ntmax <- 1#T*frac
nsmax <- 1#L*frac

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

  pdf(paste0(plotpath, "thermalization.pdf"))
  for(i in 1:nsmax){
    Wt <- wlooptocf(datapath, datafile, L, T, frac, i)
    for(j in 1:ntmax){
      timeseries <- Wt$cf[, j]
      plot(timeseries,
           xlab = TeX(r"($n_{configs}$)"),
           ylab = TeX(r"($\langle\textit{W}(r,t)\rangle$)"),
           main = TeX(sprintf(paste(r"(MC history of W( r=)", i, r"(,t=)", j, r"(), L=)", L, r"(, T=)", T, r"(, $\beta$=)", beta))))
    }
  }
  dev.off()
}
