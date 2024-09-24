library(hadron)
library(hadronhelper)
library(latex2exp)

# simulations parameters

L <- 32
T <- 16
frac <- 0.25
betaseries <- c(1.5)

# filenames

inputfilename <- function(spatialextent, temporalextent, invcoup, sizeWloops, extra){
  if(missing(extra))
    extra <- ""
  paste0("heatbathinput_", spatialextent, "_", temporalextent, "_3_", invcoup,"_nover5_planartrue_size", sizeWloops, extra)
}
datafile <- function(spatialextent, temporalextent, invcoup){
  paste0("result2p1d.u1potential.rotated.Nt", temporalextent, ".Ns", spatialextent, ".b", invcoup, ".xi1.000000.nape0.alpha1.000000finedistance")
}

# values to be used in the subsequent scripts

therm_skip <- 1000
bssamples <- 1000
blocksize_analysis <- 2
blocksize <- 100

# fit functions

exponential <- function (par, x, boot.R, ...) par[1]*exp(-par[2]*x)

