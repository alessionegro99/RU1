library(hadron)
library(latex2exp)
library(purrr)

# simulations parameters

spatialExtent <- 32
temporalExtent <- 16
sizeWLoops <- 1
betaArray <- seq(1.55, 2, 0.05)

# values to be used in the subsequent scripts

thermSkip <- 1000
bootSamples <- 500
blockSizeAnalysis <- 2
blockSize <- 100

# filenames

inputFileName <- function(spatialExtent, temporalExtent, invCoupling, sizeWLoops, extra)
{
  if(missing(extra))
    extra <- ""
  paste0("heatbathinput_", spatialExtent, "_", temporalExtent, "_3_", invCoupling,"_nover5_planartrue_size", sizeWLoops, extra)
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





