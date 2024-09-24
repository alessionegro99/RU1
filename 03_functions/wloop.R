source("~/projects/stepscaling/RU1/03_functions/header.R")

# various helpful functions

# extracts the Wilson loop(s) timeseries and puts it into a cf-type object
wlooptocf <- function (path, file, L, T, sizewloops, r, skiprows){
  # set thermalization time to 0 if not given
  if (missing(skiprows)) skiprows <- 0

  # reading the file into a matrix type object and extracting its dimension
  tmp <- as.matrix(read.table(paste(path, file, sep = "")))
  tmp <- tmp[-(1:skiprows), ]
  matrixdim <- dim(tmp)

  # number of different loops measured W(r,t)  = number of columns
  nconfs <- matrixdim[2] - 1

  # maximum of the extents of the Wilson loops measured
  nsmax = sizewloops * L
  ntmax = sizewloops * T

  # extracting <W(r,t)> for fixed r
  colindex <- seq(from = r, to = (ntmax-1)*nsmax + r, by = nsmax)

  # creating cf object
  newcf <- cf_orig(cf = tmp[, colindex])
  newcf <- cf_meta(newcf, nrObs = 1, Time = ntmax)

  return(invisible(newcf))
}
