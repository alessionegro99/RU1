# plots V(r) from the computation of the effective mass for different values
# of tMax

sMax <- 8

tMaxi <- 12
tMaxf <- 16

tMaxArr = c(tMaxi:tMaxf)

for(beta in betaArray)
{
  source("~/projects/stepscaling/RU1/03_functions/header.R")
  source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

  inputFileName <- inputFileName(spatialExtent
                                 , temporalExtent
                                 , invCoupling = beta
                                 , sizeWLoops)

  writePath <- writePath(inputFileName)
  plotPath <- plotPath(inputFileName)

  x<-c(1:sMax)

  pal <- colorRampPalette(c("red", "blue"))

  colors <- pal(tMaxf - tMaxi)

  pdf(paste0(plotPath, "mEffvstMax.pdf"))

  mEff <- readRDS(paste0(writePath, paste0("mEfftMax", tMaxi, "sMax", sMax,".rds")))

  y <- mEff$qP$res0
  dy <- mEff$qP$err0

  plotwitherror(x
                , y
                , dy
                , col = colors[1]
                , xlab = "r"
                , ylab = "mEff"
                , main = "Comparison of different tMax in AIC")

  for(tMax in (tMaxi+1):tMaxf)
  {
    mEff <- readRDS(paste0(writePath, paste0("mEfftMax", tMax, "sMax", sMax,".rds")))

    y <- mEff$qP$res0
    dy <- mEff$qP$err0

    plotwitherror(x, y, dy, col = colors[tMax-tMaxi+1], rep = TRUE)
  }
  legend(x = "bottomright",          # Position
         legend = paste0("tMax = ", tMaxArr),  # Legend texts
         col = colors,           # Line colors
         lwd = 2)                 # Line width
  dev.off()
}
