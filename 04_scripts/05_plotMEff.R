# plotting effective mass as a function of r

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

sMax <- 8

tMaxi <- 8
tMaxf <- 16

r3 <- 3
r2 <- 2
r1 <- 1

# for(beta in 2)
# {
#   inputFileName <- inputFileName(spatialExtent
#                                  , temporalExtent
#                                  , invCoupling = beta
#                                  , sizeWLoops)
#
#   writePath <- writePath(inputFileName)
#   plotPath <- plotPath(inputFileName)
#
#   x<-c(1:sMax)
#
#   pal <- colorRampPalette(c("red", "blue"))
#
#   colors <- pal(tMaxf - tMaxi)
#
#   pdf(paste0(plotPath, "mEffvstMax.pdf"))
#
#   mEff <- readRDS(paste0(writePath, paste0("mEfftMax", tMaxi, "sMax", sMax,".rds")))
#
#   y <- mEff$qP$res0
#   dy <- mEff$qP$err0
#
#   r2F <- r2F(mEff$qP, r3, r1)
#
#   plotwitherror(x, y, dy, col = colors[1], xlab = "r", ylab = "mEff")
#
#   for(tMax in (tMaxi+1):tMaxf)
#   {
#     mEff <- readRDS(paste0(writePath, paste0("mEfftMax", tMax, "sMax", sMax,".rds")))
#
#     y <- mEff$qP$res0
#     dy <- mEff$qP$err0
#
#     plotwitherror(x, y, dy, col = colors[tMax-tMaxi+1], rep = TRUE)
#   }
#   legend(x = "bottomright",          # Position
#          legend = paste0("beta = ", betaArray),  # Legend texts
#          col = colors,           # Line colors
#          lwd = 2)                 # Line width
#   dev.off()
# }





# for(beta in 1.55)
# {
#   inputFileName <- inputFileName(spatialExtent
#                                  , temporalExtent
#                                  , invCoupling = beta
#                                  , sizeWLoops)
#
#   writePath <- writePath(inputFileName)
#   plotPath <- plotPath(inputFileName)
#
#   mEff <- readRDS(paste0(writePath, "mEfftMax15sMax8.rds"))
#
#   x <- c(1:8)
#
#   y1 <- mEff$wP$res0
#   dy1 <- apply(mEff$wP$rest, 2, sd)
#
#   y2 <- mEff$qP$res0
#   dy2 <- mEff$qP$err0
#
#   pdf(paste0(plotPath, "mEffPlot.pdf"))
#   plotwitherror(x, y1, dy1, col = "red")
#   plotwitherror(x, y2, dy2, rep = TRUE, col = "blue")
#   dev.off()
# }
