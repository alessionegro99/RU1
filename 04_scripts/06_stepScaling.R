# step scaling function

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

sMax <- 8
tMax <- 16

r3 <- 3
r2 <- r3 - 1
r1 <- r2 - 1

g20 <- list()
g2Boot <- list()

g30 <- list()
g3Boot <- list()

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

  mEff <- readRDS(paste0(writePath, paste0("mEfftMax", tMax, "sMax", sMax,".rds")))

  g2 <- r2F(mEff$qP, r2, r1)
  g3 <- r2F(mEff$qP, r3, r2)

  g20[[length(g20) + 1]] <- g2$g0
  g30[[length(g30) + 1]] <- g3$g0

  g2Boot[[length(g2Boot) + 1]] <- g2$gBoot
  g3Boot[[length(g3Boot) + 1]] <- g3$gBoot
}

y1 <- unlist(g20)
y2 <- unlist(g30)
dy1 <- unlist(lapply(g2Boot, sd))
dy2 <- unlist(lapply(g3Boot, sd))

pdf(paste0(plotPath, "gvsbeta.pdf"))

plotwitherror(betaArray, y = y2, dy = dy2, ylim = c(0.05, 0.30), xlab = "beta", ylab = "g", col = "red")
plotwitherror(betaArray, y1, dy1, rep = TRUE, ylim = c(0.05, 0.30), col = "blue")

legend(x = "topright",          # Position
legend = paste0("r = ", c(3, 2)),  # Legend texts
col = c("red", "blue"),           # Line colors
lwd = 2)                 # Line width

dev.off()
