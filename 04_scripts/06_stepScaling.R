# step scaling function

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

betaMin <- 1
betaMax <- 10

sMax <- 8
tMax <- 16

r3 <- 3
r2 <- r3 - 1
r1 <- r2 - 1

g10 <- list()
g1Boot <- list()

g20 <- list()
g2Boot <- list()

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

  g1 <- r2F(mEff$qP, r1, r2)
  g2 <- r2F(mEff$qP, r2, r3)

  g10[[length(g10) + 1]] <- g1$g0
  g20[[length(g20) + 1]] <- g2$g0

  g1Boot[[length(g1Boot) + 1]] <- g1$gBoot
  g2Boot[[length(g2Boot) + 1]] <- g2$gBoot
}

y1 <- unlist(g10)
y2 <- unlist(g20)

y1Boot <- matrix(unlist(g1Boot), ncol = length(y1), byrow = TRUE)
y2Boot <- matrix(unlist(g2Boot), ncol = length(y2), byrow = TRUE)

dy1 <- unlist(lapply(g1Boot, sd))
dy2 <- unlist(lapply(g2Boot, sd))

fn <- function(par, x, boot.R, ...)
{
  par[2] / x + par[1]
}

fitResultg1 <- bootstrap.nlsfit(fn, c(0,1), y = y1[betaMin:betaMax], x = betaArray[betaMin:betaMax], y1Boot[,betaMin:betaMax], dy = dy1[betaMin:betaMax])

print(summary(fitResultg1))

pdf(paste0(plotPath, "fitg1.pdf"))
plot(fitResultg1)
dev.off()

fitResultg2 <- bootstrap.nlsfit(fn, c(0,1), y = y2[betaMin:betaMax], x = betaArray[betaMin:betaMax], y2Boot[,betaMin:betaMax], dy = dy2[betaMin:betaMax])

print(summary(fitResultg2))

pdf(paste0(plotPath, "fitg2.pdf"))
plot(fitResultg2)
dev.off()

pdf(paste0(plotPath, "gvsbeta.pdf"))

plotwitherror(betaArray, y1, dy1, col = "red", ylim = c(0, 1), xlab = "beta=1/g^2", ylab = "r^2F(r,g)")
plotwitherror(rep = TRUE, betaArray, y2, dy2, col = "blue", ylim = c(0, 1), xlab = "beta=1/g^2", ylab = "r^2F(r,g)")

legend(x = "topright",          # Position
legend = c(paste0("r1^2F(r1,g), r1 = ", r1), paste0("r2^2F(r2,g), r2 = ", r2)),  # Legend texts
col = c("red", "blue"),           # Line colors
lwd = 2)                 # Line width

dev.off()

newBeta <- fitResultg2$t0[2]/(y1[1] - fitResultg2$t0[1])
