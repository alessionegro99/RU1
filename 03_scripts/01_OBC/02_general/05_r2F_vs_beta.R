library(hadron)

plots_gen <- "/home/negro/projects/matching/RU1/02_output/plots/general/OBC/"
datas <- "/home/negro/projects/matching/RU1/02_output/data/OBC/"

## set refinement parameters
boot.R <- 200 # number of bootstrap samples (usually 200, 500 or 1000)

## set simulation parameters
tt <- 64 # temporal extent to analyse
ss <- 7
BB <- c(11, 12, 13, 14, 15, 16, 17, 18) # array of inverse couplings to analyse
R0 <- 0 # starting point (OBC related)

y1 <- rep(NA, length(BB))
dy1 <- rep(NA, length(BB))

y2 <- rep(NA, length(BB))
dy2 <- rep(NA, length(BB))

RMAX <- 9 # max length of Wloops

RR123 <- c(sqrt(8), sqrt(40), sqrt(72))
RR <- c() # available distances

for (i in seq(1, ss - 1)) {
  for (j in seq(0, i)) {
    value <- sqrt(i^2 + j^2)
    if (value < RMAX) {
      RR <- c(RR, value)
    }
  }
}

II <- match(RR123, RR)


for (bb in seq_along(BB))
{
  RR <- RR[II]
  V <- rep(NA, length(RR))
  bsV <- matrix(NA, nrow = boot.R, ncol = length(RR))

  folder <- paste0(
    "pascal_OBC_",
    bb,
    "_", ss,
    "_", tt,
    "_", R0
  )

  for (rr in seq_along(RR)) {
    tmp <- readRDS(paste0(datas, folder, "/fit.result_uncorrelated_", II[rr], ".rds"))

    V[rr] <- tmp$t0[[2]]
    bsV[, rr] <- tmp$t[, 2]
  }
  rm(tmp)

  g1 <- RR[1] * RR[1] * (V[2] - V[1]) / (RR[2] - RR[1])
  g2 <- RR[2] * RR[2] * (V[3] - V[2]) / (RR[3] - RR[2])

  bsg1 <- RR[1] * RR[1] * (bsV[, 2] - bsV[, 1]) / (RR[2] - RR[1])
  bsg2 <- RR[2] * RR[2] * (bsV[, 3] - bsV[, 2]) / (RR[3] - RR[2])

  rm(V)
  rm(bsV)

  y1 <- c(y1, g1)
  y2 <- c(y2, g2)
  dy1 <- c(dy1, sd(bsg1))
  dy2 <- c(dy2, sd(bsg2))
}

pdf(paste0(
  plots_gen,
  "g1_g2_vs_beta/g1_g2_vs_beta_L",
  ss,
  "_", round(RR123[1], 2),
  "_", round(RR123[2], 2),
  "_", round(RR123[3], 2),
  ".pdf"
))

plotwitherror(BB, y1, dy1, col = "darkblue"
              , pch = 2, cex = 1.5
              , ylab = expression(r[1]^2 * F(r[1], g)), xlab = expression(beta)
              , main = pasteo0("L = ", ss, " ,r1 = ", RR123[1], " ,r2 = ", RR123[2]))
grid()

plotwitherror(BB, y2, dy2, col = "darkblue"
              , pch = 6, cex = 1.5
              , ylab = expression(r[2]^2 * F(r[2], g)), xlab = expression(beta)
              , main = pasteo0("L = ", ss, " ,r1 = ", RR123[1], " ,r2 = ", RR123[2]))
grid()

dev.off()