library(hadron)

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/OBC/"
output <- "/home/negro/projects/matching/RU1/02_output/restricted/OBC/"
output_general <- "/home/negro/projects/matching/RU1/02_output/general/OBC/"

## set refinement parameters
boot.R <- 200 # number of bootstrap samples (usually 200, 500 or 1000)

## set simulation parameters
tt <- 32 # temporal extent
SS <- c(3, 4, 5) # array of spatial extents to analyse
BB <- c(3, 4.25, 11.8) # array of inverse couplings to analyse
R0 <- 0 # starting point (OBC related)

RMAX <- c(4, 5, 6) # max length of Wloops

RRR123 <- list(c(1, sqrt(5), sqrt(8))
               , c(sqrt(2), sqrt(10), sqrt(18))
               , c(sqrt(5), sqrt(25), sqrt(32)))

RRR <- list() # available distances
for (ss in seq_along(SS)){
  foo <- c()
  for (i in seq(1, SS[ss] - 1)) {
    for (j in seq(0, i)) {
      value <- sqrt(i^2 + j^2)
      if (value < RMAX[ss]) {
        foo <- c(foo, value)
      }
    }
  }
  RRR[[length(RRR) + 1]] <- foo
}

III <- lapply(seq_along(RRR), function(i) match(RRR123[[i]], RRR[[i]]))

CC <- c("blue", "red", "orange")

pch_lst <- list(c(2, 6), c(1, 1), c(0, 8))
cex_lst <- list(c(1.5, 1.5), c(1.5, 0.5), c(1.5, 1.5))

sigma0 <- rep(NA, length(SS))
sigmat <- matrix(NA, nrow = boot.R, ncol = length(SS))

g1 <- rep(NA, length(SS))
dg1 <- rep(NA, length(SS))
bsg1 <- matrix(NA, nrow = boot.R, ncol = length(SS))

g2 <- rep(NA, length(SS))
dg2 <- rep(NA, length(SS))
bsg2 <- matrix(NA, nrow = boot.R, ncol = length(SS))

title_str <- ""
for (ss in SS){
  title_str <- paste0(title_str, "_L", ss)
}

for (ss in seq_along(SS)) {
  II <- III[[ss]]
  
  RR <- RRR[[ss]][II]
  
  V <- rep(NA, length(RR))
  bsV <- matrix(NA, nrow = boot.R, ncol = length(RR))
  
  bb <- BB[ss]
  folder <- paste0("pascal_OBC_", bb, "_", SS[ss], "_", tt, "_", R0)
  
  for (rr in seq_along(RR)) {
    tmp <- readRDS(paste0(output, folder, "/fit_results/fit.result_uncorrelated_", II[rr], ".rds"))
    
    V[rr] <- tmp$t0[[2]]
    bsV[, rr] <- tmp$t[, 2]
  }
  rm(tmp)
  
  g1 <- RR[1] * RR[1] * (V[2] - V[1]) / (RR[2] - RR[1])
  g2 <- RR[2] * RR[2] * (V[3] - V[2]) / (RR[3] - RR[2])

  bsg1 <- RR[1] * RR[1] * (bsV[, 2] - bsV[, 1]) / (RR[2] - RR[1])
  bsg2 <- RR[2] * RR[2] * (bsV[, 3] - bsV[, 2]) / (RR[3] - RR[2])
  
  print(g1)
  print(g2)
  print(g1/g2)
  sigma0[ss] <- g2/g1
  sigmat[, ss] <- bsg2/bsg1
  
  rm(V)
  rm(bsV)
  rm(g1)
  rm(g2)
  rm(bsg1, bsg2)
}

## constant fit
const <- function(par, x, boot.r,...) par[1] + 0*x

pdf(paste0(output_general, "CL/CLSS", title_str, "_b", BB[1], "_fit.pdf"))
fit.result <- bootstrap.nlsfit(fn = const, par.guess = c(1)
                               , y = sigma0, x = 1/(SS^2), bsamples = sigmat)
summary(fit.result)
par(xpd=FALSE, mar=c(4.5, 4.5, 4, 3) + 0.1)
plot(fit.result
     , pch = c(1,1,1), cex = 1.5
     , col = c("blue", "red", "orange")
     , main = expression(Sigma[s] == r[2]^2~F(r[2],g)/r[1]^2~F(r[1],g))
     , xlab = expression((a/L)^2)
     , ylab = expression(Sigma[s])
     , col.band = "lightblue"
     , col.line = rgb(150/255, 216/255, 230/255)  )
legend(x = "topleft"
       , legend = c(paste0("const = ", round(fit.result$t0[1],4))
                    , paste0("d_const = ", round(fit.result$se[1],4))
                    , paste0("chi2/dof = ", round(fit.result$chisqr/fit.result$dof,2))))
grid()
dev.off()