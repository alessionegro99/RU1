library(hadron)

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/OBC/"
plots <- "/home/negro/projects/matching/RU1/02_output/plots/restricted/OBC/"
plots_gen <- "/home/negro/projects/matching/RU1/02_output/plots/general/OBC/"
datas <- "/home/negro/projects/matching/RU1/02_output/data/OBC/"

## set refinement parameters
boot.l <- 500 # block size
boot.R <- 500 # number of bootstrap samples (usually 200, 500 or 1000)
therm <- 500 # number of configuration to discard for thermalization

## set simulation parameters
tt <- 32 # temporal extent to analyse
SS <- c(seq(5, 20, by = 1), seq(22, 32, by  = 2)) # array of spatial extents to analyse
bb <- 11.8 # array of inverse couplings to analyse
R0 <- 0 # starting point (OBC related)

CC <- c("purple", "orange")

RMAX <- 6 # max length of Wloops

RR123 <- c(sqrt(5), sqrt(25), sqrt(32))
RRR <- list() # available distances
for (ss in seq_along(SS)){
  foo <- c()
  for (i in seq(1, SS[ss] - 1)) {
    for (j in seq(0, i)) {
      value <- sqrt(i^2 + j^2)
      if (value < RMAX) {
        foo <- c(foo, value)
      }
    }
  }
  RRR[[length(RRR) + 1]] <- foo
}

III <- lapply(seq_along(RRR), function(i) match(RR123, RRR[[i]]))

pdf(paste0(plots_gen
           , "FSE/FSE_r2F_"
           , round(RR123[1],2)
           , "_", round(RR123[2],2)
           , "_", round(RR123[3],2)
           , "_", bb
           , ".pdf"))
for (ss in seq_along(SS))
{
  II <- III[[ss]]
  
  RR <- RRR[[ss]][II]
  
  V <- rep(NA, length(RR))
  bsV <- matrix(NA, nrow = boot.R, ncol = length(RR))
  
  folder <- paste0("pascal_OBC_"
                   , bb
                   , "_", SS[ss]
                   , "_", tt
                   , "_", R0)
  
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
  
  if(ss == 1){
    par(xpd = TRUE, mar = c(4.5, 4.5, 4, 3))
    plotwitherror(SS[ss], g1, sd(bsg1)
                  , col = CC[1], pch = 2, cex = 1.5
                  , ylab = expression(r^2*F(r,g)), xlab = "L/a"
                  , xlim = c(SS[1], SS[length(SS)])
                  , ylim = c(g1*0.85, g2))
    par(xpd = FALSE)
    grid()
    par(xpd = TRUE)
  }
  else{
    plotwitherror(SS[ss], g1, sd(bsg1)
                  , col = CC[1], pch = 2, cex = 1.5
                  , rep = TRUE)
  }

  legend(
    x = "topright",
    legend = c(bquote(r[1] == .(round(RR[1],2)))
               , bquote(r[2] == .(round(RR[2],2)))
               , bquote(r[3] == .(round(RR[3],2)))
               , bquote("beta" == .(bb))),
    ncol = 2
  )
  legend(
    x = "top",
    inset = c(0, -0.1),
    bty = "n",
    legend = c(
      expression(r^2*F(r[1],g)),
      expression(r^2*F(r[2],g))
    ),
    pch = c(2, 6, 0), 
    pt.cex = 1.5, 
    col = CC,
    
    y.intersp = 1.5, 
    ncol = 2
  )
  plotwitherror(SS[ss], g2, sd(bsg2)
                , col = CC[2], pch = 6, cex = 1.5, rep = TRUE)
}
dev.off()


