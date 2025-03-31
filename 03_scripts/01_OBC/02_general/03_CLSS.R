library(hadron)

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/OBC/"
output <- "/home/negro/projects/matching/RU1/02_output/restricted/OBC/"
output_general <- "/home/negro/projects/matching/RU1/02_output/general/OBC/"

## set refinement parameters
boot.R <- 200 # number of bootstrap samples (usually 200, 500 or 1000)

## set simulation parameters
tt <- 64 # temporal extent
SS <- c(8, 10, 12) # array of spatial extents to analyse
BB <- c(3, 3.75, 7.25) # array of inverse couplings to analyse
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

pdf(paste0(output_general, "CL/CL", title_str, "_b", BB[1], ".pdf"))
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
  
  g1[ss] <- RR[1] * RR[1] * (V[2] - V[1]) / (RR[2] - RR[1])
  g2[ss] <- RR[2] * RR[2] * (V[3] - V[2]) / (RR[3] - RR[2])
  
  bsg1[, ss] <- RR[1] * RR[1] * (bsV[, 2] - bsV[, 1]) / (RR[2] - RR[1])
  bsg2[, ss] <- RR[2] * RR[2] * (bsV[, 3] - bsV[, 2]) / (RR[3] - RR[2])
  
  rm(V)
  rm(bsV)
  
  dg1[ss] <- sd(bsg1[,ss])
  dg2[ss] <- sd(bsg2[,ss])
  
  ## plotting
  cc <- CC[ss]
  pch <- pch_lst[[ss]]
  cex <- cex_lst[[ss]]
  
  if (ss == 1) {
    plotwitherror(bb, g1[ss], dg1[ss],
                  col = cc, pch = pch[1], cex = cex[1],
                  xlab = "beta=1/g²", ylab = "r²F(r,g)",
                  xlim = c(3, 12), ylim = c(0.04, 0.23), xaxt = "n"
    )
    segments(x0 = bb, x1 = bb, y0 = g1[1], y1 = g2[ss], col = cc, lty = "dashed", lwd = 2)
    segments(x0 = bb, x1 = bb, y0 = 0, y1 = g1[ss], col = "lightgray", lty = "dotted", lwd = 2)
    axis(1, at = bb, labels = bb)
    grid()
  } else {
    plotwitherror(bb, g1[ss], dg1[ss],
                  col = cc, pch = pch[1], cex = cex[1],
                  rep = TRUE
    )
    segments(x0 = bb, x1 = bb
             , y0 = g1[ss], y1 = g2[ss], col = cc, lty = "dashed", lwd = 2)
    segments(x0 = bb, x1 = bb
             , y0 = 0, y1 = g1[ss]
             , col = "lightgray"
             , lty = "dotted"
             , lwd = 2)
    axis(1, at = bb, labels = bb)
  }
  plotwitherror(bb, g2[ss], dg2[ss],
                col = cc, pch = pch[2], cex = cex[2],
                rep = TRUE
  )
}

par(xpd = TRUE, mar = c(5, 4, 6, 2))

legend(
  x = "top",
  inset = c(0, -0.15),
  bty = "n",
  legend = c(
    bquote(r[1]^2 ~ F(r[1] == .(round(RRR[[1]][III[[1]]][1],2)), g) ~ .(SS[1])~"x"~.(SS[1])),
    bquote(r[2]^2 ~ F(r[2] == .(round(RRR[[1]][III[[1]]][2],2)), g) ~ .(SS[1])~"x"~.(SS[1])),
    bquote(r[1]^2 ~ F(r[1] == .(round(RRR[[2]][III[[2]]][1],2)), g) ~ .(SS[2])~"x"~.(SS[2])),
    bquote(r[2]^2 ~ F(r[2] == .(round(RRR[[2]][III[[2]]][2],2)), g) ~ .(SS[2])~"x"~.(SS[2])),
    bquote(r[1]^2 ~ F(r[1] == .(round(RRR[[3]][III[[3]]][1],2)), g) ~ .(SS[3])~"x"~.(SS[3])),
    bquote(r[2]^2 ~ F(r[2] == .(round(RRR[[3]][III[[3]]][2],2)), g) ~ .(SS[3])~"x"~.(SS[3]))
  ),
  pch = c(2, 6, 1, 1, 0, 8), # Point types
  pt.cex = c(1.5, 1.5, 1.5, 0.5, 1.5, 1.5), # Point sizes
  col = c("blue", "blue", "red", "red", "orange", "orange"), # Colors
  y.intersp = 1.5, # Vertical spacing
  ncol = 3
)
dev.off()

## constant fit
const <- function(par, x, boot.r,...) par[1] + 0*x

## linear fits for g1 and g2
fn <- function(par, x, boot.r,...) par[1] + par[2]*x

title_str <- ""
for (ss in SS){
  title_str <- paste0(title_str, "_L", ss)
}
pdf(paste0(output_general, "CL/CL", title_str, "_b", BB[1], "_fit.pdf"))
fit.result <- bootstrap.nlsfit(fn = fn, par.guess = c(1,1)
                               , y = g1, x = 1/(SS^2), bsamples = bsg1)
summary(fit.result)
par(xpd=FALSE, mar=c(4.5, 4.5, 4, 3) + 0.1)
plot(fit.result
     , pch = c(2,1,0), cex = 1.5
     , col = c("blue", "red", "orange")
     , xlab = expression(1/r[latt]^2)
     , ylab = expression(r[latt]^2 ~ F(r[latt], g))
     , col.band = "lightblue"
     , col.line = rgb(150/255, 216/255, 230/255)  )
legend(x = "topleft"
       , legend = c(paste0("a = ", round(fit.result$t0[1],4))
                    , paste0("da = ", round(fit.result$se[1],4))
                    , paste0("b = ", round(fit.result$t0[2],3))
                    , paste0("db = ", round(fit.result$se[2],4))))
grid()

fit.result <- bootstrap.nlsfit(fn = fn, par.guess = c(1,1)
                               , y = g2, x = 1/(SS^2), bsamples = bsg2)
summary(fit.result)
par(xpd=FALSE, mar=c(4.5, 4.5, 4, 3) + 0.1)
plot(fit.result
     , pch = c(2,1,0), cex = 1.5
     , col = c("blue", "red", "orange")
     , xlab = expression(1/r[latt]^2)
     , ylab = expression(r[latt]^2 ~ F(r[latt], g))
     , col.band = "lightblue"
     , col.line = rgb(150/255, 200/255, 215/255)  )
legend(x = "topleft"
       , legend = c(paste0("a = ", round(fit.result$t0[1],4))
                    , paste0("da = ", round(fit.result$se[1],4))
                    , paste0("b = ", round(fit.result$t0[2],3))
                    , paste0("db = ", round(fit.result$se[2],4))))
grid()

fit.result <- bootstrap.nlsfit(fn = const, par.guess = c(1)
                               , y = g2, x = 1/(SS^2), bsamples = bsg2)
summary(fit.result)
par(xpd=FALSE, mar=c(4.5, 4.5, 4, 3) + 0.1)
plot(fit.result
     , pch = c(2,1,0), cex = 1.5
     , col = c("blue", "red", "orange")
     , xlab = expression(1/r[latt]^2)
     , ylab = expression(r[latt]^2 ~ F(r[latt], g))
     , col.band = "lightblue"
     , col.line = rgb(150/255, 216/255, 230/255)  )
legend(x = "topleft"
       , legend = c(paste0("a = ", round(fit.result$t0[1],4))
                    , paste0("da = ", round(fit.result$se[1],4))))
grid()

dev.off()