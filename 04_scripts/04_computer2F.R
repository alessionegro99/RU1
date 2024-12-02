# Computing r2F

source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

distances <- c(1, 2, 3)

yLst <- rep(NA, length(spatialExtentArray))
bsyLst <- matrix(NA, nrow = bootSamples, ncol = length(spatialExtentArray))

V1vec <- rep(NA, length(spatialExtentArray))
bsV1vec <- matrix(NA, nrow = bootSamples, ncol = length(spatialExtentArray))

for (beta in betaArray){
  for(i in seq_along(spatialExtentArray)){
    spatialExtent <- spatialExtentArray[i]
      cat("Now L = ", spatialExtent, " beta = ", beta, "\n")

      V <- c(rep(NA,length(distances)))
      bsV <- matrix(data = NA, nrow = bootSamples, ncol = length(distances))

      for (r in seq_along(distances)){
        mEff <- readRDS(paste0("/home/negro/projects/stepscaling/RU1/02_output/data/heatbath_"
                              , spatialExtent
                              , "_"
                              , temporalExtent
                              ,"_"
                              , beta
                              ,"_0/plateauData_r"
                              , r,".rds"))
        V[r] <- mEff$effmassfit$t0[[1]]
        bsV[,r] <- mEff$effmassfit$t[,1]
      }

      Fr1 <- V[2] - V[1]
      Fr2 <- V[3] - V[2]

      g1 <- distances[1]^2*Fr1
      g2 <- distances[2]^2*Fr2

      bsg1 <- distances[1]^2*(bsV[,2]-bsV[,1])
      bsg2 <- distances[2]^2*(bsV[,3]-bsV[,2])

      cat(g1, sd(bsg1), "\n")
      cat(g2, sd(bsg2), "\n")
      cat("\n")

      yLst[i] <- g1
      yLst[i] <- g1
      bsyLst[, i] <- bsg1

      V1vec[i] <- V[1]
      bsV1vec[, i] <- bsV[, 1]
  }
  const <- function(par, x, bootSamples, ...) par[1] + 0*x

  # fitr2F <- bootstrap.nlsfit(const, c(1), x = spatialExtentArray, y = yLst, bsamples = bsyLst)
  # pdf(paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/stepScalingPlots/r2F_r_1_beta_", beta, ".pdf"))
  # plot(fitr2F, xlab = "L/a", ylab = "r²F(r,beta)", main = paste0("r²F(r,beta)=r²(V(r+1,beta)-V(r,beta)) for beta = ", beta, " r = 1"))
  # dev.off()


  pdf(paste0("/home/negro/projects/stepscaling/RU1/02_output/plots/stepScalingPlots/mEff_beta_", beta, ".pdf"))
  VvsL <- bootstrap.nlsfit(const, c(1), x = spatialExtentArray, y = V[i], bsamples = bsV[,i])
  plot(VvsL, xlab = "L/a", ylab = "mEff(L)", main = paste0("mEff(r,) for beta = ", beta, " r = 1"))

  dev.off()
}




