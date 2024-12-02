source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

# plotting the MC history for Polyakov loops

arrayBeta <- seq(0.5, 6.0, 0.1)

skip <- 10

polyModule <- matrix(NA, nrow=1000-skip, ncol=length(arrayBeta))
polyModulesq <- matrix(NA, nrow=1000-skip, ncol=length(arrayBeta))

for(index in seq_along(arrayBeta)){
  i <- sprintf("%.1f", arrayBeta[index])
  dataFile <- paste0("/home/negro/projects/stepscaling/RU1/01_rawdata/metropolis/metropolis_3_16_", i, "_0/omeas/spatial_polyakov/polyakov.dat")
  tmp <- as.matrix(read.table(dataFile))
  tmp1 <- as.numeric(tmp[(2+skip):nrow(tmp), 2])
  tmp2 <- as.numeric(tmp[(2+skip):nrow(tmp), 3])

  polyModule[, index] <- sqrt(tmp1^2+tmp2^2)
  polyModulesq[, index] <- tmp1^2+tmp2^2
}

polyModulecf <- cf_orig(cf=polyModule)
polyModulecf <- cf_meta(polyModulecf, nrObs=1, Time=length(arrayBeta))

polyModulecf <- bootstrap.cf(polyModulecf, boot.R = 500, boot.l = 10, seed=1442556)

polyModulesqcf <- cf_orig(cf=polyModulesq)
polyModulesqcf <- cf_meta(polyModulesqcf, nrObs=1, Time=length(arrayBeta))

polyModulesqcf <- bootstrap.cf(polyModulesqcf, boot.R = 500, boot.l = 10, seed=1442556)

chiPt0 <- polyModulesqcf$cf.tsboot$t0 - polyModulecf$cf.tsboot$t0^2

chiPt <- polyModulesqcf$cf.tsboot$t - polyModulecf$cf.tsboot$t^2

errchiPt0 <- apply(chiPt, 2, sd)

x <- arrayBeta
y <- chiPt0
dy <- errchiPt0

pdf("/home/negro/projects/stepscaling/RU1/02_output/plots/pLoop/chiPloop.pdf")

plotwitherror <- function(x, y, dy) {
  plot(x=x, y=y, xlab = "beta", ylab = "<|L|²>-<|L|>²", main = "Polyakov loop susceptibility")
  arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01,
         angle=90, code=3)
}
plotwitherror(x=x, y=y, dy=dy)
dev.off()

pdf("/home/negro/projects/stepscaling/RU1/02_output/plots/pLoop/Ploop.pdf")
plot(polyModulecf, xlab = "beta", ylab = "<|L|>", main = "Polyakov loop", xaxt="n")
axis(1, at = seq_along(arrayBeta) - 1,  labels = arrayBeta)
dev.off()

## bootstrapped fit

initfit <- 19
endfit <- 31

fn <- function (par, x, boot.r, ...) par[1]^2/(par[2]^2+(x-par[3])^2)

value <- y[initfit:endfit]
dvalue <- dy[initfit:endfit]
xval <- x[initfit:endfit]
boot.R <- 500

bsamples <- chiPt[,initfit:endfit]

fit.result <- bootstrap.nlsfit(fn, c(1, 1, 3), value, xval, bsamples)
summary(fit.result)
plot(fit.result, xlim = c(0.5,6), ylim = c(0.00005, 0.00016))
points(x, y, pch = 10)
arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01, angle=90, code=3)


