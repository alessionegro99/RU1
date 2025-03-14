value <- c(0.0843, 0.0972, 0.1147)
x <- c(2.1, 2.0, 1.9)
dvalue <- c(0.0006, 0.0005, 0.0004)
boot.R <- 1500

#dx <- c(0, 0, 0)
dx <- c(0.01, 0.01, 0.01)
bsamples <- parametric.bootstrap(boot.R, c(value, x), c(dvalue, dx))


#bsamples <- array(sample(x=value, size=length(value)*boot.R, replace=TRUE), dim=c(boot.R,length(value)))


fn <- function (par, x, boot.r, ...) par[1] + par[2] / x ^2


fit.result <- bootstrap.nlsfit(fn, c(-0.05, 0.6), value, x, bsamples, use.minpack.lm = TRUE)

summary(fit.result)
plot(fit.result, main = 'Ribbon on top')

a <- fit.result$t0[1]
b <- fit.result$t0[2]

cat(sqrt(b/(-a + 0.0979)))

# source("~/projects/stepscaling/RU1/04_scripts/01_preliminaries.R")
# source("~/projects/stepscaling/RU1/04_scripts/02_Wrt.R")
# source("~/projects/stepscaling/RU1/04_scripts/03_computemeffWrt.R")
