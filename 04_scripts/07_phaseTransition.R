source("~/projects/stepscaling/RU1/03_functions/header.R")
source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

x <- c(4, 6, 8)

y <- c(1.83, 2.08, 2.14)
dy <- c(0.02, 0.02, 0.05)

weights <- 1 / (dy^2)

fit <- nls(y ~ ( a + b * log(c * x)), start = list(a = 0, b = 0.5, c = 1), weights = weights)
summary(fit)

# Plot data points
plot(x, y, main = "Weighted Nonlinear Fit with Error Bars", xlab = "Predictor (x)", ylab = "Response (y)", pch = 16, ylim = c(1.5,2.5), xlim = c(4, 16))

# Add error bars
arrows(x, y - dy, x, y + dy, angle = 90, code = 3, length = 0.1, col = "gray")

# Add the fitted curve
curve(coef(fit)["a"] + coef(fit)["b"] * log(coef(fit)["c"]*x), add = TRUE, col = "blue")

