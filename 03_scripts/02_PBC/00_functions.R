setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fn <- function(par, x, boot.r, ...) par[1] * exp(-par[2] * x)

V_wt_ansatz <- function(par, x, boot.r, ...) par[1] + par[2]*x
V_Cornell <- function(par, x, boot.r, ...) par[1]*x + par[2]*1/x + par[3]

random_pastel_color <- function() {
  r <- runif(1, min = 0.5, max = 1)
  g <- runif(1, min = 0.5, max = 1)
  b <- runif(1, min = 0.5, max = 1)
  
  # Convert to hex
  rgb(r, g, b, maxColorValue = 1)
}

valid_indices <- function(vec, x) {
  foo <- which(vec < x)
  
  xi <- foo[1]
  xf <- foo[length(foo)-1]
  
  return(seq(xi:xf))
}