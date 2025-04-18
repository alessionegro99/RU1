setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fn <- function(par, x, boot.r, ...) par[1] * exp(-par[2] * x)

V_wt_ansatz <- function(par, x, boot.r, ...) par[1] + par[2]*x
V_Cornell <- function(par, x, boot.r, ...) par[1]*x + par[2]*1/x + par[3]

random_pastel_color <- function() {
  # Generate RGB values closer to white (pastel)
  r <- runif(1, min = 0, max = 0.6)
  g <- runif(1, min = 0, max = 0.6)
  b <- runif(1, min = 0, max = 0.6)
  
  # Convert to hex
  rgb(r, g, b, maxColorValue = 1)
}

valid_indices <- function(vec, x) {
  foo <- which(vec < x)
  
  xi <- foo[1]
  xf <- foo[length(foo)-1]
  
  return(seq(xi:xf))
}