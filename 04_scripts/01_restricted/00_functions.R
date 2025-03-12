setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

fn <- function(par, x, boot.r, ...) par[1] * exp(-par[2] * x)

valid_indices <- function(vec, x) {
  foo <- which(vec < x)
  
  xi <- foo[1]
  xf <- foo[length(foo)-1]
  
  return(seq(xi:xf))
}