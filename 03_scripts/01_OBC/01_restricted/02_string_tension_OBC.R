setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("00_functions.R")
library(hadron)

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/OBC/"
output <- "/home/negro/projects/matching/RU1/02_output/restricted/OBC/"

## set refinement parameters
boot.l <- 200 # block size
boot.R <- 200 # number of bootstrap samples (usually 200, 500 or 1000)
therm <- 500 # number of configuration to discard for thermalization

## set simulation parameters
TT <- c(64) # array of temporal extents to analyse
SS <- c(3) # array of spatial extents to analyse
BB <- c(3) # array of inverse couplings to analyse
RR <- seq(0,7) # number of replicas
R0 <- 0 # starting point (OBC related)

WT <- 10 # max T of Wloops 
WS <- 10 # max R of Wloops

non_planar <- TRUE

seq_max <- 10 # max plotting in wt

for (tt in TT) {
  for (ss in SS) {
    for (bb in BB) {
      
      r_i <- list()
      
      if(non_planar){
        for (i in seq(1, ss - 1)) {
          for (j in seq(0, i)) {
            if ((sqrt(i^2 + j^2)) < WS) {
              r_i[[length(r_i) + 1]] <- sqrt(i^2 + j^2)
            }
          }
        }
      } else {
        r_i <- seq(1, WS)
      }
      
      r_i <- unlist(r_i)
      
      random_pastel_colors <- c()
      
      for(c in seq(1:max(WT, WS))){
        random_pastel_colors <- c(random_pastel_colors, random_pastel_color())
      }
      
      data <- NULL
      
      for(rr in RR){
        folder <- paste0("pascal_OBC_", bb, "_", ss, "_", tt, "_", R0,"_",rr)
        tmp <- as.matrix(read.table(paste0(rawdata, folder,".dat"), skip = 27 + therm))
        if (is.null(data)) data <- tmp else data <- rbind(tmp, data)
      }
      folder <- paste0("pascal_OBC_", bb, "_", ss, "_", tt, "_", R0)
      
      if (!dir.exists(paste0(output, folder))) {
        dir.create(paste0(output, folder))
      }
      
      message(sprintf(
        "Analyzing data for:\n- Temporal extent: %s\n- Spatial extent: %s\n- Inverse coupling: %s",
        tt, ss, bb
      ))
      
      message("thermalization...")
      ## MC history plot for every distance in r_i
      pdf(paste0(output, folder, "/thermalization.pdf"))
      for (i in seq_along(r_i)) {
        W <- data[, seq(4 + i, ncol(data), length(r_i))]
        for (j in c(1)) {
          time_series <- W[, j]
          plot(
            x = seq(1,length(time_series), by = 500),
            y = time_series[seq(1, length(time_series), by = 500)],
            type = "l",  # Line plot instead of scatter
            xlab = "n_configs",
            ylab = "W(ws,wt)",
            main = paste0("W(ws,wt) for T = ", j, ", R = ", round(r_i[i],2))
          )
        }
      }
      rm(W)
      rm(time_series)
      dev.off()
      
      message("bootstrap analysis...")
      ## bootstrap analysis for every distance in r_i
      pdf(paste0(output, folder, "/bootstrap_analysis.pdf"))
      
      for (i in seq_along(r_i)) {
        tmp <- data[, seq(4 + i, ncol(data), length(r_i))]
        Time <- ncol(tmp)
        
        tmp <- cf_orig(cf = tmp)
        tmp <- cf_meta(tmp, nrObs = 1, Time = Time)
        
        for (j in c(1)){
          suppressMessages(bootstrap.analysis(tmp$cf[, j], boot.R = boot.R, boot.l = 2, pl = TRUE))
        }
      }
      rm(tmp)
      dev.off()
      
      message("bootstrapping Wilson loops data...")
      ## bootstrapping for every time
      W <- list()
      
      for (i in seq(1, WS)) {
        tmp <- data[, seq(4 + 1 + (i-1)*length(r_i), 4 + i*length(r_i))]
        Time <- ncol(tmp)
        
        tmp <- cf_orig(cf = tmp)
        tmp <- cf_meta(tmp, nrObs = 1, nrStypes = 0, Time = Time)
        
        W[[length(W) + 1]] <- bootstrap.cf(cf = tmp, boot.R = boot.R, boot.l = boot.l, sim = "fixed", endcorr = TRUE, seed = 1234567)
      }
      
      # plotting V vs wr for every wt ##########################################
      
      V_wt <- list()
      bs_V_wt <- list()
      
      for(wt in seq(1,seq_max)){
        V_wt[[wt]] <- -1/wt*log(W[[wt]]$cf0)
        bs_V_wt[[wt]] <- -1/wt*log(W[[wt]]$cf.tsboot$t)
      }
      
      legend <- c()
      
      plotwitherror(r_i, V_wt[[1]], apply(bs_V_wt[[1]], 2, sd, na.rm = TRUE), col = random_pastel_colors[1], pch=1, xlab ="wr", ylab ="V(wr)", xlim = c(0, max(r_i)), ylim = c(0,max(V_wt[[1]])))
      legend <- c(legend, paste0("wt = ", 1))
      
      for(i in seq(2,seq_max)){
        plotwitherror(r_i, V_wt[[i]], apply(bs_V_wt[[i]], 2, sd, na.rm = TRUE), col = random_pastel_colors[i], rep = TRUE, pch=1)
        legend <- c(legend, paste0("wt = ", i))
      }
      grid()
      legend("topleft", legend = legend, col = random_pastel_colors[1:length(legend)], pch=1)
      dev.off()
      
      ##########################################################################
    }
  }
}