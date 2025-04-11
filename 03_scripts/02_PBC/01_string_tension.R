setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("00_functions.R")

library(hadron)

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/PBC/"
output <- "/home/negro/projects/matching/RU1/02_output/restricted/PBC/"

## set refinement parameters
boot.l <- 500 # block size
boot.R <- 200 # number of bootstrap samples (usually 200, 500 or 1000)
therm <- 1000 # number of configuration to discard for thermalization

## set simulation parameters
TT <- c(42) # array of temporal extents to analyse
SS <- c(42) # array of spatial extents to analyse
BB <- c(1.4) # array of inverse couplings to analyse
RR <- seq(0, 7) # number of replicas

WT <- 10 # max T of Wloops
WS <- 10 # max R of Wloops

non_planar <- FALSE

seq_min <- 2 # min fitting in wt
seq_max <- 3 # max plotting in wt

random_pastel_colors <- c()

for (i in seq_along(1:max(WT, WS))) {
  random_pastel_colors <- c(random_pastel_colors, random_pastel_color())
}

for (tt in TT) {
  for (ss in SS) {
    for (bb in BB) {
      r_i <- list()

      if (non_planar) {
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

      data <- NULL

      for (rr in RR) {
        folder <- paste0("pascal_PBC_", bb, "_", ss, "_", tt, "_", rr)
        tmp <- as.matrix(read.table(paste0(rawdata, folder, ".dat"), skip = 27 + therm))
        if (is.null(data)) data <- tmp else data <- rbind(tmp, data)
      }
      folder <- paste0("pascal_PBC_", bb, "_", ss, "_", tt)

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
            x = seq(1, length(time_series), by = 250),
            y = time_series[seq(1, length(time_series), by = 250)],
            type = "l", # Line plot instead of scatter
            xlab = "n_configs",
            ylab = "W(ws,wt)",
            main = paste0("W(ws,wt) for T = ", j, ", R = ", round(r_i[i], 2))
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

        for (j in c(1)) {
          suppressMessages(bootstrap.analysis(tmp$cf[, j], boot.R = boot.R, boot.l = 2, pl = TRUE))
        }
      }
      rm(tmp)
      dev.off()

      message("bootstrapping Wilson loops data...")
      ## bootstrapping for every time
      W <- list()

      for (i in seq_along(r_i)) {
        tmp <- data[, seq(4 + 1 + (i - 1) * length(r_i), 4 + i * length(r_i))]
        Time <- ncol(tmp)

        tmp <- cf_orig(cf = tmp)
        tmp <- cf_meta(tmp, nrObs = 1, nrStypes = 0, Time = Time)

        W[[length(W) + 1]] <- bootstrap.cf(cf = tmp, boot.R = boot.R, boot.l = boot.l, sim = "fixed", endcorr = TRUE, seed = 1234567)
      }

      # plotting V vs wr for every wt ##########################################

      V_wt <- list()
      bs_V_wt <- list()

      for (wt in seq_along(1:seq_max)) {
        V_wt[[wt]] <- -1 / wt * log(W[[wt]]$cf0)
        bs_V_wt[[wt]] <- -1 / wt * log(W[[wt]]$cf.tsboot$t)
      }

      legend <- c()

      pdf(paste0(output, folder, "/W_wt(wr).pdf"))

      for (i in seq_along(1:seq_max)) {
        if (i == 1) {
          plotwitherror(r_i, V_wt[[i]], apply(bs_V_wt[[i]], 2, sd, na.rm = TRUE), col = random_pastel_colors[i], pch = 1, xlab = "wr", ylab = "V(wr)", xlim = c(0, max(r_i)), ylim = c(0, max(V_wt[[i]])))
        } else {
          plotwitherror(r_i, V_wt[[i]], apply(bs_V_wt[[i]], 2, sd, na.rm = TRUE), col = random_pastel_colors[i], pch = 1, rep = TRUE)
        }
        legend <- c(legend, paste0("wt = ", i))
      }
      grid()
      legend("topleft", legend = legend, col = random_pastel_colors[1:length(legend)], pch = 1)
      dev.off()

      ##########################################################################

      # plotting V vs wt for every wr ##########################################

      V_wr <- list()
      bs_V_wr <- list()

      for (i in seq_along(r_i)) {
        tmp <- c()
        bs_tmp <- list()
        for (j in seq(1, seq_max)) {
          tmp <- c(tmp, V_wt[[j]][i])
          bs_tmp[[j]] <- bs_V_wt[[j]][, i]
        }
        V_wr[[i]] <- tmp
        bs_V_wr[[i]] <- bs_tmp
      }

      mask <- rep(list(seq(seq_min, seq_max)), length(r_i))

      legend <- c()

      pdf(paste0(output, folder, "/W_wr(wt).pdf"))

      for (i in seq_along(r_i)) {
        if (i == 1) {
          plotwitherror(1 / seq(1:seq_max), V_wr[[i]],
            apply(matrix(unlist(bs_V_wr[[i]]), nrow = boot.R, ncol = seq_max), 2, sd),
            col = random_pastel_colors[i], pch = 1, xlab = "1/wt", ylab = "V(wt)", ylim = c(min(unlist(V_wr)), max(unlist(V_wr))), xlim = c(0, 1)
          )
        } else {
          plotwitherror(1 / seq(1, seq_max), V_wr[[i]], apply(matrix(unlist(bs_V_wr[[i]]), nrow = boot.R, ncol = seq_max), 2, sd), col = random_pastel_colors[i], pch = 1, rep = TRUE)
        }

        legend <- c(legend, paste0("wr = ", i))
      }
      grid()
      legend("topleft", legend = legend, col = random_pastel_colors[1:length(legend)], pch = 1)

      for (i in seq_along(r_i)) {
        if (i == 1) {
          fit.result <- bootstrap.nlsfit(fn = V_wt_ansatz, par.guess = c(1, 1), x = 1 / seq(1, seq_max), y = V_wr[[1]], bsamples = matrix(unlist(bs_V_wr[[1]]), nrow = boot.R, ncol = seq_max), mask = mask[[1]])
          cat(c(fit.result$t0, fit.result$se, fit.result$chisqr / fit.result$dof), file = paste0(output, folder, "/fit.csv"), fill = TRUE, append = FALSE)
          cat(c(fit.result$t[, 1], "\n"), file = paste0(output, folder, "/fit_boot.csv"), append = FALSE, fill = FALSE)
          cat(c(fit.result$t[, 2], "\n"), file = paste0(output, folder, "/fit_boot.csv"), append = TRUE, fill = FALSE)
          plot(fit.result, col = random_pastel_colors[i], pch = 1, xlab = "1/wt", ylab = "V(wt)", ylim = c(min(unlist(V_wr)), max(unlist(V_wr))), xlim = c(0, 1), col.line = random_pastel_colors[i], col.band = random_pastel_colors[i], plot.range = c(0, 1))
        } else {
          fit.result <- bootstrap.nlsfit(fn = V_wt_ansatz, par.guess = c(1, 1), x = 1 / seq(1, seq_max), y = V_wr[[i]], bsamples = matrix(unlist(bs_V_wr[[i]]), nrow = boot.R, ncol = seq_max), mask = mask[[i]])
          plot(fit.result, col = random_pastel_colors[i], pch = 1, col.line = random_pastel_colors[i], col.band = random_pastel_colors[i], plot.range = c(0, 1), rep = TRUE)
          cat(c(fit.result$t0, fit.result$se, fit.result$chisqr / fit.result$dof), file = paste0(output, folder, "/fit.csv"), append = TRUE, fill = TRUE)
          cat(c(fit.result$t[, 1], "\n"), file = paste0(output, folder, "/fit_boot.csv"), append = TRUE, fill = FALSE)
          cat(c(fit.result$t[, 2], "\n"), file = paste0(output, folder, "/fit_boot.csv"), append = TRUE, fill = FALSE)
        }
      }
      dev.off()

      ##########################################################################

      data <- as.matrix(read.table(paste0(output, folder, "/fit.csv"), header = FALSE, sep = " "))
      bs_data <- t(as.matrix(read.table(paste0(output, folder, "/fit_boot.csv"))))

      y <- data[, 1]
      bs_y <- bs_data[, seq(1, ncol(bs_data), 2)]

      fit.result_Cornell <- bootstrap.nlsfit(fn = V_Cornell, par.guess = c(0.1, 1, 1), x = r_i, y = y, bsamples = bs_y, mask = seq(3, 9))

      pdf(paste0(output, folder, "/V_Cornell.pdf"))
      plot(fit.result_Cornell, col = random_pastel_colors[1], pch = 1, col.line = random_pastel_colors[1], col.band = random_pastel_colors[1], plot.range = c(1, max(r_i)), xlab = "wr", ylab = "V_wr(wt=inf)")
      grid()
      legend("topleft", legend = c(paste0("sigma = ", round(fit.result_Cornell$t0[1], 4)), paste0("d_sigma = ", round(fit.result_Cornell$se[1], 4)), paste0("chi2/dof = ", round(fit.result_Cornell$chisqr / fit.result_Cornell$dof, 2))))
      dev.off()
    }
  }
}
