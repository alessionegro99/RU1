setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
getwd()

library(hadron)

source("00_functions.R")

rawdata <- "/home/negro/projects/matching/RU1/01_rawdata/OBC/"
plots <- "/home/negro/projects/matching/RU1/02_output/plots/restricted/OBC/"
datas <- "/home/negro/projects/matching/RU1/02_output/data/OBC/"

## set refinement parameters
boot.l <- 500 # block size
boot.R <- 500 # number of bootstrap samples (usually 200, 500 or 1000)
therm <- 500 # number of configuration to discard for thermalization

## set simulation parameters
TT <- c(32) # array of temporal extents to analyse
SS <- c(19, 20, 22, 24, 26, 28, 30, 32) # array of spatial extents to analyse
BB <- c(11.8) # array of inverse couplings to analyse
R0 <- 0 # starting point (OBC related)

RMAX <- 6 # max length of Wloops

## array of distances (simulation dependent)

for (temporal_extent in TT) {
  for (spatial_extent in SS) {
    for (inv_coupling in BB) {
      r_i <- list()
      for(i in seq(1,spatial_extent-1)){
        for(j in seq(0,i)){
          if((sqrt(i^2 + j^2)) < RMAX)
            r_i[[length(r_i) + 1]] <- sqrt(i^2 + j^2)
        }
      }

      #######################################
      folder <- paste0("pascal_OBC_", inv_coupling, "_", spatial_extent, "_", temporal_extent, "_", R0)
      if (!dir.exists(paste0(plots, folder))) {
        dir.create(paste0(plots, folder))
      }
      if (!dir.exists(paste0(datas, folder))) {
        dir.create(paste0(datas, folder))
      }
      
      message(sprintf(
        "Analyzing data for:\n- Temporal extent: %s\n- Spatial extent: %s\n- Inverse coupling: %s",
        temporal_extent, spatial_extent, inv_coupling
      ))

      ## skipping the header
      data <- as.matrix(read.table(paste0(rawdata, folder, "/", folder, ".dat"), skip = 20 ))
      #######################################
      
      message("thermalization...")
      ## MC history plot for every distance in r_i
      pdf(paste0(plots, folder, "/thermalization.pdf"))
      for (i in seq_along(r_i)) {
        W <- data[, seq(4 + i, ncol(data), length(r_i))]
        time_series <- W[, 1]
        plot(x = c(1:length(time_series)), y = time_series, xlab = "n_configs", ylab = "W(R,T)", main = paste0("W(R,T) for T = 1, R = ", r_i[i]))

        rm(W)
        rm(time_series)
      }
      dev.off()

      data <- data[-(1:therm),]
      #######################################
      
      message("bootstrap analysis...")
      ## bootstrap analysis for every distance in r_i
      pdf(paste0(plots, folder, "/bootstrap_analysis.pdf"))
      for (i in seq_along(r_i)) {
        tmp <- data[, seq(5 + (i - 1), ncol(data), length(r_i))]
        Time <- ncol(tmp)

        tmp <- cf_orig(cf = tmp)
        tmp <- cf_meta(tmp, nrObs = 1, Time = Time)

        suppressMessages(bootstrap.analysis(tmp$cf[, 1], boot.R = boot.R, boot.l = 2, pl = TRUE))
        rm(tmp)
      }
      dev.off()
      #######################################
      
      message("bootstrapping Wilson loops data...")
      ## bootstrapping for every distance in r_i
      W <- list()
      for (i in seq_along(r_i)) {
        tmp <- data[, seq(5 + (i - 1), ncol(data), length(r_i))]
        Time <- ncol(tmp)

        tmp <- cf_orig(cf = tmp)
        tmp <- cf_meta(tmp, nrObs = 1, nrStypes = 0, Time = Time)

        W[[length(W) + 1]] <- bootstrap.cf(cf = tmp, boot.R = boot.R, boot.l = boot.l, sim = "fixed", endcorr = TRUE, seed = 1234567)
      }
      #######################################
      
      message("saving bootstrapped Wilson loops data...")
      ## saving wloop data
      saveRDS(W, paste0(datas, folder, "/Wloops.rds"))
      #######################################

      message("plotting Wilson loops...")
      ## plotting wloop data
      pdf(paste0(plots, folder, "/NPWL.pdf"))
      for (i in seq_along(r_i)) {
        plot(W[[i]], main = paste0("Non planar Wilson loop for R = ", r_i[i]), ylab = "W(R,T)", xlab = "t")
        grid()
        plot(W[[i]], log = "y", main = paste0("Non planar Wilson loop for R = ", r_i[i], ", logscale"), ylab = "W(R,T)", xlab = "t")
        grid()
      }
      dev.off()
      #######################################
      
      message("computing effective mass...")
      ## computing effective mass
      if (!exists("W")) {
        W <- readRDS(paste0(datas, folder, "/Wloops.rds"))
      }
      pdf(paste0(plots, folder, "/EMASS.pdf"))
      for (i in seq_along(r_i)) {
        e_mass <- bootstrap.effectivemass(W[[i]], type = "log")
        ulim <- e_mass$effMass[1] + e_mass$deffMass[1]
        llim <- e_mass$effMass[8] - 3*e_mass$deffMass[8]
        lims <- c(llim, ulim)
        plot(e_mass, main = paste0("Effective mass for r = ", r_i[i]), ylab = expression(m[eff]), xlab = "t", xlim = c(0, 16), ylim = lims)
        grid()
      }
      dev.off()
      #######################################
      
      message("performing uncorrelated fit to the data...")
      ## performing uncorrelated fit to Wloop in a range specified by mask
      if (!exists("W")) {
        W <- readRDS(paste0(datas, folder, "/Wloops.rds"))
      }
      
      mask <- rep(list(seq(6, 16)), length(r_i))
                  
      range <- seq(1, temporal_extent - 1)

      pdf(paste0(plots, folder, "/uncorrelated_fit.pdf"))
      for (i in seq_along(r_i)) {
        value <- W[[i]]$cf.tsboot$t0
        bsamples <- W[[i]]$cf.tsboot$t

        fit.result_uncorrelated <- bootstrap.nlsfit(fn = fn, par.guess = c(0.05, 0.02), y = value, x = range, bsamples = bsamples, mask = mask[[i]], use.minpack.lm = FALSE)

        plot(fit.result_uncorrelated, main = paste0("Fit with a*exp(-E_0*t) with r = ", r_i[i]), ylab = "W(R,T)", xlab = "t", xlim = c(mask[[i]][1], mask[[i]][length(mask[[i]])]), ylim = rev(c(value[mask[[i]][1]], value[mask[[i]][length(mask[[i]])]])), log = "y")
        legend(x = "topright", legend = c(paste0("chi2 = ", round(fit.result_uncorrelated$chisqr, 2)), paste0("dof = ", fit.result_uncorrelated$dof), paste0("a = ", round(fit.result_uncorrelated$t0[1], 5)), paste0("da = ", round(fit.result_uncorrelated$se[1], 5)), paste0("E_0 = ", round(fit.result_uncorrelated$t0[2], 5)), paste0("dE_0 = ", round(fit.result_uncorrelated$se[2], 5))))

        plot(fit.result_uncorrelated, main = paste0("Fit with a*exp(-E_0*t) with r = ", r_i[i]), ylab = "W(R,T)", xlab = "t", log = "y")
        legend(x = "topright", legend = c(paste0("chi2 = ", round(fit.result_uncorrelated$chisqr, 2)), paste0("dof = ", fit.result_uncorrelated$dof), paste0("a = ", round(fit.result_uncorrelated$t0[1], 5)), paste0("da = ", round(fit.result_uncorrelated$se[1], 5)), paste0("E_0 = ", round(fit.result_uncorrelated$t0[2], 5)), paste0("dE_0 = ", round(fit.result_uncorrelated$se[2], 5))))

        dof <- fit.result_uncorrelated$dof
        chi2 <- fit.result_uncorrelated$chisqr

        print(paste0("R = ", r_i[i]))
        print("result")
        print(fit.result_uncorrelated$t0)
        print("error")
        print(fit.result_uncorrelated$se)
        print("chi2/dof")
        print(chi2 / dof)
        cat("\n")

        saveRDS(fit.result_uncorrelated, paste0(datas, folder, "/fit.result_uncorrelated_", i, ".rds"))
      }
      dev.off()
      #######################################
      
      message("plotting the potential...")
      ## performing uncorrelated fit to Wloop in a range specified by mask
      V <- rep(NA, length(r_i))
      bsV <- matrix(NA, nrow = boot.R, ncol = length(r_i))
      
      for (rr in seq_along(r_i))
      {
        tmp <- readRDS(paste0(datas, folder, "/fit.result_uncorrelated_", rr, ".rds"))
        
        V[rr] <- tmp$t0[[2]]
        bsV[, rr] <- tmp$t[, 2]
      }
      pdf(paste0(plots, folder, "/V_r.pdf"))
      plotwitherror(unlist(r_i), V, apply(bsV, 2, sd)
                    , ylab = "V(r)", xlab = "r"
                    , pch = 2, cex = 0.75
                    , main = paste0("V(r) for L = ", spatial_extent, ", beta = ", inv_coupling))
      grid()
      dev.off()
      rm(tmp)
    }
  }
}