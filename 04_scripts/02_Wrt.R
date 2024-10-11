# plotting W_r(t) for every r and t


for(beta in betaArray){

  print(beta)

  source("~/projects/stepscaling/RU1/03_functions/header.R")
  source("~/projects/stepscaling/RU1/03_functions/myFunctions.R")

  tmp <- inputFileName(spatialExtent, temporalExtent, beta, sizeWLoops)
  writePath <- writePath(tmp)

  if (!dir.exists(writePath)) {
    dir.create(writePath, recursive = TRUE)
    cat("Directory created:", writePath, "\n")
  }

  wr<- plotWrt(spatialExtent
          , temporalExtent
          , beta
          , sizeWLoops
          , thermSkip
          , performFit = FALSE)


  saveRDS(wr, paste0(writePath, "wrt.rds"))

}
# print(wr)
# print(wr[[3]]$cf.tsboot$t0[[1]])
# phi1 <- acos(wr[[3]]$cf.tsboot$t0[[1]])
# phi2 <- acos(wr[[1]]$cf.tsboot$t0[[15]])
# phi3 <- acos(wr[[2]]$cf.tsboot$t0[[1]])
#
# for(i in c(1,2)){
#   for(j in c(1,2)){
#     for(k in c(1,2)){
#       a <- (-1) ** i * phi1 + (-1) ** j * phi2 + (-1) ** k * phi3
#       print(a/pi)
#   }}}

x1 <- c(1:16)
x2 <- c(31:16)

y1 <- wr[[1]]$cf.tsboot$t0
y2 <- wr[[2]]$cf.tsboot$t0

dy1 <- apply(wr[[1]]$cf.tsboot$t, 2, sd)
dy2 <- apply(wr[[2]]$cf.tsboot$t, 2, sd)

labello <- c(1, 4, 7, 10, 13, "T = 16", 13, 10, 7, 4, 1)

pdf("L3T16WilsonLoopR1R2JointZoom.pdf")
plotwitherror(x1, y = y1, dy = dy1, xlim = c(12,24), ylim = c(0.01,0.04), xaxt = 'n',col = "red", ylab = "<W(r,t)>", xlab = "t", main = "<W> for r=1, r=2, L=3, T=16, beta=1.65, D=2+1 pure gauge QED")
axis(1, at = seq(1, 31, 3), labels = labello)
plotwitherror(x2, y = y2, dy = dy2,xlim = c(12,24), ylim = c(0.01, 0.04), rep = TRUE, col = "blue")

legend(x = "topright",          # Position
      legend = c("W(r=1)", "W(r=2)"),  # Legend texts
      col = c("red", "blue"),           # Line colors
      lwd = 2)                 # Line width

dev.off()


