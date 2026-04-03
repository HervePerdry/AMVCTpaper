require(AMVCTpaper)
source("fig.evol.r")
source("fig.nu.r")

#-----------------------------------------------------
fig2 <- function() {
  par(cex = 2.2, mar = c(5,5,1,1)) 
  fig.nu.rho(h2.0 = 0.2) 
}

tiff("fig-2-nu-rho.tiff", width = 1500, height = 1000, compression = "jpeg")
fig2()
dev.off()

pdf("fig-2-nu-rho.pdf", width = 1500/72, height = 1000/72)
fig2()
dev.off()

#-----------------------------------------------------
fig3 <- function() {
  par(cex = 2.2, mar = c(5,5,1,1)) 
  fig.nu.h2snp(h2.0 = 0.2)
}

tiff("fig-3-nu-h2snp.tiff", width = 1500, height = 1000, compression = "jpeg")
fig3()
dev.off()

pdf("fig-3-nu-h2snp.pdf", width = 1500/72, height = 1000/72)
fig3()
dev.off()

#-----------------------------------------------------
fig4 <- function() {
  par(cex = 2.2, mar = c(5,5,3,1)) 
  fig.nu.variance.comp(h2.0 = 0.2, r = 0.6)

  legend(0.01, 1.2, pch = c(15,15,22), col = c("gray50", "gray70", "black"), horiz = TRUE, 
        legend = c(expression(frac(a^2,sigma^2)), expression(paste(frac(2*rho*a*e,sigma^2), "   ")), expression(frac(e^2, sigma^2))) , 
        bty = "n", pt.cex = 2, bg = "black")
  legend(0.26, 1.16, lty = 2, lwd = 2, legend = expression({h^2}[SNP]), bty = "n")
}

tiff("fig-4-nu-varcomp.tiff", width = 1500, height = 1000, compression = "jpeg")
fig4()
dev.off()

pdf("fig-4-nu-varcomp.pdf", width = 1500/72, height = 1000/72)
fig4()
dev.off()

#-----------------------------------------------------
fig5 <- function() {
  par(cex = 2.2, mar = c(5,5,1,1)) 
  fig.nu.gamma(h2.0 = 0.2)
}

tiff("fig-5-nu-gamma.tiff", width = 1500, height = 1000, compression = "jpeg")
fig5()
dev.off()

pdf("fig-5-nu-gamma.pdf", width = 1500/72, height = 1000/72)
fig5()
dev.off()


# ----------------------------------------------------
# reading simulation results for two next figures
A.25  <- readRDS("simus-50-N1000-pop25k.rds")
A.100 <- readRDS("simus-50-N1000-pop100k.rds")

# equilibrium values
eq <- AMVCT(h2.0 = 0.2, r = 0.6, nu = 0.6)

#-----------------------------------------------------
fig6 <- function() {
  par(mfrow = c(1,2), cex = 2.2, mar = c(5, 5, 4, 1))

  fig.evol(A.25, "rho", main = "population size 25,000", ylim = c(0, 0.25))
  axis(1, at = 0:15, rep("", 16))
  abline(h = eq$rho, lty = 2)

  fig.evol(A.100, "rho", main = "population size 100,000", ylim = c(0, 0.25))
  axis(1, at = 0:15, rep("", 16))
  abline(h = eq$rho, lty = 2)
}

tiff("fig-6-evol-rho.tiff", width = 2000, height = 1000, compression = "jpeg")
fig6()
dev.off()

pdf("fig-6-evol-rho.pdf", width = 2000/72, height = 1000/72)
fig6()
dev.off()

#-----------------------------------------------------
fig7 <- function() {
  par(mfrow = c(1,3), cex = 2, mar = c(5, 5, 2, 1))

  fig.evol(A.25, "r.ga")
  axis(1, at = 0:15, rep("", 16))
  abline(h = eq$r.ga, lty = 2)

  fig.evol(A.25, "a2")
  axis(1, at = 0:15, rep("", 16))
  abline(h = eq$a^2, lty = 2)

  fig.evol(A.25, "sigma2")
  axis(1, at = 0:15, rep("", 16))
  abline(h = eq$sigma2, lty = 2)
}

tiff("fig-7-evol-rga-a2-sigma2.tiff", width = 2100, height = 700, compression = "jpeg")
fig7()
dev.off()

pdf("fig-7-evol-rga-a2-sigma2.pdf", width = 2100/72, height = 700/72)
fig7()
dev.off()

