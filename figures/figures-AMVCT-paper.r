require(AMVCTpaper)
source("fig.evol.r")
source("fig.nu.r")

A.25  <- readRDS("simus-50-N1000-pop25k.rds")
A.100 <- readRDS("simus-50-N1000-pop100k.rds")
eq <- AMVCT(h2.0 = 0.2, r = 0.6, nu = 0.6)

#-----------------------------------------------------
tiff("fig-2-evol-rho.tiff", width = 2000, height = 1000, compression = "jpeg")
par(mfrow = c(1,2), cex = 2.2, mar = c(5, 5, 4, 1))
fig.evol(A.25, "rho", main = "population size 25,000", ylim = c(0, 0.25))
axis(1, at = 0:15, rep("", 16))
abline(h = eq$rho, lty = 2)
fig.evol(A.100, "rho", main = "population size 100,000", ylim = c(0, 0.25))
axis(1, at = 0:15, rep("", 16))
abline(h = eq$rho, lty = 2)
dev.off()

#-----------------------------------------------------
tiff("fig-3-evol-rga-a2-sigma2.tiff", width = 2100, height = 700, compression = "jpeg")
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

dev.off()

#-----------------------------------------------------
tiff("fig-4-nu-rho.tiff", width = 1500, height = 1000, compression = "jpeg")
par(cex = 2.2, mar = c(5,5,1,1)) 
fig.nu.rho(h2.0 = 0.2) 
dev.off()

#-----------------------------------------------------
tiff("fig-5-nu-h2snp.tiff", width = 1500, height = 1000, compression = "jpeg")
par(cex = 2.2, mar = c(5,5,1,1)) 
fig.nu.h2snp(h2.0 = 0.2)
dev.off()

#-----------------------------------------------------
tiff("fig-6-nu-gamma.tiff", width = 1500, height = 1000, compression = "jpeg")
par(cex = 2.2, mar = c(5,5,1,1)) 
fig.nu.gamma(h2.0 = 0.2)
dev.off()

#-----------------------------------------------------
tiff("fig-7-nu-varcomp.tiff", width = 1500, height = 1000, compression = "jpeg")
par(cex = 2.2, mar = c(5,5,3,1)) 
fig.nu.variance.comp(h2.0 = 0.2, r = 0.6)

legend(0.01, 1.2, pch = c(15,15,22), col = c("gray50", "gray70", "black"), horiz = TRUE, 
      legend = c(expression(frac(a^2,sigma^2)), expression(paste(frac(2*rho*a*e,sigma^2), "   ")), expression(frac(e^2, sigma^2))) , 
      bty = "n", pt.cex = 2, bg = "black")
legend(0.26, 1.16, lty = 2, lwd = 2, legend = expression({h^2}[SNP]), bty = "n")

dev.off()
