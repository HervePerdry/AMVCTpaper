# gamma is the inflation factor beta^/ beta
# R = values of r.ho for which the plot is done
fig.nu.gamma <- function(h2.0 = 0.2, g0 = sqrt(h2.0/2), e = sqrt(1-h2.0), R = seq(0, 0.8, by = 0.2), ... ) {

  L0 <- list(g0 = g0, e = e, r.ho = NA, nu = NA)
  N <- 101
  NU <- seq(0, 1, length = N)

  R <- sort(R, decreasing = TRUE)
  for(j in seq_along(R)) {
    L0$r.ho <- R[j];
    GAMMA <- numeric(N)
    for(i in 1:N) {
      L0$nu <- NU[i]
      L <- do.call(AMVCT, L0)
      GAMMA[i] <- with(L, (1 + rho * e / a)*(1 + r.ga)/(1 - r.ga))
    }
    par(xpd = TRUE) # to write outside the frame
    if(j == 1)
      plot(NU, GAMMA, type = "l", xlab = expression(nu), ylab = expression(gamma), ylim = c(1, max(GAMMA, na.rm = TRUE)), bty = "l", ...)
    else
      lines(NU, GAMMA, type = "l")

    k <- which(is.na(GAMMA))[1]-1
    points(NU[k], GAMMA[k],pch=16)
    text(NU[k], GAMMA[k], bquote(r[ho] == .(L$r.ho)), pos = 4, cex = 0.9)
  }
}

# R = values of r.ho for which the plot is done
fig.nu.h2snp <- function(h2.0 = 0.2, g0 = sqrt(h2.0/2), e = sqrt(1-h2.0), R = seq(0, 0.8, by = 0.2), ylim = c(0,1), ... ) {

  L0 <- list(g0 = g0, e = e, r.ho = NA, nu = NA)
  N <- 101
  NU <- seq(0, 1, length = N)

  for(j in seq_along(R)) {
    L0$r.ho <- R[j];
    H2SNP <- numeric(N)
    A <- numeric(N)
    for(i in 1:N) {
      L0$nu <- NU[i]
      L <- do.call(AMVCT, L0)
      H2SNP[i] <- with(L, (a + rho*e)^2/sigma2)
    }
    par(xpd = TRUE) # to write outside the frame
    if(j == 1)
      plot(NU, H2SNP, type = "l", xlab = expression(nu), ylab = expression(h[SNP]^2), ylim = ylim, bty = "l", ...)
    else
      lines(NU, H2SNP, type = "l")

    k <- which(is.na(H2SNP))[1]-1
    points(NU[k], H2SNP[k],pch=16)
    text(NU[k], H2SNP[k], bquote(r[ho] == .(L$r.ho)), pos = 4, cex = 0.9)
  }
}


# R = values of r.ho for which the plot is done
fig.nu.rho <- function(h2.0 = 0.2, g0 = sqrt(h2.0/2), e = sqrt(1-h2.0), R = seq(0, 0.8, by = 0.2), ylim = c(0,1), ... ) {

  L0 <- list(g0 = g0, e = e, r.ho = NA, nu = NA)
  N <- 101
  NU <- seq(0, 1, length = N)

  for(j in seq_along(R)) {
    L0$r.ho <- R[j];
    RHO <- numeric(N)
    for(i in 1:N) {
      L0$nu <- NU[i]
      L <- do.call(AMVCT, L0)
      RHO[i] <- L$rho
    }
    par(xpd = TRUE) # to write outside the frame
    if(j == 1)
      plot(NU, RHO, type = "l", xlab = expression(nu), ylab = expression(rho), ylim = ylim, bty = "l", ...)
    else
      lines(NU, RHO, type = "l")

    k <- which(is.na(RHO))[1]-1
    points(NU[k], RHO[k],pch=16)
    text(NU[k], RHO[k], bquote(r[ho] == .(L$r.ho)), pos = 4, cex = 0.9)
  }
}



fig.nu.variance.comp <- function(h2.0 = 0.2, g0 = sqrt(h2.0/2), e = sqrt(1-h2.0), r.ho = 0.6, ... ) {

  L0 <- list(g0 = g0, e = e, r.ho = r.ho, nu = NA)
  N <- 101
  NU <- seq(0, 1, length = N)

  H2SNP <- numeric(N)
  VAR <- matrix(0, nrow = 3, ncol = N)
  for(i in 1:N) {
    L0$nu <- NU[i]
    L <- do.call(AMVCT, L0)
    H2SNP[i] <- with(L, (a + rho*e)^2/sigma2)
    VAR[,i] <- L$decompo/L$sigma2
  }
  n <- which(is.na(H2SNP))[1] - 1
  NU <- NU[1:n]
  H2SNP <- H2SNP[1:n]
  VAR <- VAR[,1:n]

  par(xpd = TRUE) # to write outside the frame
  plot(NU, VAR[1,], type = "n", xlab = expression(nu), ylab = "", ylim = c(0,1), bty = "n", xaxt = "n", ...)
 
  polygon( c(NU, rev(NU)), c(VAR[1,], rep(0,n)), col = "gray50", border = NA)
  polygon( c(NU, rev(NU)), c(VAR[1,], rev(VAR[1,]+VAR[2,])), col = "gray70", border = NA)
  polygon( c(NU, rev(NU)), rep(0:1, each = n))

  lines(NU, H2SNP, lty = 2, lwd = 2)
  ti <- seq(0, max(NU), by = 0.2)
  axis(1, c(ti, max(NU)), c(ti,""))
}
