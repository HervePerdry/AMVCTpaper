

next.gen <- function(x) {
  next.N.kappa <- with(x, 0.5*(1 + N.kappa + (1 - 1/N) * r.ga * N.kappa))
  next.g <- sqrt(next.N.kappa) * x$g0
  next.rga <- with(x, 0.5*r.ho/sigma2*(a + rho*e)^2 * g^2 / next.g^2 * (1 + r.ga))
  next.a <- sqrt(2*(1 + next.rga))*next.g
  next.rho <- with(x, nu * a/next.a * (rho + r.ho/sigma2*(rho*a + e)*(a + rho*e)) / (1 + r.ho/sigma2*(rho*a + e)^2))
  next.sigma2 <- with(x, next.a^2 + 2 * next.rho * next.a * e + e^2)
  list(g0 = x$g0, e = x$e, r.ho = x$r.ho, nu = x$nu, N = x$N, N.kappa = next.N.kappa, g = next.g, r.ga = next.rga, rho = next.rho, a = next.a, sigma2 = next.sigma2)
}
