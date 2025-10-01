
pop <- function(g0, e, r.ho, nu, N) {
  a <- sqrt(2)*g0
  sigma2 <- a^2 + e^2
  list(g0 = g0, e = e, r.ho = r.ho, nu = nu, N = N, N.kappa = 1, g = g0, r.ga = 0, rho = 0, a = a, sigma2 = sigma2)
}


