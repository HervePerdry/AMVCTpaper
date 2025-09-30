# Solving equation 41 of supplementary material

a <- function(rho, g0, e, r.AM) {
  if(rho == 1) 
    return( sqrt(2*g0^2 / (1 - r.AM)) )

  f <- function(a) { 
    r.A1.A2 <- r.AM * (a + rho*e)^2 / (a^2 + 2 * rho*a*e + e^2)
    a^2 * (1 - r.A1.A2) - 2 * g0^2 
  }
  ur <- uniroot(f, c(0, sqrt(2*g0^2/(1 - r.AM)) + 1e-10))
  ur$root
}

