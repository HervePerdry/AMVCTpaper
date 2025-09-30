rho <- function(g0, e, r.AM, nu) {
  f <- function(rho) {
    a <- a(rho, g0, e, r.AM)
    sigma2 <- a^2 + 2*rho*a*e + e^2
    if(rho == 1) 
      (nu - 1) * (sigma2 + r.AM * (a + e)**2) # get rid of rounding error
    else
     nu * (rho * sigma2 + r.AM * (rho * a + e)*(rho * e + a)) - rho * (sigma2 + r.AM * (rho*a + e)**2)
  }  
  ur <- uniroot(f, c(0, 1))
  ur$root
}
