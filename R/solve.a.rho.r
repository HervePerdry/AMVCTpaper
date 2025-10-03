#' Solve a non-linear system of equation for 'a' and 'rho'
#'
#' @description given the natural parameters \eqn{g_0}{g0}, \eqn{e}, \eqn{r_{AM}}{r.ho} and \eqn{\nu}{nu},
#' computes \eqn{a} and \eqn{rho}.
#'
#' @param g0     Standard deviation of gametic value in a population without assortative mating
#' @param e      Standard deviation of environmental effects
#' @param r.ho   Correlation between mates 
#' @param nu     Correlation between parents' and offspring environmental effects
#' 
#' @details In the article text, 'r.ho' is denoted by \eqn{r_{ho}}{r_ho}.
#' The system formed by the equations 42 and 43 from the supplementary material is solved
#' to get the values of \eqn{a} and \eqn{\rho}{rho} at equilibrium. 
#'
#' @details The equations are solved numerically using the 'uniroot' function with default
#' parameters: the values obtained have a precision of roughly 1e-4.
#'
#' @details The system may have solutions even if the inequality ensuring that the environmental
#' variance matrix of a trio is positive definite (equation 26). This inequality is tested, and
#' if it is not verified, a warning is emitted and NAs are returned.
#'
#' @return a vector with named components 'a' and 'rho'.
#'
#' @examples solve.a.rho(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4)
#'
#' @export solve.a.rho
solve.a.rho <- function(g0, e, r.ho, nu) {
  rho <- rho(g0, e, r.ho, nu) 
  a <- a(rho, g0, e, r.ho)

  sigma2 <- a^2 + 2*rho*a*e + e^2
  if(nu^2 > 0.5*(1 + r.ho  * (rho*a + e)^2/sigma2)) {
    warning("'nu' is too high")
    a <- rho <- NA
  }
  c(a = a, rho = rho)
}
