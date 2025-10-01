#' Population evolution
#'
#' @description Compute all parameters of the evolution of a population which is initially 
#' at gametic equilibrium, without assortative mating.
#' 
#' @param g0     Standard deviation of gametic value in a population without assortative mating
#' @param e      Standard deviation of environmental effects
#' @param r.ho   Correlation between mates 
#' @param nu     Correlation between parents' and offspring environmental effects
#' @param N      Number of causal variants 
#' @param nb.gen Number of generations
#' 
#' @details The details of the computation are given in the supplementary material of the paper,
#' in particular in appendix B.
#' 
#' @return A data frame with class "pop.evolution", with columns
#' 't' (the generation, from 0 to nb.gen), 'N.kappa' (the value of \eqn{N \overline\kappa(t)}{N kappa bar(t)}),
#' 'g' (the standard deviation of the gametic value, \eqn{g(t)}), 'r.ga' (the gametic correlation, \eqn{r_{ga}(t)}{r_ga(t)}), 
#' 'rho' (the gene-environment correlation, \eqn{\rho(t)}{rho(t)}), 'a' (the standard deviation of the genetic value, \eqn{g(t)}),
#' and 'sigma2' (the variance of the phenotype, \eqn{sigma^2(t)}{sigma^2(t)}).
#' 
#' @return The parameters 'g0', 'e', 'r.ho', 'nu' and 'N' are stored in an attribute named 'parameters'.
#'
#' @examples 
#' pop.evolution(1, 1, 0.6, 0.2, Inf, 10)
#'
#' @export
pop.evolution <- function(g0, e, r.ho, nu, N, nb.gen) {
  L <- pop(g0, e, r.ho, nu, N)
  R <- as.data.frame(L)
  for(i in 1:nb.gen) {
    L <- next.gen(L)
    R <- rbind(R, L)
  }
  R <- cbind(t = 0:nb.gen, R[ , -(1:5) ])
  attr(R, "parameters") <- c(g0 = g0, e = e, r.ho = r.ho, nu = nu, N = N)
  class(R) <- c("pop.evolution", "data.frame")
  R
}


