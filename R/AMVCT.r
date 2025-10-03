#' Assortative Mating with Vertical Cultural Transmission
#'
#' @param h2.0   Heritability in a population without assortative mating
#' @param g0     Standard deviation of gametic value in a population without assortative mating
#' @param e      Standard deviation of environmental effects
#' @param r.ho   Correlation between mates 
#' @param nu     Correlation between parents' and offspring environmental effects
#'
#' @details The user must provide 'r.ho' and 'nu', and either 'g0' and 'e', or 'h2.0'. 
#' The function will solve the non-linear system of equations with 'solve.a.rho', and
#' compute the principal quantities of interest defined in the paper. A method is
#' defined to have them printed nicely.
#'
#' @return an object of class 'AMVCT', which is a list with components
#' 'g0', 'e', 'r.ho', 'nu', 'rho', 'a', 'sigma2', 'r.ga', 'g', 'h2.SNP',
#' 'decompose.sigma2', 'cor.mates' (corresponding to table 1 of the supplementary material)
#' 'cor.parent.offspring' (corresponding to table 2 of the supplementary material)
#' 'cor.pheno.parent.offspring' (correlation of the phenotypes of parent and offspring).
#'
#' @examples AMVCT(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4)
#'
#' @export
AMVCT <- function(h2.0, g0 = sqrt(h2.0/2), e = sqrt(1 - h2.0), r.ho, nu) {
  if(!missing(h2.0)) {
    if(g0 != sqrt(h2.0/2) | e != sqrt(1-h2.0))
      stop("specify either h2.0 or g0, e")
  }  
  x <- solve.a.rho(g0, e, r.ho, nu)
  rho <- x["rho"]
  a <- x["a"]

  sigma2 <- a^2 + 2*rho*a*e + e^2

  decompose.sigma2 <-  c(a2 = a**2, deux.rho.a.e = 2*rho*a*e, e2 = e**2)

  # correlation of mates  
  C1 <- matrix( c(1, rho, rho, 1), 2, 2)
  C2 <- matrix( c(r.ho*(a + rho*e)**2 / sigma2,        r.ho*(a + rho*e)*(rho*a + e)/sigma2, 
                  r.ho*(a + rho*e)*(rho*a + e)/sigma2, r.ho*(rho*a + e)**2 / sigma2), 
               2, 2)
  cor.mates <- rbind( cbind(C1,C2), cbind(C2, C1))
  rownames(cor.mates) <- colnames(cor.mates) <- c("A1", "E1", "A2", "E2")

  # covariance of mates
  cov.mates <-  cor.mates * tcrossprod( c(a,e,a,e) )

  # gametic correlation
  r.ga <- cor.mates[1,3] / (2 - cor.mates[1,3])

  # gametic sd
  g <- g0/sqrt(1 - r.ga)
 
  # SNP heritability
  h2.SNP <- (a + rho * e)**2/sigma2
  
  # correlations parent offspring
  cor.parent.offspring <- matrix( c( .5*(1 + cor.mates[1,3]), rho, 0.5*(cor.mates[2,1] + cor.mates[2,3]), nu) , 2,2 , byrow = TRUE)
  rownames(cor.parent.offspring) <- c("A1", "E1")
  colnames(cor.parent.offspring) <- c("A3", "E3")

  # phenotype correlation parent offspring
  cov.pheno.parent.offspring <- cor.parent.offspring["A1", "A3"] * a^2 + cor.parent.offspring["A1", "E3"] * a*e + 
                                cor.parent.offspring["E1", "A3"] * a*e + cor.parent.offspring["E1", "E3"] * e^2
  cor.pheno.parent.offspring <- cov.pheno.parent.offspring / sigma2

  L <- list(g0 = g0, e = e, r.ho = r.ho, nu = nu, rho = rho, a = a, sigma2 = sigma2, r.ga = r.ga, g = g, h2.SNP = h2.SNP,
            decompose.sigma2 = decompose.sigma2, cor.mates = cor.mates, cor.parent.offspring = cor.parent.offspring,
            cor.pheno.parent.offspring = cor.pheno.parent.offspring)
  class(L) <- "AMVCT"
  L
}
