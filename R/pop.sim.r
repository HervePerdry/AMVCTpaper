#'  Simple population simulations
#'
#' @description Simulate the evolution of a population with Assortative Mating and Vertical Cultural Transmission, with unlinked causal variants
#' 
#' @param g0           Standard deviation of gametic value in a population without assortative mating
#' @param e            Standard deviation of environmental effects
#' @param r.ho         Correlation between mates 
#' @param nu           Correlation between parents' and offspring environmental effects
#' @param N            Number of causal variants 
#' @param nb.gen       Number of generations
#' @param pop.size     Population size
#' @param digest       Logical. Set to 'TRUE' to send back only a digest of the results
#' @param keep.N.kappa Logical. Set to 'TRUE' to keep trace of \eqn{N\overline\kappa(t)}{N kappa bar}.
#'
#' @details Will simulate a population under Assortative Mating and Vertical Cultural Transmission with the parameters
#' supplied. Generations are non overlapping.
#'
#' @details All SNPs have maf = 0.5 and are unlinked. Mate pairs are formed and eachh pair has two offsprings.
#'
#' @details Beware: setting 'keep.N.kappa' to 'TRUE' results
#' in lengthy computations.
#'
#' @return if 'digest' is 'TRUE', a data frame similar to the result of 'pop.evolution', with an additional
#' column for 'e' which contain the sd of the environmental components at each generation (this fluctuates
#' slightly around the value of the corresponding parameter, due to random sampling). If 'digest' is 'FALSE',
#' a list with member 'digest' (the digest as described before), 'G' (the matrix of genotypes at the last
#' generation), 'A' (the genetic value at the last generation), 'E' (the environmental value at the last 
#' generation), 'P' (the phenotypes at the last generation).
#'
#' @examples # a simulation with N = 100 loci only
#' R <- pop.sim(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4, pop.size = 25000, N = 100, nb.gen = 20, TRUE, TRUE)
#' # theoretical evolution
#' ev <- pop.evolution(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4, N = 100, nb.gen = 20) 
#' # theoretical equilibrium values
#' limits <- AMVCT(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4)
#'
#' par(mfrow=c(2,2))
#' # plotting evolution of rho
#' plot(R$t, R$rho, type = "l", xlab = "t", ylab = expression(rho))
#' lines(ev$t, ev$rho, col = "red")
#' abline(h = limits$rho, col = "red", lty = 3)
#'
#' # plotting evolution of N kappa bar
#' plot(R$t, R$N.kappa, type = "o", xlab = "t", ylab = expression(N * bar(kappa)))
#' lines(ev$t, ev$N.kappa, col = "red")
#' abline(h = 1/(1-limits$r.ga), col = "red", lty = 3)
#'
#' # plotting evolution of a
#' plot(R$t, R$a, type = "o", xlab = "t", ylab = "a") 
#' lines(ev$t, ev$a, col = "red")
#' abline(h = limits$a, col = "red", lty = 3)
#'
#' # plotting evolution of r.ga
#' plot(R$t, R$r.ga, type = "o", xlab = "t", ylab = expression(r[ga]))
#' lines(ev$t, ev$r.ga, col =" red")
#' abline(h = limits$r.ga, col = "red", lty = 3)
#'
#' @export
pop.sim <- function(g0, e, r.ho, nu, N = 100, nb.gen = 10, pop.size = 25000, digest = TRUE, keep.N.kappa = FALSE) {
  # size must be even
  pop.size <- (as.integer(pop.size) %/% 2L) * 2L

  # betas for equal variance (MAF = 0.5)
  beta <- rep( 1/sqrt(0.25*N), N ) * g0

  # initial population [one individual by column, all SNPs have MAF = 0.5]
  # paternal and maternal haplotypes
  Hp <- matrix( rbinom(pop.size*N, 1, 0.5), ncol = pop.size )
  Hm <- matrix( rbinom(pop.size*N, 1, 0.5), ncol = pop.size )
  # genotypes
  G <- Hp + Hm
  # genomic values
  A <- genomic.value(G, beta) # variance 2g0^2
  # environmental values
  E <- rnorm(pop.size, sd = e)
  # phenotypes
  P <- A + E

  # gametic value of parental gametes
  gvHp <- genomic.value(Hp, beta)
  gvHm <- genomic.value(Hm, beta)
  if(keep.N.kappa) {
    DEL <- cor( rbind(t(Hp), t(Hm)) )
    N.kappa <- sum(DEL)/N
  } else {
    N.kappa <- NA
  }

  R <- data.frame(e = sd(E), N.kappa = N.kappa, g = sd(A)/sqrt(2), r.ga = cor(gvHp, gvHm), rho = cor(A, E), a = sd(A), sigma2 = var(P))

  for(i in 1:nb.gen) {
    # make mate pairs
    pairs <- mate.pairs(P, r.ho)

    # gametes for 1st offspring
    H1 <- gametes(G)
    H1p <- H1[, pairs$I1, drop = FALSE] 
    H1m <- H1[, pairs$I2, drop = FALSE]

    # gametes for 2nd offspring
    H2 <- gametes(G)
    H2p <- H2[, pairs$I1, drop = FALSE] 
    H2m <- H2[, pairs$I2, drop = FALSE]

    # new matrix of genotypes
    G <- cbind(H1p + H1m, H2p + H2m )

    # Environment of offspring
    Ep <- E[pairs$I1]  # E of father
    Em <- E[pairs$I2]  # E of mother
    cpm <- cor(Ep, Em)
    cx <- nu/(1 + cpm)
    sdres <- sqrt(1 - 2*nu**2/(1 + cpm))*e
    E1 <- cx*(Ep + Em) + rnorm(pop.size/2, sd = sdres) # environment of offspring 1
    E2 <- cx*(Ep + Em) + rnorm(pop.size/2, sd = sdres) # environment of offspring 2

    # new vector of environments.
    Eo <- c(E1, E2)

    # concluding with new genomic value and new phenotype
    Ao <- genomic.value(G, beta)
    Po <- Ao + Eo
    
    # this is now the population
    A <- Ao; P <- Po; E <- Eo
    
    # the genomic values
    gvH1p <- genomic.value(H1p, beta)
    gvH1m <- genomic.value(H1m, beta)

    if(keep.N.kappa) {
      # DEL <- cor(t(H1)) # in fact could also be cor(rbind(t(H1),t(H2))), but this is lighter
      DEL <- cor(rbind(t(H1), t(H2)))
      N.kappa <- sum(DEL)/N
    } else {
      N.kappa <- NA
    }

    R <- rbind(R, data.frame(e = sd(E), N.kappa = N.kappa, g = sd(c(gvH1p, gvH1m)), r.ga = cor(gvH1p, gvH1m), rho = cor(A, E), a = sd(A), sigma2 = var(P)))
  } 

  R <- cbind(t = 0:nb.gen, R)
  attr(R, "parameters") <- c(g0 = g0, e = e, r.ho = r.ho, nu = nu, N = N, pop.size = pop.size)
  class(R) <- c("pop.evolution", "data.frame")

  if(digest) return(R)

  list(digest = R, G = G, A = A, E = E, P = P)
}


# calcule les genomic values (A = matrice de génomes *ou* d'haplotypes)
genomic.value <- function(A, beta) {
  N <- nrow(A)
  gv <- as.vector(beta %*% A)  
  gv - mean(gv)
}


mating.structure <- function(pop.size, r) {
  M <- MASS::mvrnorm(pop.size/2, c(0,0), matrix( c(1, r, r, 1), 2, 2 ))
  R <- matrix(rank(M), ncol = 2)
}

# renvoie deux listes d'indices qui correspondent aux deux mates d'une paire
mate.pairs <- function(P, r) {
  ms <- mating.structure(length(P), r)

  o <- order(P)
  I1 <- o[ms[,1]] 
  I2 <- o[ms[,2]]
  list(I1 = I1, I2 = I2)
}

# gametes. Tous les locis sont indépendants.
# G : une matrice de genotypes
gametes <- function(G) {
  N <- nrow(G)
  ps <- ncol(G)
  matrix(rbinom(N*ps, size = 1, p = as.vector(G/2)), nrow = N, ncol = ps)
}

