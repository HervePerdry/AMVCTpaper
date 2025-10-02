require(AMVCTpaper)
require(parallel)

# parameters h2.0, g0, e, r.ho, nu, pop.size, N, nb.gen, keep.N.kappa = as in pop.sim and pop.evolution
# B = number of simulations to run
# n.core = number of cores to use (defaults to half of available cores at most)
#
# returns a list with components 'simulations' (itself a list of the results of pop.sim)
# and 'evolution' (the result of pop.evolution)

multiple.pop.sim <- function(h2.0 = 0.2, g0 = sqrt(h2.0/2), e = sqrt(1-h2.0), r.ho = 0.6, nu = 0.6, 
                             pop.size = 25000, N = 1000, nb.gen = 15, keep.N.kappa = FALSE, B = 20, n.cores) {
  if(missing(n.cores)) {
    n.cores <- parallel::detectCores() %/% 2
    n.cores <- min(B, n.cores)
  }

  # force evaluation of all arguments before defining f  
  force(g0)
  force(e)
  force(r.ho)
  force(nu)
  force(pop.size)
  force(N)
  force(nb.gen)
  force(keep.N.kappa)

  f <- function(i) {
    cat("start", i, "\n")
    S <- AMVCTpaper::pop.sim(g0 = g0, e = e, r.ho = r.ho, nu = nu, pop.size = pop.size, N = N, nb.gen = nb.gen, keep.N.kappa = keep.N.kappa)
    cat("end", i, "\n")
    S
  }

  cl <- makePSOCKcluster(n.cores, outfile = "")
  clusterExport(cl, "f", environment())
  clusterSetRNGStream(cl) 
  R <- parLapply(cl, 1:B, f)
  stopCluster(cl)

  ev <- pop.evolution(g0 = g0, e = e, r.ho = r.ho, nu = nu, N = N, nb.gen = nb.gen)
  list(simulations = R, evolution = ev)
}

