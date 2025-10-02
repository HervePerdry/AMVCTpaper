source("multiple.pop.sim.r")

N <- 1000
for(pop in c(25,100)*1e3) {
  set.seed(1)
  si <- multiple.pop.sim(h2.0 = 0.2, r.ho = 0.6, nu = 0.6, pop.size = pop, N = N, nb.gen = 15, B = 50)
  name <- sprintf("simus-50-N%d-pop%dk.rds", N, pop/1000)
  saveRDS(si, file = name)
  cat(name, "- ok \n")
}

