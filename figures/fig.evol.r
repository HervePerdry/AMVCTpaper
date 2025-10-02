require(AMVCTpaper)

# this function will take the result of multiple.pop.sim and make a figure
# (cf examples in figures-AMVCT-paper.r)
fig.evol <- function(x, par = c("rho", "r.ga", "g", "a", "a2", "sigma2"), ylim, xlab = "t", ylab, ...) {
  par <- match.arg(par)
  if(par == "a2")
    return(fig.evol.a2(x, ylim, xlab, ylab, ...))

  ev <- x$evolution
  R <- x$simulations
  B <- length(R)
  gens <- seq(0, nrow(R[[1]]) - 1)
  if(missing(ylim)) {
    mi <- min(sapply(R, function(x) min(x[[par]])))
    ma <- max(sapply(R, function(x) max(x[[par]])))
    ylim <- c(mi, ma)
  }
  if(missing(ylab)) {
    if(par == "rho") ylab <- expression(rho(t))
    else if(par == "r.ga") ylab <- expression(r[ga](t))
    else if(par == "sigma2") ylab <- expression(sigma^2(t))
    else ylab <- par
  }
  plot(gens, ev[[par]], ylim = ylim, type = "n", xlab = xlab, ylab = ylab, ...)
  for(i in 1:B) {
    lines(gens, R[[i]][[par]], type = "l", col = "gray70")
  }
  lines(gens, ev[[par]], lwd = 2)
}

# there's no column for a2 so we need a special function
fig.evol.a2 <- function(x, ylim, xlab = "t", ylab, ...)  {
  ev <- x$evolution
  R <- x$simulations
  B <- length(R)
  gens <- seq(0, nrow(R[[1]]) - 1)
  if(missing(ylim)) {
    mi <- min(sapply(R, function(x) min(x$a^2)))
    ma <- max(sapply(R, function(x) max(x$a^2)))
    ylim <- c(mi, ma)
  }
  if(missing(ylab)) {
    ylab <- expression({a^2}(t))
  }
  plot(gens, ev$a^2, ylim = ylim, type = "n", xlab = xlab, ylab = ylab, ...)
  for(i in 1:B) {
    lines(gens, R[[i]]$a^2, type = "l", col = "gray70")
  }
  lines(gens, ev$a^2, lwd = 2)
}


