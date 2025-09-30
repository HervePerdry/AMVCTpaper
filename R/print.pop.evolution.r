#'
#' @exportS3Method print pop.evolution
print.pop.evolution <- function(x) {
  p <- attr(x, "parameters")
  cat("Evolution of a population with parameters :\n")
  print(p)
  cat("\n")
  NextMethod()
}
