#' Expit function
#'
#' \eqn{\frac{\exp\{x\}}{1 + \exp\{x\}}}
#'
#' @param x A numeric value or vector to compute the expit function on.
#'
#' @return \code{expit} returns the result of the function
#'   \eqn{f(x) = \frac{\exp\{x\}}{1 + \exp\{x\}}} for a given \code{x}.
#'
#' @examples \dontrun{
#' x <- c(1, 2, 3)
#'
#' expit(x)
#'
#' exp(x) / (1 + exp(x))
#' }
expit <- function(x){
  expit_return = exp(x) / (1 + exp(x))
  return(expit_return)
  }
