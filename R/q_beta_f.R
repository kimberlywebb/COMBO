#' M-Step Expected Log-Likelihood with respect to Beta
#'
#' Objective function of the form:
#' \eqn{ Q_\beta = \sum_{i = 1}^N \Bigl[ \sum_{j = 0}^1 w_{ij} \text{log} \{ \pi_{ij} \}\Bigr]}.
#' Used to obtain estimates of \eqn{\beta} parameters.
#'
#' @param beta A numeric vector of regression parameters for the
#'   \code{Y} (true outcome) ~ \code{X} (predictor matrix of interest).
#' @param X A numeric design matrix.
#' @param w_mat Matrix of E-step weights obtained from \code{w_j}.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   can take.
#'
#' @return \code{q_beta_f} returns the negative value of the expected log-likelihood function,
#'   \eqn{ Q_\beta = \sum_{i = 1}^N \Bigl[ \sum_{j = 1}^2 w_{ij} \text{log} \{ \pi_{ij} \}\Bigr]},
#'   at the provided inputs.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_j.R
#'
#' @importFrom stats rnorm rgamma rmultinom optim
#'
q_beta_f <- function(beta, X, w_mat,
                     sample_size, n_cat){

  beta_mat = matrix(beta, ncol = 1)
  pi_terms = pi_compute(beta_mat, X, n = sample_size, n_cat)
  pi_terms = ifelse(pi_terms == 0, pi_terms + 0.00001, pi_terms)

  w_pi = c(w_mat) * log(c(pi_terms))
  result = -sum(sum_every_n(w_pi, n = sample_size))
  return(result)
}
