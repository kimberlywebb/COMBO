#' Compute Probability of Each True Outcome, for Every Subject
#'
#' @param beta A numeric column matrix of regression parameters for the
#'   \code{Y} (true outcome) ~ \code{X} (predictor matrix of interest).
#' @param X A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   can take.
#'
#' @return \code{pi_compute} returns a matrix of probabilities,
#'   \eqn{P(Y_i = j | X_i) = \frac{\exp(X_i \beta)}{1 + \exp(X_i \beta)}}
#'   for each of the \eqn{i = 1, \dots,} \code{n} subjects. Rows of the matrix
#'   correspond to each subject. Columns of the matrix correspond to the true outcome
#'   categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include sum_every_n1.R
#'
#' @importFrom stats rnorm
#'
pi_compute <- function(beta, X, n, n_cat){
  exp_xb = exp(X %*% beta)
  pi_result = exp_xb[,1] / rep(sum_every_n1(exp_xb[,1], n), n_cat - 1)

  pi_matrix = matrix(c(pi_result, 1 - pi_result),
                     ncol = n_cat, nrow = n,
                     byrow = FALSE)
  return(pi_matrix)
}
