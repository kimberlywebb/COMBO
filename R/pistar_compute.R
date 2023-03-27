#' Compute Conditional Probability of Each Observed Outcome Given Each True Outcome, for Every Subject
#'
#' @param gamma A numeric matrix of regression parameters for the observed
#'   outcome mechanism, \code{Y* | Y}
#'   (observed outcome, given the true outcome) ~ \code{Z} (misclassification
#'   predictor matrix). Rows of the matrix correspond to parameters for the \code{Y* = 1}
#'   observed outcome, with the dimensions of \code{Z}.
#'   Columns of the matrix correspond to the true outcome categories
#'   \eqn{j = 1, \dots,} \code{n_cat}.
#' @param Z A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{pistar_compute} returns a matrix of conditional probabilities,
#'   \eqn{P(Y_i^* = k | Y_i = j, Z_i) = \frac{\text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}{1 + \text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}}
#'   for each of the \eqn{i = 1, \dots,} \code{n} subjects. Rows of the matrix
#'   correspond to each subject and observed outcome. Specifically, the probability
#'   for subject \eqn{i} and observed category $1$ occurs at row \eqn{i}. The probability
#'   for subject \eqn{i} and observed category $2$ occurs at row \eqn{i +} \code{n}.
#'   Columns of the matrix correspond to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include sum_every_n.R
#' @include sum_every_n1.R
#'
#' @importFrom stats rnorm
#'
pistar_compute <- function(gamma, Z, n, n_cat){

  exp_zg = exp(Z %*% gamma)
  pi_denominator = apply(exp_zg, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result = exp_zg / rbind(pi_denominator)

  pistar_matrix = rbind(pi_result,
                        1 - apply(pi_result,
                                  FUN = sum_every_n, n = n,
                                  MARGIN = 2))

  return(pistar_matrix)
}
