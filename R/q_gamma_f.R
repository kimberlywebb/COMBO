#' M-Step Expected Log-Likelihood with respect to Gamma
#'
#' Objective function of the form:
#' \eqn{Q_{\gamma} = \sum_{i = 1}^N \Bigl[\sum_{j = 1}^2 \sum_{k = 1}^2 w_{ij} y^*_{ik} \text{log} \{ \pi^*_{ikj} \}\Bigr]}.
#' Used to obtain estimates of \eqn{\gamma} parameters.
#'
#' @param gamma_v A numeric vector of regression parameters for the observed
#'   outcome mechanism, \code{Y* | Y}
#'   (observed outcome, given the true outcome) ~ \code{Z} (misclassification
#'   predictor matrix). In matrix form, the gamma parameter matrix rows
#'   correspond to parameters for the \code{Y* = 0}
#'   observed outcome, with the dimensions of \code{Z}.
#'   In matrix form, the gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{j = 1, \dots,} \code{n_cat}. The numeric vector \code{gamma_v} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_v <- c(gamma_matrix)}.
#' @param Z A numeric design matrix.
#' @param obs_Y_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param w_mat Matrix of E-step weights obtained from \code{w_j}.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{q_beta_f} returns the negative value of the expected log-likelihood function,
#'   \eqn{Q_{\gamma} = \sum_{i = 1}^N \Bigl[\sum_{j = 1}^2 \sum_{k = 1}^2 w_{ij} y^*_{ik} \text{log} \{ \pi^*_{ikj} \}\Bigr]},
#'   at the provided inputs.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_j.R
#'
#' @importFrom stats rnorm rgamma rmultinom optim
#'
q_gamma_f <- function(gamma_v, Z, obs_Y_matrix, w_mat,
                      sample_size, n_cat){

  gamma_mat = matrix(gamma_v, ncol = n_cat, byrow = FALSE)

  pistar_terms_v = pistar_compute(gamma_mat, Z, sample_size, n_cat)

  big_Ystar_matrix = matrix(rep(c(obs_Y_matrix), n_cat),
                            ncol = n_cat, byrow = FALSE)
  big_w_matrix = do.call(rbind, replicate(n_cat, w_mat, simplify = FALSE))

  summand = big_Ystar_matrix * big_w_matrix * log(pistar_terms_v)
  result = -sum(summand)
  return(result)
}
