#' M-Step Expected Log-Likelihood with respect to Delta
#'
#' Objective function of the form:
#' \eqn{Q_{\delta} = \sum_{i = 1}^N \Bigl[\sum_{j = 1}^2 \sum_{k = 1}^2 \sum_{\ell = 1}^2 w_{ij} y^*_{ik} \tilde{y}_{i \ell} \text{log} \{ \tilde{\pi}_{i \ell kj} \}\Bigr]}.
#' Used to obtain estimates of \eqn{\delta} parameters.
#'
#' @param delta_v A numeric array of regression parameters for the second-stage observed
#'   outcome mechanism, \eqn{\tilde{Y} | Y^*, Y}
#'   (second-stage observed outcome, given the first-stage observed outcome and the true outcome) ~ \code{V} (misclassification
#'   predictor matrix). The \eqn{\delta} vector is obtained from the array form. In array form,
#'   the first dimension (matrix rows) of \code{delta}
#'   corresponds to parameters for the \eqn{\tilde{Y} = 1}
#'   second-stage observed outcome, with the dimensions of the \code{V}
#'   The second dimension (matrix columns) correspond to the first-stage
#'   observed outcome categories \eqn{Y^* \in \{1, 2\}}. The third dimension of
#'   \code{delta_start} corresponds to to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}. The numeric vector \eqn{\delta} is obtained by
#'   concatenating the delta array, i.e. \code{delta_v <- c(delta_array)}.
#' @param V A numeric design matrix.
#' @param obs_Ystar_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param obs_Ytilde_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \eqn{\tilde{Y}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param w_mat Matrix of E-step weights obtained from \code{w_j_2stage}.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{V}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcomes can take.
#'
#' @return \code{q_beta_f} returns the negative value of the expected log-likelihood function,
#'   \eqn{Q_{\delta} = \sum_{i = 1}^N \Bigl[\sum_{j = 1}^2 \sum_{k = 1}^2 \sum_{\ell = 1}^2 w_{ij} y^*_{ik} \tilde{y}_{i \ell} \text{log} \{ \tilde{\pi}_{i \ell kj} \}\Bigr]},
#'   at the provided inputs.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include pitilde_compute.R
#' @include w_j_2stage.R
#'
#' @importFrom stats rnorm rgamma rmultinom optim
#'
#' @examples
q_delta_f <- function(delta_v, V, obs_Ystar_matrix, obs_Ytilde_matrix, w_mat,
                      sample_size, n_cat){

  delta_mat = array(delta_v, dim = c(ncol(V), 2, 2))

  pitilde_terms_array = pitilde_compute(delta_mat, V, sample_size, n_cat)

  big_Ystar_matrix = array(c(rep(obs_Ystar_matrix[,1], 2),
                             rep(obs_Ystar_matrix[,2], 2),
                             rep(obs_Ystar_matrix[,1], 2),
                             rep(obs_Ystar_matrix[,2], 2)),
                           dim = c(2 * sample_size, 2, 2))

  big_Ytilde_matrix = array(rep(c(obs_Ytilde_matrix[,1],
                                  obs_Ytilde_matrix[,2]), 4),
                            dim = c(2 * sample_size, 2, 2))

  big_w_matrix = array(c(rep(w_mat[,1], 4),
                         rep(w_mat[,2], 4)),
                       dim = c(2 * sample_size, 2, 2))

  summand = big_Ystar_matrix * big_Ytilde_matrix * big_w_matrix * log(pitilde_terms_array)
  result = -sum(summand)
  return(result)
}
