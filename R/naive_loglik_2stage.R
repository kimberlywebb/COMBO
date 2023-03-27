#' Observed Data Log-Likelihood Function for Estimation of the Naive Two-Stage Misclassification Model
#'
#' @param param_current A numeric vector of regression parameters, in the order
#'   \eqn{\beta, \delta}. The \eqn{\delta} vector is obtained from the matrix form.
#'   In matrix form, the gamma parameter matrix rows
#'   correspond to parameters for the \eqn{\tilde{Y} = 1}
#'   observed outcome, with the dimensions of \code{V}.
#'   In matrix form, the gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{j = 1, \dots,} \code{n_cat}. The numeric vector \code{delta_v} is
#'   obtained by concatenating the delta matrix, i.e. \code{delta_v <- c(delta_matrix)}.
#' @param X A numeric design matrix for the first-stage observed mechanism.
#' @param V A numeric design matrix for the second-stage observed mechanism.
#' @param obs_Ystar_matrix A numeric matrix of indicator variables (0, 1) for the first-stage observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param obs_Ytilde_matrix A numeric matrix of indicator variables (0, 1) for the second-stage observed
#'   outcome \eqn{\tilde{Y}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param sample_size Integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X} or \code{V}.
#' @param n_cat The number of categorical values that the first- and second-stage outcomes,
#'   \eqn{Y^*} and \eqn{\tilde{Y}}, can take.
#'
#' @return \code{naive_loglik_2stage} returns the negative value of the observed data log-likelihood function,
#'   \eqn{ \sum_{i = 1}^N \Bigl[ \sum_{k = 1}^2 \sum_{k = 1}^2 \sum_{\ell = 1}^2 y^*_{ik} \tilde{y_i} \text{log} \{ P(\tilde{Y}_{i} = \ell, Y^*_i = k | x_i, v_i) \}\Bigr]},
#'   at the provided inputs.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
naive_loglik_2stage <- function(param_current,
                                X, V,
                                obs_Ystar_matrix, obs_Ytilde_matrix,
                                sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  delta_current = matrix(c(param_current[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(V)))]),
                         ncol = n_cat, byrow = FALSE)

  pi_k_obs = pi_compute(beta_current, X, sample_size, n_cat)
  pi_kl_obs = pistar_compute(delta_current, V, sample_size, n_cat)

  k1l1 = ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,1] == 1,
                1, 0)
  k1l2 = ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,2] == 1,
                1, 0)
  k2l1 = ifelse(obs_Ystar_matrix[,2] == 1 & obs_Ytilde_matrix[,1] == 1,
                1, 0)
  k2l2 = ifelse(obs_Ystar_matrix[,2] == 1 & obs_Ytilde_matrix[,2] == 1,
                1, 0)

  l1_index <- 1:sample_size
  l2_index <- (sample_size + 1):(2 * sample_size)

  k1l1_sum = k1l1 * log((pi_kl_obs[l1_index, 1] * pi_k_obs[,1]) +
                          (pi_kl_obs[l1_index, 2] * pi_k_obs[,2]))
  k1l2_sum = k1l2 * log((pi_kl_obs[l2_index, 1] * pi_k_obs[,1]) +
                          (pi_kl_obs[l2_index, 2] * pi_k_obs[,2]))
  k2l1_sum = k2l1 * log((pi_kl_obs[l1_index, 1] * pi_k_obs[,1]) +
                          (pi_kl_obs[l1_index, 2] * pi_k_obs[,2]))
  k2l2_sum = k2l2 * log((pi_kl_obs[l2_index, 1] * pi_k_obs[,1]) +
                          (pi_kl_obs[l2_index, 2] * pi_k_obs[,2]))

  loglik_return = -sum(k1l1_sum + k1l2_sum + k2l1_sum + k2l2_sum)

  return(loglik_return)
}
