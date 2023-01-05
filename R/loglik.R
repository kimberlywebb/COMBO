#' Expected Complete Data Log-Likelihood Function for Estimation of the Misclassification Model
#'
#' @param param_current A numeric vector of regression parameters, in the order
#'   \eqn{\beta, \gamma}. The \eqn{\gamma} vector is obtained from the matrix form.
#'   In matrix form, the gamma parameter matrix rows
#'   correspond to parameters for the \code{Y* = 1}
#'   observed outcome, with the dimensions of \code{Z}.
#'   In matrix form, the gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{j = 1, \dots,} \code{n_cat}. The numeric vector \code{gamma_v} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_v <- c(gamma_matrix)}.
#' @param obs_Y_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param X A numeric design matrix for the true outcome mechanism.
#' @param Z A numeric design matrix for the observation mechanism.
#' @param sample_size Integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X} or \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{loglik} returns the negative value of the expected log-likelihood function,
#'   \eqn{ Q = \sum_{i = 1}^N \Bigl[ \sum_{j = 1}^2 w_{ij} \text{log} \{ \pi_{ij} \} + \sum_{j = 1}^2 \sum_{k = 1}^2 w_{ij} y^*_{ik} \text{log} \{ \pi^*_{ikj} \}\Bigr]},
#'   at the provided inputs.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_j.R
#' @include q_beta_f.R
#' @include q_gamma_f.R
#' @include em_function.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#'
#' my_data <- COMBO_data(sample_size = n,
#'                       x_mu = x_mu, x_sigma = x_sigma,
#'                       z_shape = z_shape,
#'                       beta = true_beta, gamma = true_gamma)
#'
#' obs_Y_matrix = my_data[["obs_Y_matrix"]]
#' X = my_data[["x_design_matrix"]]
#' Z = my_data[["z_design_matrix"]]
#'
#' starting_values <- rnorm(6)
#'
#' loglik_value <- loglik(starting_values,
#'                        obs_Y_matrix = obs_Y_matrix,
#'                        X = X, Z = Z,
#'                        sample_size = n, n_cat = 2)
#' loglik_value
#' }
loglik <- function(param_current,
                   obs_Y_matrix, X, Z,
                   sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  gamma_current = matrix(c(param_current[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))]),
                         ncol = n_cat, byrow = FALSE)

  pi_terms_v = pi_compute(beta_current, X, sample_size, n_cat)
  pistar_terms_v = pistar_compute(gamma_current, Z, sample_size, n_cat)

  weights = w_j(obs_Y_matrix, pistar_terms_v, pi_terms_v, sample_size, n_cat)

  loglikelihood = sum(
    (q_beta_f(beta_current, X = X, w_mat = weights,
              sample_size = sample_size, n_cat = n_cat)) +
      (q_gamma_f(c(gamma_current), Z = Z,
                 obs_Y_matrix = obs_Y_matrix,
                 w_mat = weights,
                 sample_size = sample_size, n_cat = n_cat)))

  return(loglikelihood)
}
