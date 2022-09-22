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
#' @examples \dontrun{
#' set.seed(123)
#' n <- 1000
#' n_cat <- 2
#'
#' ones <- rep(1, n)
#' x <- rnorm(n)
#' X <- matrix(c(ones, x), nrow = n, byrow = FALSE)
#' z <- rgamma(n, shape = 1)
#' Z <- matrix(c(ones, z), nrow = n, byrow = FALSE)
#'
#' beta <- matrix(c(1, -2), ncol = 1)
#' gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#'
#' probabilities <- pi_compute(beta, X, n, n_cat = 2)
#' conditional_probabilities <- pistar_compute(gamma, Z, n, n_cat = 2)
#'
#' true_Y <- rep(NA, n)
#' for(i in 1:n){
#'   true_Y[i] = which(rmultinom(1, 1, probabilities[i,]) == 1)
#' }
#' pistar_matrix_labels <- rep(1:n_cat, each = n)
#' obs_Y <- rep(NA, n)
#' for(i in 1:n){
#'  true_j = true_Y[i]
#'  obs_Y[i] = which(rmultinom(1, 1,
#'                             conditional_probabilities[c(i, n + i), true_j]) == 1)
#'  }
#'
#' obs_Y_reps <- matrix(rep(obs_Y, n_cat), nrow = n, byrow = FALSE)
#' category_matrix <- matrix(rep(1:n_cat, each = n), nrow = n, byrow = FALSE)
#' obs_Y_matrix <- 1 * (obs_Y_reps == category_matrix)
#'
#' e_step_weights <- w_j(ystar_matrix = obs_Y_matrix,
#'                       pistar_matrix = conditional_probabilities,
#'                       pi_matrix = probabilities,
#'                       sample_size = n, n_cat = n_cat)
#'
#' start_values <- c(1, 1)
#'
#' estimate_beta <- optim(par = start_values, fn = q_beta_f,
#'                        X = X, w_mat = e_step_weights,
#'                        sample_size = n, n_cat = n_cat,
#'                        method = "BFGS",
#'                        control = list(maxit = 5))
#' estimate_beta$par
#'
#' q_beta_f(beta = estimate_beta$par,
#'          X = X, w_mat = e_step_weights,
#'          sample_size = n, n_cat = n_cat)
#' }
q_beta_f <- function(beta, X, w_mat,
                     sample_size, n_cat){

  beta_mat = matrix(beta, ncol = 1)
  pi_terms = pi_compute(beta_mat, X, n = sample_size, n_cat)
  pi_terms = ifelse(pi_terms == 0, pi_terms + 0.00001, pi_terms)

  w_pi = c(w_mat) * log(c(pi_terms))
  result = -sum(sum_every_n(w_pi, n = sample_size))
  return(result)
}
