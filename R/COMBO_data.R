#' Generate Data to use in COMBO Functions
#'
#' @param sample_size An integer specifying the sample size of the generated data set.
#' @param x_mu A numeric value specifying the mean of \code{x} predictors
#' generated from a Normal distribution.
#' @param x_sigma A positive numeric value specifying the standard deviation of
#'   \code{x} predictors generated from a Normal distribution.
#' @param z_shape A positive numeric value specifying the shape parameter of
#'   \code{z} predictors generated from a Gamma distribution.
#' @param beta A column matrix of \eqn{\beta} parameter values (intercept, slope)
#'   to generate data under in the true outcome mechanism.
#' @param gamma A numeric matrix of \eqn{\gamma} parameters
#'   to generate data under in the observation mechanism.
#'   In matrix form, the \code{gamma} matrix rows correspond to intercept (row 1)
#'   and slope (row 2) terms. The gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}.
#'
#' @return \code{COMBO_data} returns a list of generated data elements:
#'   \item{obs_Y}{A vector of observed outcomes.}
#'   \item{true_Y}{A vector of true outcomes.}
#'   \item{obs_Y_matrix}{A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row contains
#'   exactly one 0 entry and exactly one 1 entry.}
#'   \item{x}{A vector of generated predictor values in the true outcome
#'   mechanism, from the Normal distribution.}
#'   \item{z}{A vector of generated predictor values in the observation
#'   mechanism from the Gamma distribution.}
#'   \item{x_design_matrix}{The design matrix for the \code{x} predictor.}
#'   \item{z_design_matrix}{The design matrix for the \code{z} predictor.}
#'
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
#' @examples
#' set.seed(123)
#' n <- 500
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
#' table(my_data[["obs_Y"]], my_data[["true_Y"]])
COMBO_data <- function(sample_size,
                       x_mu, x_sigma,
                       z_shape,
                       beta, gamma){

  n_cat <- 2
  x <- rnorm(sample_size, x_mu, x_sigma)
  x_matrix <- matrix(c(rep(1, sample_size),
                       x),
                     nrow = sample_size, byrow = FALSE)

  z <- rgamma(sample_size, z_shape)
  z_matrix <- matrix(c(rep(1, sample_size),
                       z),
                     nrow = sample_size, byrow = FALSE)

  pi_matrix <- pi_compute(beta, x_matrix, sample_size, n_cat)

  true_Y <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_Y[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }

  pistar_matrix <- pistar_compute(gamma, z_matrix, sample_size, n_cat)

  obs_Y <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_Y[i]
    obs_Y[i] = which(rmultinom(1, 1,
                               pistar_matrix[c(i,sample_size + i),
                                             true_j]) == 1)
  }

  obs_Y_reps <- matrix(rep(obs_Y, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Y_matrix <- 1 * (obs_Y_reps == category_matrix)

  data_output <- list(obs_Y = obs_Y,
                      true_Y = true_Y,
                      obs_Y_matrix = obs_Y_matrix,
                      x = x,
                      z = z,
                      x_design_matrix = x_matrix,
                      z_design_matrix = z_matrix)

  return(data_output)

}
