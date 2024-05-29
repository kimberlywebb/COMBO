#' Generate data to use in two-stage COMBO Functions
#'
#' @param sample_size An integer specifying the sample size of the generated data set.
#' @param x_mu A numeric value specifying the mean of \code{x} predictors
#' generated from a Normal distribution.
#' @param x_sigma A positive numeric value specifying the standard deviation of
#'   \code{x} predictors generated from a Normal distribution.
#' @param z1_shape A positive numeric value specifying the shape parameter of
#'   \code{z1} predictors generated from a Gamma distribution.
#' @param z2_shape A positive numeric value specifying the shape parameter of
#'   \code{z2} predictors generated from a Gamma distribution.
#' @param beta A column matrix of \eqn{\beta} parameter values (intercept, slope)
#'   to generate data under in the true outcome mechanism.
#' @param gamma1 A numeric matrix of \eqn{\gamma^{(1)}} parameters
#'   to generate data under in the first-stage observation mechanism.
#'   In matrix form, the \code{gamma1} matrix rows correspond to intercept (row 1)
#'   and slope (row 2) terms. The \code{gamma1} parameter matrix columns correspond to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}.
#' @param gamma2 A numeric array of \eqn{\gamma^{(2)}} parameters to generate data under
#'   the second-stage observation mechanism. In array form, the \code{gamma2} matrix rows
#'   correspond to intercept (row 1) and slope (row 2) terms. The matrix columns correspond
#'   to first-stage observed outcome categories. The third dimension of the \code{gamma2}
#'   array is indexed by the true outcome categories.
#'
#' @return \code{COMBO_data_2stage} returns a list of generated data elements:
#'   \item{obs_Ystar1}{A vector of first-stage observed outcomes.}
#'   \item{obs_Ystar2}{A vector of second-stage observed outcomes.}
#'   \item{true_Y}{A vector of true outcomes.}
#'   \item{obs_Ystar1_matrix}{A numeric matrix of indicator variables (0, 1) for the first-stage observed
#'   outcome \eqn{Y^{*(1)}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row contains
#'   exactly one 0 entry and exactly one 1 entry.}
#'   \item{obs_Ystar2_matrix}{A numeric matrix of indicator variables (0, 1) for the second-stage observed
#'   outcome \eqn{Y^{*(2)}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row contains
#'   exactly one 0 entry and exactly one 1 entry.}
#'   \item{x}{A vector of generated predictor values in the true outcome
#'   mechanism, from the Normal distribution.}
#'   \item{z1}{A vector of generated predictor values in the first-stage observation
#'   mechanism from the Gamma distribution.}
#'   \item{z2}{A vector of generated predictor values in the second-stage observation
#'   mechanism from the Gamma distribution.}
#'   \item{x_design_matrix}{The design matrix for the \code{x} predictor.}
#'   \item{z1_design_matrix}{The design matrix for the \code{z1} predictor.}
#'   \item{z2_design_matrix}{The design matrix for the \code{z2} predictor.}
#'
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include pitilde_compute.R
#'
#' @examples
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z1_shape <- 1
#' z2_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma1 <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#' true_gamma2 <- array(c(1.5, 1, .5, .5, -.5, 0, -1, -1), dim = c(2, 2, 2))
#'
#' my_data <- COMBO_data_2stage(sample_size = n,
#'                              x_mu = x_mu, x_sigma = x_sigma,
#'                              z1_shape = z1_shape, z2_shape = z2_shape,
#'                              beta = true_beta, gamma1 = true_gamma1, gamma2 = true_gamma2)
#' table(my_data[["obs_Ystar2"]], my_data[["obs_Ystar1"]], my_data[["true_Y"]])
COMBO_data_2stage <- function(sample_size,
                              x_mu, x_sigma,
                              z1_shape, z2_shape,
                              beta, gamma1, gamma2){

  z_shape <- z1_shape
  v_shape <- z2_shape
  gamma <- gamma1
  delta <- gamma2

  n_cat <- 2
  x <- rnorm(sample_size, x_mu, x_sigma)
  x_matrix <- matrix(c(rep(1, sample_size),
                       x),
                     nrow = sample_size, byrow = FALSE)

  z <- rgamma(sample_size, z_shape)
  z_matrix <- matrix(c(rep(1, sample_size),
                       z),
                     nrow = sample_size, byrow = FALSE)

  v <- rgamma(sample_size, v_shape)
  v_matrix <- matrix(c(rep(1, sample_size),
                       v),
                     nrow = sample_size, byrow = FALSE)

  pi_matrix <- pi_compute(beta, x_matrix, sample_size, n_cat)

  true_Y <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_Y[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }

  pistar_matrix <- pistar_compute(gamma, z_matrix, sample_size, n_cat)

  obs_Ystar <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_Y[i]
    obs_Ystar[i] = which(rmultinom(1, 1,
                                   pistar_matrix[c(i,sample_size + i),
                                                 true_j]) == 1)
  }

  obs_Ystar_reps <- matrix(rep(obs_Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ystar_matrix <- 1 * (obs_Ystar_reps == category_matrix)

  pitilde_array <- pitilde_compute(delta, v_matrix, sample_size, n_cat)

  obs_Ytilde <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_Y[i]
    obs_k = obs_Ystar[i]
    obs_Ytilde[i] = which(rmultinom(1, 1,
                                   pitilde_array[c(i,sample_size + i),
                                                 obs_k, true_j]) == 1)
  }

  obs_Ytilde_reps <- matrix(rep(obs_Ytilde, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ytilde_matrix <- 1 * (obs_Ytilde_reps == category_matrix)

  data_output <- list(obs_Ystar1 = obs_Ystar,
                      obs_Ystar2 = obs_Ytilde,
                      true_Y = true_Y,
                      obs_Ystar1_matrix = obs_Ystar_matrix,
                      obs_Ystar2_matrix = obs_Ytilde_matrix,
                      x = x,
                      z1 = z,
                      z2 = v,
                      x_design_matrix = x_matrix,
                      z1_design_matrix = z_matrix,
                      z2_design_matrix = v_matrix)

  return(data_output)

}
