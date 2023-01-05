#' Set up a Naive Logistic Regression \code{jags.model} Object for a Given Prior
#'
#' @param prior  character string specifying the prior distribution for the naive
#'   \eqn{\beta} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#' @param sample_size An integer value specifying the number of observations in the sample.
#' @param dim_x An integer specifying the number of columns of the design matrix of the true outcome mechanism, \code{X}.
#' @param n_cat An integer specifying the number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#' @param Ystar A numeric vector of indicator variables (1, 2) for the observed
#'   outcome \code{Y*}. The reference category is 2.
#' @param X A numeric design matrix for the true outcome mechanism.
#' @param beta_prior_parameters A numeric list of prior distribution parameters
#'   for the \eqn{\beta} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain a matrix of location, lower bound, mean, or shape parameters,
#'   respectively, for \eqn{\beta} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain a matrix of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for \eqn{\beta} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   a matrix of the degrees of freedom for \eqn{\beta} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All matrices in the list should have dimensions \code{dim_x} X \code{n_cat}, and all
#'   elements in the \code{n_cat} column should be set to \code{NA}.
#' @param number_MCMC_chains An integer specifying the number of MCMC chains to compute.
#' @param naive_model_file A .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#'
#' @return \code{naive_jags_picker} returns a \code{jags.model} object for a naive
#'   logistic regression model predicting the potentially misclassified \code{Y*}
#'   from the predictor matrix \code{x}. The object includes the specified
#'   prior distribution, model, number of chains, and data.
#'
#' @importFrom stats rnorm rmultinom
#' @importFrom rjags jags.model
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 100
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
#' obs_Y = my_data[["obs_Y"]]
#' X = my_data[["x_design_matrix"]]
#'
#' unif_lower_beta <- matrix(c(-5, -5, NA, NA), nrow = 2, byrow = TRUE)
#' unif_upper_beta <- matrix(c(5, 5, NA, NA), nrow = 2, byrow = TRUE)
#'
#' beta_prior_parameters <- list(lower = unif_lower_beta, upper = unif_upper_beta)
#'
#' modelstring = naive_model_picker(prior = "uniform")
#' temp_model_file = tempfile()
#' tmps = file(temp_model_file, "w")
#' cat(modelstring, file = tmps)
#' close(tmps)
#'
#' jags_model_object <- naive_jags_picker(prior = "uniform",
#'                                        sample_size = n,
#'                                        dim_x = ncol(X),
#'                                        n_cat = 2,
#'                                        Ystar = obs_Y, X = X,
#'                                        beta_prior_parameters = beta_prior_parameters,
#'                                        number_MCMC_chains = 1,
#'                                        naive_model_file = temp_model_file)
#' }

naive_jags_picker_2stage <- function(prior, sample_size, dim_x, dim_v, n_cat,
                                     Ystar, Ytilde, X, V,
                                     beta_prior_parameters,
                                     delta_prior_parameters,
                                     number_MCMC_chains,
                                     naive_model_file){
  if (prior == "t") {
    jags_object <- jags.model(
      naive_model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_v = dim_v,
                  n_cat = n_cat,
                  Y_star = Ystar,
                  Y_tilde = Ytilde,
                  x = X,
                  v = V,
                  t_mu_beta = beta_prior_parameters[[1]],
                  t_tau_beta = beta_prior_parameters[[2]],
                  t_df_beta = beta_prior_parameters[[3]],
                  t_mu_delta = delta_prior_parameters[[1]],
                  t_tau_delta = delta_prior_parameters[[2]],
                  t_df_delta = delta_prior_parameters[[3]]),
      n.chains = number_MCMC_chains)
  } else if (prior == "uniform") {
    jags_object <- jags.model(
      naive_model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_v = dim_v,
                  n_cat = n_cat,
                  Y_star = Ystar,
                  Y_tilde = Ytilde,
                  x = X,
                  v = V,
                  unif_l_beta = beta_prior_parameters[[1]],
                  unif_u_beta = beta_prior_parameters[[2]],
                  unif_l_delta = delta_prior_parameters[[1]],
                  unif_u_delta = delta_prior_parameters[[2]]),
      n.chains = number_MCMC_chains)
  } else if (prior == "normal") {
    jags_object <- jags.model(
      naive_model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X,
                  normal_mu_beta = beta_prior_parameters[[1]],
                  normal_sigma_beta = beta_prior_parameters[[2]]),
      n.chains = number_MCMC_chains)
  } else if (prior == "dexp") {
    jags_object <- jags.model(
      naive_model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X,
                  dexp_mu_beta = beta_prior_parameters[[1]],
                  dexp_b_beta = beta_prior_parameters[[2]]),
      n.chains = number_MCMC_chains)
  } else { print("Please select a model.")}

  return(jags_object)
}
