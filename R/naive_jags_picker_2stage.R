#' Set up a Naive Two-Stage Regression \code{jags.model} Object for a Given Prior
#'
#' @param prior  character string specifying the prior distribution for the naive
#'   \eqn{\beta} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#' @param sample_size An integer value specifying the number of observations in the sample.
#' @param dim_x An integer specifying the number of columns of the design matrix of the first-stage outcome mechanism, \code{X}.
#' @param dim_v An integer specifying the number of columns of the design matrix of the second-stage outcome mechanism, \code{V}.
#' @param n_cat An integer specifying the number of categorical values that
#' the observed outcomes can take.
#' @param Ystar A numeric vector of indicator variables (1, 2) for the first-stage observed
#'   outcome \code{Y*}. The reference category is 2.
#' @param Ytilde A numeric vector of indicator variables (1, 2) for the second-stage observed
#'   outcome \eqn{\tilde{Y}}. The reference category is 2.
#' @param X A numeric design matrix for the true outcome mechanism.
#' @param V A numeric design matrix for the second-stage outcome mechanism.
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
#' @param delta_prior_parameters A numeric list of prior distribution parameters
#'   for the naive \eqn{\delta} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain an array of location, lower bound, mean, or shape parameters,
#'   respectively, for \eqn{\delta} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain an array of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for \eqn{\delta} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   an array of the degrees of freedom for \eqn{\delta} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All arrays in the list should have dimensions \code{n_cat} X \code{n_cat} X \code{dim_v},
#'   and all elements in the \code{n_cat} row should be set to \code{NA}.
#' @param number_MCMC_chains An integer specifying the number of MCMC chains to compute.
#' @param naive_model_file A .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#'
#' @return \code{naive_jags_picker_2stage} returns a \code{jags.model} object for a naive
#'   two-stage regression model predicting the potentially misclassified \code{Y*}
#'   from the predictor matrix \code{x} and the potentially misclassified \eqn{\tilde{Y} | Y^*}
#'   from the predictor matrix \code{v}. The object includes the specified
#'   prior distribution, model, number of chains, and data.
#'
#' @importFrom stats rnorm rmultinom
#' @importFrom rjags jags.model
#'
#' @examples \dontrun{
#' set.seed(123)
#' sample_size <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#' v_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#' true_delta <- array(c(1.5, 1, .5, .5, -.5, 0, -1, -1), dim = c(2, 2, 2))
#'
#' my_data <- COMBO_data_2stage(sample_size = sample_size,
#'                              x_mu = x_mu, x_sigma = x_sigma,
#'                              z_shape = z_shape, v_shape = v_shape,
#'                              beta = true_beta, gamma = true_gamma, delta = true_delta)
#'
#' obs_Ystar = my_data[["obs_Ystar"]]
#' obs_Ytilde = my_data[["obs_Ytilde"]]
#' X = my_data[["x_design_matrix"]]
#' Z = my_data[["z_design_matrix"]]
#' V = my_data[["v_design_matrix"]]
#'
#' unif_lower_beta <- matrix(c(-5, -5, NA, NA), nrow = 2, byrow = TRUE)
#' unif_upper_beta <- matrix(c(5, 5, NA, NA), nrow = 2, byrow = TRUE)
#'
#' unif_lower_delta <- array(data = c(-5, NA, -5, NA, -5, NA, -5, NA),
#'                           dim = c(2,2,2))
#' unif_upper_delta <- array(data = c(5, NA, 5, NA, 5, NA, 5, NA),
#'                           dim = c(2,2,2))
#'
#' beta_prior_parameters <- list(lower = unif_lower_beta, upper = unif_upper_beta)
#' delta_prior_parameters <- list(lower = unif_lower_delta, upper = unif_upper_delta)
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
#'                                        dim_v = ncol(V),
#'                                        n_cat = 2,
#'                                        Ystar = obs_Ystar, Ytilde = obs_Ytilde,
#'                                        X = X, V = V,
#'                                        beta_prior_parameters = beta_prior_parameters,
#'                                        delta_prior_parameters = delta_prior_parameters,
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
                  dim_v = dim_v,
                  n_cat = n_cat,
                  Y_star = Ystar,
                  Y_tilde = Ytilde,
                  x = X,
                  v = V,
                  normal_mu_beta = beta_prior_parameters[[1]],
                  normal_sigma_beta = beta_prior_parameters[[2]],
                  normal_mu_delta = delta_prior_parameters[[1]],
                  normal_sigma_delta = delta_prior_parameters[[2]]),
      n.chains = number_MCMC_chains)
  } else if (prior == "dexp") {
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
                  dexp_mu_beta = beta_prior_parameters[[1]],
                  dexp_b_beta = beta_prior_parameters[[2]],
                  dexp_mu_delta = delta_prior_parameters[[1]],
                  dexp_b_delta = delta_prior_parameters[[2]]),
      n.chains = number_MCMC_chains)
  } else { print("Please select a model.")}

  return(jags_object)
}
