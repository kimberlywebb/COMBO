#' Set up a Binary Outcome Misclassification \code{jags.model} Object for a Given Prior
#'
#' @param prior A character string specifying the prior distribution for the
#'   \eqn{\beta} and \eqn{\gamma} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#' @param sample_size An integer value specifying the number of observations in the sample.
#' @param dim_x An integer specifying the number of columns of the design matrix of the true outcome mechanism, \code{X}.
#' @param dim_z An integer specifying the number of columns of the design matrix of the observation mechanism, \code{Z}.
#' @param n_cat An integer specifying the number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#' @param Ystar A numeric vector of indicator variables (1, 2) for the observed
#'   outcome \code{Y*}. The reference category is 2.
#' @param X A numeric design matrix for the true outcome mechanism.
#' @param Z A numeric design matrix for the observation mechanism.
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
#' @param gamma_prior_parameters A numeric list of prior distribution parameters
#'   for the \eqn{\gamma} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain an array of location, lower bound, mean, or shape parameters,
#'   respectively, for \eqn{\gamma} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain an array of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for \eqn{\gamma} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   an array of the degrees of freedom for \eqn{\gamma} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All arrays in the list should have dimensions \code{n_cat} X \code{n_cat} X \code{dim_z},
#'   and all elements in the \code{n_cat} row should be set to \code{NA}.
#' @param number_MCMC_chains An integer specifying the number of MCMC chains to compute.
#' @param model_file A .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#' @param display_progress A logical value specifying whether messages should be
#'   displayed during model compilation. The default is \code{TRUE}.
#'
#' @return \code{jags_picker} returns a \code{jags.model} object for a binary
#'   outcome misclassification model. The object includes the specified
#'   prior distribution, model, number of chains, and data.
#'
#' @importFrom stats rnorm rmultinom optim
#' @importFrom rjags jags.model
#'
jags_picker <- function(prior, sample_size, dim_x, dim_z, n_cat,
                        Ystar, X, Z,
                        beta_prior_parameters, gamma_prior_parameters,
                        number_MCMC_chains,
                        model_file, display_progress = TRUE){

  quiet_argument <- !display_progress

  if (prior == "t") {
    jags_object <- jags.model(
      model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_z = dim_z,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X, z = Z,
                  t_mu_beta = beta_prior_parameters[[1]],
                  t_tau_beta = beta_prior_parameters[[2]],
                  t_df_beta = beta_prior_parameters[[3]],
                  t_mu_gamma = gamma_prior_parameters[[1]],
                  t_tau_gamma = gamma_prior_parameters[[2]],
                  t_df_gamma = gamma_prior_parameters[[3]]),
      n.chains = number_MCMC_chains,
      quiet = quiet_argument)
  } else if (prior == "uniform") {
    jags_object <- jags.model(
      model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_z = dim_z,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X, z = Z,
                  unif_l_beta = beta_prior_parameters[[1]],
                  unif_u_beta = beta_prior_parameters[[2]],
                  unif_l_gamma = gamma_prior_parameters[[1]],
                  unif_u_gamma = gamma_prior_parameters[[2]]),
      n.chains = number_MCMC_chains,
      quiet = quiet_argument)
  } else if (prior == "normal") {
    jags_object <- jags.model(
      model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_z = dim_z,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X, z = Z,
                  normal_mu_beta = beta_prior_parameters[[1]],
                  normal_sigma_beta = beta_prior_parameters[[2]],
                  normal_mu_gamma = gamma_prior_parameters[[1]],
                  normal_sigma_gamma = gamma_prior_parameters[[2]]),
      n.chains = number_MCMC_chains,
      quiet = quiet_argument)
  } else if (prior == "dexp") {
    jags_object <- jags.model(
      model_file,
      data = list(sample_size = sample_size,
                  dim_x = dim_x,
                  dim_z = dim_z,
                  n_cat = n_cat,
                  obs_Y = Ystar,
                  x = X, z = Z,
                  dexp_mu_beta = beta_prior_parameters[[1]],
                  dexp_b_beta = beta_prior_parameters[[2]],
                  dexp_mu_gamma = gamma_prior_parameters[[1]],
                  dexp_b_gamma = gamma_prior_parameters[[2]]),
      n.chains = number_MCMC_chains,
      quiet = quiet_argument)
  } else { stop("Please select a model.")}

  return(jags_object)
}
