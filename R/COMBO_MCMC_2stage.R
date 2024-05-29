#' MCMC Estimation of the Two-Stage Binary Outcome Misclassification Model
#'
#' Jointly estimate \eqn{\beta}, \eqn{\gamma^{(1)}}, and \eqn{\gamma^{(2)}} parameters from the true outcome
#' first-stage observation, and second-stage observation mechanisms, respectively,
#' in a two-stage binary outcome misclassification model.
#'
#' @param Ystar1 A numeric vector of indicator variables (1, 2) for the observed
#'   outcome \eqn{Y^{*(1)}}. The reference category is 2.
#' @param Ystar2 A numeric vector of indicator variables (1, 2) for the second-stage
#'   observed outcome \eqn{Y^{*(2)}}. There should be no \code{NA} terms. The
#'   reference category is 2.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept.
#' @param z1_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z1_matrix} should not contain an intercept.
#' @param z2_matrix A numeric matrix of covariates in the second-stage observation mechanism.
#'   \code{z2_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param prior A character string specifying the prior distribution for the
#'   \eqn{\beta}, \eqn{\gamma^{(1)}}, and \eqn{\gamma^{(2)}} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
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
#'   All matrices in the list should have dimensions \code{n_cat} X \code{dim_x}, and all
#'   elements in the \code{n_cat} row should be set to \code{NA}.
#' @param gamma1_prior_parameters A numeric list of prior distribution parameters
#'   for the \eqn{\gamma^{(1)}} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain an array of location, lower bound, mean, or shape parameters,
#'   respectively, for \eqn{\gamma^{(1)}} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain an array of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for \eqn{\gamma^{(1)}} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   an array of the degrees of freedom for \eqn{\gamma^{(1)}} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All arrays in the list should have dimensions \code{n_cat} X \code{n_cat} X \code{dim_z1},
#'   and all elements in the \code{n_cat} row should be set to \code{NA}.
#' @param gamma2_prior_parameters A numeric list of prior distribution parameters
#'   for the \eqn{\gamma^{(2)}} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain an array of location, lower bound, mean, or shape parameters,
#'   respectively, for \eqn{\gamma^{(2)}} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain an array of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for \eqn{\gamma^{(2)}} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   an array of the degrees of freedom for \eqn{\gamma^{(2)}} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All arrays in the list should have dimensions \code{n_cat} X
#'   \code{n_cat} X \code{n_cat} X \code{dim_z2},
#'   and all elements in the \code{n_cat} row should be set to \code{NA}.
#' @param naive_gamma2_prior_parameters A numeric list of prior distribution parameters
#'   for the naive model \eqn{\gamma^{(2)}} terms. For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the first element of the
#'   list should contain an array of location, lower bound, mean, or shape parameters,
#'   respectively, for naive \eqn{\gamma^{(2)}} terms.
#'   For prior distributions \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"}, the second element of the
#'   list should contain an array of shape, upper bound, standard deviation, or scale parameters,
#'   respectively, for naive \eqn{\gamma^{(2)}} terms.
#'   For prior distribution \code{"t"}, the third element of the list should contain
#'   an array of the degrees of freedom for naive \eqn{\gamma^{(2)}} terms.
#'   The third list element should be empty for all other prior distributions.
#'   All arrays in the list should have dimensions \code{n_cat} X \code{n_cat} X \code{dim_z2},
#'   and all elements in the \code{n_cat} row should be set to \code{NA}.
#'   Note that prior distributions for the naive \eqn{\beta} terms are inherted
#'   from the \code{beta_prior_parameters} argument.
#' @param number_MCMC_chains An integer specifying the number of MCMC chains to compute.
#'   The default is \code{4}.
#' @param MCMC_sample An integer specifying the number of MCMC samples to draw.
#'   The default is \code{2000}.
#' @param burn_in An integer specifying the number of MCMC samples to discard
#'   for the burn-in period. The default is \code{1000}.
#' @param display_progress A logical value specifying whether messages should be
#'   displayed during model compilation. The default is \code{TRUE}.
#'
#' @return \code{COMBO_MCMC_2stage} returns a list of the posterior samples and posterior
#'   means for both the binary outcome misclassification model and a naive logistic
#'   regression of the observed outcome, \code{Y*}, predicted by the matrix \code{x}.
#'   The list contains the following components:
#'   \item{posterior_sample_df}{A data frame containing three columns. The first
#'   column indicates the chain from which a sample is taken, from 1 to \code{number_MCMC_chains}.
#'   The second column specifies the parameter associated with a given row. \eqn{\beta}
#'   terms have dimensions \code{dim_x} X \code{n_cat}. The \eqn{\gamma^{(1)}} terms
#'   have dimensions \code{n_cat} X \code{n_cat} X \code{dim_z1}, where the first
#'   index specifies the first-stage observed outcome category and the second index specifies
#'   the true outcome category. The \eqn{\gamma^{(2)}} terms
#'   have dimensions \code{n_cat} X \code{n_cat} X \code{n_cat} X \code{dim_z2}, where the first
#'   index specifies the second-stage observed outcome category, the second index specifies
#'   the first-stage observed outcome category, and the third index specifies the true outcome category.
#'   The final column provides the MCMC sample.}
#'   \item{posterior_means_df}{A data frame containing three columns. The first
#'   column specifies the parameter associated with a given row. Parameters are
#'   indexed as in the \code{posterior_sample_df}. The second column provides
#'   the posterior mean computed across all chains and all samples. The final column
#'   provides the posterior median computed across all chains and all samples.}
#'   \item{naive_posterior_sample_df}{A data frame containing three columns.
#'   The first column indicates the chain from which a sample is taken, from
#'   1 to \code{number_MCMC_chains}. The second column specifies the parameter
#'   associated with a given row. Naive \eqn{\beta} terms have dimensions
#'   \code{dim_x} X \code{n_cat}. The final column provides the MCMC sample.}
#'   \item{naive_posterior_means_df}{A data frame containing three columns. The first
#'   column specifies the naive parameter associated with a given row. Parameters are
#'   indexed as in the \code{naive_posterior_sample_df}. The second column provides
#'   the posterior mean computed across all chains and all samples. The final column
#'   provides the posterior median computed across all chains and all samples.}
#'
#' @export
#'
#' @include model_picker_2stage.R
#' @include jags_picker.R
#' @include naive_model_picker.R
#' @include naive_jags_picker.R
#' @include pistar_by_chain.R
#' @include check_and_fix_chains.R
#' @include sum_every_n.R
#' @include sum_every_n1.R
#'
#' @importFrom stats rnorm rgamma rmultinom median
#' @importFrom rjags coda.samples jags.model
#' @importFrom dplyr select filter `%>%` mutate group_by ungroup summarise all_of
#' @importFrom tidyr gather
#'
#' @examples \donttest{
#'
#' # Helper functions
#' sum_every_n <- function(x, n){
#' vector_groups = split(x,
#'                       ceiling(seq_along(x) / n))
#' sum_x = Reduce(`+`, vector_groups)
#'
#' return(sum_x)
#' }
#'
#' sum_every_n1 <- function(x, n){
#' vector_groups = split(x,
#'                       ceiling(seq_along(x) / n))
#' sum_x = Reduce(`+`, vector_groups) + 1
#'
#' return(sum_x)
#' }
#'
#' # Example
#'
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
#' x_matrix = matrix(rnorm(n, x_mu, x_sigma), ncol = 1)
#' X = matrix(c(rep(1, n), x_matrix[,1]), ncol = 2, byrow = FALSE)
#' z1_matrix = matrix(rgamma(n, z1_shape), ncol = 1)
#' Z1 = matrix(c(rep(1, n), z1_matrix[,1]), ncol = 2, byrow = FALSE)
#' z2_matrix = matrix(rgamma(n, z2_shape), ncol = 1)
#' Z2 = matrix(c(rep(1, n), z2_matrix[,1]), ncol = 2, byrow = FALSE)
#'
#' exp_xb = exp(X %*% true_beta)
#' pi_result = exp_xb[,1] / (exp_xb[,1] + 1)
#' pi_matrix = matrix(c(pi_result, 1 - pi_result), ncol = 2, byrow = FALSE)
#'
#' true_Y <- rep(NA, n)
#' for(i in 1:n){
#'     true_Y[i] = which(stats::rmultinom(1, 1, pi_matrix[i,]) == 1)
#' }
#'
#' exp_z1g1 = exp(Z1 %*% true_gamma1)
#' pistar1_denominator = matrix(c(1 + exp_z1g1[,1], 1 + exp_z1g1[,2]),
#'                              ncol = 2, byrow = FALSE)
#' pistar1_result = exp_z1g1 / pistar1_denominator
#'
#' pistar1_matrix = matrix(c(pistar1_result[,1], 1 - pistar1_result[,1],
#'                           pistar1_result[,2], 1 - pistar1_result[,2]),
#'                         ncol = 2, byrow = FALSE)
#'
#'
#' obs_Y1 <- rep(NA, n)
#' for(i in 1:n){
#'     true_j = true_Y[i]
#'     obs_Y1[i] = which(rmultinom(1, 1,
#'                       pistar1_matrix[c(i, n + i),
#'                                      true_j]) == 1)
#'  }
#'
#' Ystar1 <- obs_Y1
#'
#' exp_z2g2_1 = exp(Z2 %*% true_gamma2[,,1])
#' exp_z2g2_2 = exp(Z2 %*% true_gamma2[,,2])
#'
#' pi_denominator1 = apply(exp_z2g2_1, FUN = sum_every_n1, n, MARGIN = 2)
#' pi_result1 = exp_z2g2_1 / rbind(pi_denominator1)
#'
#' pi_denominator2 = apply(exp_z2g2_2, FUN = sum_every_n1, n, MARGIN = 2)
#' pi_result2 = exp_z2g2_2 / rbind(pi_denominator2)
#'
#' pistar2_matrix1 = rbind(pi_result1,
#'                         1 - apply(pi_result1,
#'                                   FUN = sum_every_n, n = n,
#'                                   MARGIN = 2))
#'
#' pistar2_matrix2 = rbind(pi_result2,
#'                         1 - apply(pi_result2,
#'                                   FUN = sum_every_n, n = n,
#'                                   MARGIN = 2))
#'
#' pistar2_array = array(c(pistar2_matrix1, pistar2_matrix2),
#'                    dim = c(dim(pistar2_matrix1), 2))
#'
#' obs_Y2 <- rep(NA, n)
#' for(i in 1:n){
#'     true_j = true_Y[i]
#'     obs_k = Ystar1[i]
#'     obs_Y2[i] = which(rmultinom(1, 1,
#'                                   pistar2_array[c(i,n+ i),
#'                                                 obs_k, true_j]) == 1)
#' }
#'
#' Ystar2 <- obs_Y2
#'
#' unif_lower_beta <- matrix(c(-5, -5, NA, NA), nrow = 2, byrow = TRUE)
#' unif_upper_beta <- matrix(c(5, 5, NA, NA), nrow = 2, byrow = TRUE)
#'
#' unif_lower_gamma1 <- array(data = c(-5, NA, -5, NA, -5, NA, -5, NA),
#'                           dim = c(2,2,2))
#' unif_upper_gamma1 <- array(data = c(5, NA, 5, NA, 5, NA, 5, NA),
#'                           dim = c(2,2,2))
#'
#' unif_upper_gamma2 <- array(rep(c(5, NA), 8), dim = c(2,2,2,2))
#' unif_lower_gamma2 <- array(rep(c(-5, NA), 8), dim = c(2,2,2,2))
#'
#' unif_lower_naive_gamma2 <- array(data = c(-5, NA, -5, NA, -5, NA, -5, NA),
#'                                 dim = c(2,2,2))
#' unif_upper_naive_gamma2 <- array(data = c(5, NA, 5, NA, 5, NA, 5, NA),
#'                                 dim = c(2,2,2))
#'
#' beta_prior_parameters <- list(lower = unif_lower_beta, upper = unif_upper_beta)
#' gamma1_prior_parameters <- list(lower = unif_lower_gamma1, upper = unif_upper_gamma1)
#' gamma2_prior_parameters <- list(lower = unif_lower_gamma2, upper = unif_upper_gamma2)
#' naive_gamma2_prior_parameters <- list(lower = unif_lower_naive_gamma2,
#'                                       upper = unif_upper_naive_gamma2)
#'
#' MCMC_results <- COMBO_MCMC_2stage(Ystar1, Ystar2,
#'                                   x_matrix = x_matrix, z1_matrix = z1_matrix,
#'                                   z2_matrix = z2_matrix,
#'                                   prior = "uniform",
#'                                   beta_prior_parameters = beta_prior_parameters,
#'                                   gamma1_prior_parameters = gamma1_prior_parameters,
#'                                   gamma2_prior_parameters = gamma2_prior_parameters,
#'                                   naive_gamma2_prior_parameters = naive_gamma2_prior_parameters,
#'                                   number_MCMC_chains = 2,
#'                                   MCMC_sample = 200, burn_in = 100)
#' MCMC_results$posterior_means_df}
COMBO_MCMC_2stage <- function(Ystar1, Ystar2, x_matrix, z1_matrix, z2_matrix,
                              prior,
                              beta_prior_parameters,
                              gamma1_prior_parameters,
                              gamma2_prior_parameters,
                              naive_gamma2_prior_parameters,
                              number_MCMC_chains = 4,
                              MCMC_sample = 2000,
                              burn_in = 1000,
                              display_progress = TRUE){

  Ystar <- Ystar1
  Ytilde <- Ystar2

  x <- x_matrix
  z <- z1_matrix
  v <- z2_matrix

  gamma_prior_parameters <- gamma1_prior_parameters
  delta_prior_parameters <- gamma2_prior_parameters
  naive_delta_prior_parameters <- naive_gamma2_prior_parameters

  # Define global variables to make the "NOTES" happy.
  chain_number <- NULL
  parameter_name <- NULL

  if (!is.numeric(Ystar) || !is.vector(Ystar))
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ystar))) != 0)
    stop("'Ystar' must be coded 1/2, where the reference category is 2.")

  sample_size = length(Ystar)
  n_cat = 2

  # FIX THIS! DIMENSIONS DON'T WORK FOR x, z WITH MORE THAN ONE COL
  X = cbind(matrix(1, nrow = sample_size, ncol = 1), x)
  Z = cbind(matrix(1, nrow = sample_size, ncol = 1), z)
  V = cbind(matrix(1, nrow = sample_size, ncol = 1), v)

  dim_x = ncol(X)
  dim_z = ncol(Z)
  dim_v = ncol(V)

  modelstring = model_picker_2stage(prior)

  temp_model_file = tempfile()
  tmps = file(temp_model_file, "w")
  cat(modelstring, file = tmps)
  close(tmps)

  jags <- jags_picker_2stage(prior, sample_size, dim_x, dim_z, dim_v, n_cat,
                             Ystar, Ytilde, X, Z, V,
                             beta_prior_parameters, gamma_prior_parameters,
                             delta_prior_parameters,
                             number_MCMC_chains,
                             model_file = temp_model_file,
                             display_progress = display_progress)

  display_progress_bar <- ifelse(display_progress == TRUE, "text", "none")

  posterior_sample = coda.samples(jags,
                                  c('beta', 'gamma', 'delta'),
                                  MCMC_sample,
                                  progress.bar = display_progress_bar)

  pistarjj = pistar_by_chain(n_chains = number_MCMC_chains,
                             chains_list = posterior_sample,
                             Z = Z, n = sample_size, n_cat = n_cat)

  pitildejjj = pitilde_by_chain(n_chains = number_MCMC_chains,
                                chains_list = posterior_sample,
                                V = V, n = sample_size, n_cat = n_cat)

  posterior_sample_fixed = check_and_fix_chains_2stage(n_chains = number_MCMC_chains,
                                                       chains_list = posterior_sample,
                                                       pistarjj_matrix = pistarjj,
                                                       pitildejjj_matrix = pitildejjj,
                                                       dim_x, dim_z, dim_v,
                                                       n_cat)

  posterior_sample_df <- do.call(rbind.data.frame, posterior_sample_fixed)
  posterior_sample_df$chain_number <- rep(1:number_MCMC_chains, each = MCMC_sample)
  posterior_sample_df$sample <- rep(1:MCMC_sample, number_MCMC_chains)

  ##########################################
  naive_modelstring = naive_model_picker_2stage(prior)

  naive_temp_model_file = tempfile()
  tmps = file(naive_temp_model_file, "w")
  cat(naive_modelstring, file = tmps)
  close(tmps)
  naive_jags <- naive_jags_picker_2stage(prior, sample_size, dim_x, dim_v, n_cat,
                                         Ystar, Ytilde, X, V,
                                         beta_prior_parameters,
                                         naive_delta_prior_parameters,
                                         number_MCMC_chains,
                                         naive_model_file = naive_temp_model_file,
                                         display_progress = display_progress)

  naive_posterior_sample = coda.samples(naive_jags,
                                        c('beta', 'delta'),
                                        MCMC_sample,
                                        progress.bar = display_progress_bar)

  naive_posterior_sample_df <- do.call(rbind.data.frame, naive_posterior_sample)
  naive_posterior_sample_df$chain_number <- rep(1:number_MCMC_chains, each = MCMC_sample)
  naive_posterior_sample_df$sample <- rep(1:MCMC_sample, number_MCMC_chains)
  ###########################################

  beta_names <- paste0("beta[1,", 1:dim_x, "]")
  gamma_names <- paste0("gamma1[1,", rep(1:n_cat, dim_z), ",", rep(1:dim_z, each = n_cat), "]")
  delta_names <- paste0("gamma2[1,",
                        rep(1:n_cat, dim_v*dim_v), ",",
                        rep(rep(1:n_cat, each = dim_v), dim_v), ",",
                        rep(1:dim_v, each = n_cat * n_cat), "]")
  naive_delta_names <- paste0("gamma2[1,", rep(1:n_cat, dim_v), ",", rep(1:dim_z, each = n_cat), "]")

  naive_posterior_sample_burn <- naive_posterior_sample_df %>%
    dplyr::select(dplyr::all_of(beta_names), dplyr::all_of(naive_delta_names),
                  chain_number, sample) %>%
    dplyr::filter(sample > burn_in) %>%
    tidyr::gather(parameter_name, sample, beta_names[1]:naive_delta_names[length(naive_delta_names)],
                  factor_key = TRUE) %>%
    dplyr::mutate(parameter_name = paste0("naive_", parameter_name))

  naive_posterior_means <- naive_posterior_sample_burn %>%
    dplyr::group_by(parameter_name) %>%
    dplyr::summarise(posterior_mean = mean(sample),
                     posterior_median = stats::median(sample)) %>%
    dplyr::ungroup()

  posterior_sample_burn <- posterior_sample_df %>%
    dplyr::select(dplyr::all_of(beta_names), dplyr::all_of(gamma_names),
                  dplyr::all_of(delta_names),
                  chain_number, sample) %>%
    dplyr::filter(sample > burn_in) %>%
    tidyr::gather(parameter_name, sample,
                  beta_names[1]:delta_names[length(delta_names)], factor_key = TRUE)

  posterior_means <- posterior_sample_burn %>%
    dplyr::group_by(parameter_name) %>%
    dplyr::summarise(posterior_mean = mean(sample),
                     posterior_median = stats::median(sample)) %>%
    dplyr::ungroup()


  results = list(posterior_sample_df = posterior_sample_burn,
                 posterior_means_df = posterior_means,
                 naive_posterior_sample_df = naive_posterior_sample_burn,
                 naive_posterior_means_df = naive_posterior_means)

  return(results)
}
