#' Compute Conditional Probability of Each Observed Outcome Given Each True Outcome for a given MCMC Chain, for Every Subject
#'
#' @param chain_colMeans A numeric vector containing the posterior means for all
#'   sampled parameters in a given MCMC chain. \code{chain_colMeans} must be a named
#'   object (i.e. each parameter must be named as \code{gamma[k,j,p]}).
#' @param Z A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{pistar_compute_for_chains} returns a matrix of conditional probabilities,
#'   \eqn{P(Y_i^* = k | Y_i = j, Z_i) = \frac{\text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}{1 + \text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}}
#'   for each of the \eqn{i = 1, \dots,} \code{n} subjects. Rows of the matrix
#'   correspond to each subject and observed outcome. Specifically, the probability
#'   for subject \eqn{i} and observed category $0$ occurs at row \eqn{i}. The probability
#'   for subject \eqn{i} and observed category $1$ occurs at row \eqn{i +} \code{n}.
#'   Columns of the matrix correspond to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include sum_every_n.R
#'
#' @importFrom stats rnorm
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 100
#' ones <- rep(1, n)
#' z <- rnorm(n)
#' Z <- matrix(c(ones, z), nrow = n, byrow = FALSE)
#' gamma <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = FALSE)
#' example_MCMC_chain <- matrix(c(rnorm(n, mean = gamma[1,1]),
#'                                rnorm(n, mean = gamma[2,1]),
#'                                rnorm(n, mean = gamma[1,2]),
#'                                rnorm(n, mean = gamma[2,2])),
#'                              nrow = n, byrow = FALSE)
#' colnames(example_MCMC_chain) <- c("gamma[1,1,1]", "gamma[1,1,2]",
#'                                   "gamma[1,2,1]", "gamma[1,2,2]")
#' chain_colMeans <- colMeans(example_MCMC_chain)
#' conditional_probabilities <- pistar_compute_for_chains(chain_colMeans, Z, n, n_cat = 2)
#' head(conditional_probabilities)
#' }
pistar_compute_for_chains <- function(chain_colMeans, Z, n, n_cat){

  dim_z = ncol(Z)
  gamma_names <- paste0("gamma[1,", rep(1:n_cat, dim_z), ",", rep(1:dim_z, each = n_cat), "]")
  gamma = matrix(chain_colMeans[gamma_names],
                 nrow = 2, byrow = TRUE)

  exp_zg_nobaseline = exp(Z %*% gamma)
  exp_zg_baseline = matrix(1, nrow = n, ncol = n_cat)
  exp_zg = rbind(exp_zg_nobaseline, exp_zg_baseline)

  pi_denominator = apply(exp_zg, FUN = sum_every_n, n, MARGIN = 2)
  pi_result = exp_zg / do.call("rbind", rep(list(pi_denominator),
                                            n_cat))

  return(pi_result)
}
