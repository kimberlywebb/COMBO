#' Compute Conditional Probability of Each Observed Outcome Given Each True Outcome for a given MCMC Chain, for Every Subject
#'
#' @param chain_colMeans A numeric vector containing the posterior means for all
#'   sampled parameters in a given MCMC chain. \code{chain_colMeans} must be a named
#'   object (i.e. each parameter must be named as \code{delta[l,k,j,p]}).
#' @param V A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{V}.
#' @param n_cat The number of categorical values that the true outcome, \eqn{Y},
#'   the first-stage observed outcome, \eqn{Y^*}, and the second-stage
#'   observed outcome, \eqn{\tilde{Y}},\ can take.
#'
#' @return \code{pitilde_compute_for_chains} returns a matrix of conditional probabilities,
#'   \eqn{P(\tilde{Y}_i = \ell | Y^*_i = k, Y_i = j, V_i) = \frac{\text{exp}\{\delta_{\ell kj0} + \delta_{\ell kjV} V_i\}}{1 + \text{exp}\{\delta_{\ell kj0} + \delta_{\ell kjV} V_i\}}}
#'   corresponding to each subject and observed outcome. Specifically, the probability
#'   for subject \eqn{i} and second-stage observed category $1$ occurs at row \eqn{i}. The probability
#'   for subject \eqn{i} and second-stage observed category $2$ occurs at row \eqn{i +} \code{n}.
#'   Columns of the matrix correspond to the first-stage outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'   The third dimension of the array corresponds to the true outcome categories,
#'   \eqn{j = 1, \dots,} \code{n_cat}.
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
pitilde_compute_for_chains <- function(chain_colMeans, V, n, n_cat){

  dim_v = ncol(V)
  delta_names <- paste0("delta[1,",
                        rep(1:n_cat, dim_v*dim_v), ",",
                        rep(rep(1:n_cat, each = dim_v), dim_v), ",",
                        rep(1:dim_v, each = n_cat * n_cat), "]")
  delta_j1_names <- which(substring(delta_names, 11, 11) == 1)
  delta_j2_names <- which(substring(delta_names, 11, 11) == 2)
  delta_j1 <- matrix(chain_colMeans[delta_names[delta_j1_names]],
                     ncol = 2, byrow = TRUE)
  delta_j2 <- matrix(chain_colMeans[delta_names[delta_j2_names]],
                     ncol = 2, byrow = TRUE)
  delta <- array(c(delta_j1, delta_j2), dim = c(dim_v, n_cat, n_cat))

  exp_vd1 = exp(V %*% delta[,,1])
  exp_vd2 = exp(V %*% delta[,,2])

  pi_denominator1 = apply(exp_vd1, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result1 = exp_vd1 / rbind(pi_denominator1)

  pi_denominator2 = apply(exp_vd2, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result2 = exp_vd2 / rbind(pi_denominator2)

  pitilde_matrix1 = rbind(pi_result1,
                          1 - apply(pi_result1,
                                    FUN = sum_every_n, n = n,
                                    MARGIN = 2))

  pitilde_matrix2 = rbind(pi_result2,
                          1 - apply(pi_result2,
                                    FUN = sum_every_n, n = n,
                                    MARGIN = 2))


  pitilde_array = array(c(pitilde_matrix1, pitilde_matrix2),
                        dim = c(dim(pitilde_matrix1), 2))

  return(pitilde_array)
}
