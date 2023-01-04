#' Compute the Mean Conditional Probability of Correct Classification, by True Outcome Across all Subjects for each MCMC Chain
#'
#' @param n_chains An integer specifying the number of MCMC chains to compute over.
#' @param chains_list A numeric list containing the samples from \code{n_chains}
#'   MCMC chains.
#' @param Z A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{pistar_by_chain} returns a numeric matrix of the average
#'   conditional probability \eqn{P(Y^* = j | Y = j, Z)} across all subjects for
#'   each MCMC chain. Rows of the matrix correspond to MCMC chains, up to \code{n_chains}.
#'   The first column contains the conditional probability \eqn{P(Y^* = 1 | Y = 1, Z)}.
#'   The second column contains the conditional probability \eqn{P(Y^* = 2 | Y = 2, Z)}.
#'
#' @include pistar_compute_for_chains.R
#' @include mean_pistarjj_compute.R
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
#' example_chain1 <- matrix(c(rnorm(n, mean = gamma[1,1]),
#'                            rnorm(n, mean = gamma[2,1]),
#'                            rnorm(n, mean = gamma[1,2]),
#'                            rnorm(n, mean = gamma[2,2])),
#'                            nrow = n, byrow = FALSE)
#' colnames(example_chain1) <- c("gamma[1,1,1]", "gamma[1,1,2]",
#'                               "gamma[1,2,1]", "gamma[1,2,2]")
#' example_chain2 <- matrix(c(rnorm(n, mean = gamma[1,1]),
#'                            rnorm(n, mean = gamma[2,1]),
#'                            rnorm(n, mean = gamma[1,2]),
#'                            rnorm(n, mean = gamma[2,2])),
#'                            nrow = n, byrow = FALSE)
#' colnames(example_chain2) <- c("gamma[1,1,1]", "gamma[1,1,2]",
#'                               "gamma[1,2,1]", "gamma[1,2,2]")
#' chains_list <- list(example_chain1, example_chain2)
#' pistar_by_chain(n_chains = 2, chains_list = chains_list,
#'                 Z = Z, n = n, n_cat = 2)
#' }
pitilde_by_chain <- function(n_chains, chains_list, V, n, n_cat){

  colmeans_by_chain = lapply(chains_list, colMeans)

  pitilde_results = matrix(NA, nrow = n_chains, ncol = n_cat)
  for(i in 1:n_chains){
    chain_means_i = colmeans_by_chain[[i]]
    pitilde_array_i = pitilde_compute_for_chains(chain_means_i, V, n, n_cat)
    pitilde_results[i, 1] = mean(pitilde_array_i[1:n, 1, 1])
    pitilde_results[i, 2] = mean(pitilde_array_i[(n + 1):(2 * n), 2, 2])
  }

  return(pitilde_results)
}
