#' Compute the Mean Conditional Probability of Correct Classification, by True Outcome Across all Subjects for each MCMC Chain for a 2-stage model
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
pistar_by_chain_2stage <- function(n_chains, chains_list, Z, n, n_cat){

  colmeans_by_chain = lapply(chains_list, colMeans)

  pistar_results = matrix(NA, nrow = n_chains, ncol = n_cat)
  for(i in 1:n_chains){
    chain_means_i = colmeans_by_chain[[i]]
    pistar_matrix_i = pistar_compute_for_chains_2stage(chain_means_i, Z, n, n_cat)
    pistar_results[i, 1] = mean_pistarjj_compute(pistar_matrix_i, 1, n)
    pistar_results[i, 2] = mean_pistarjj_compute(pistar_matrix_i, 2, n)
  }
  return(pistar_results)
}
