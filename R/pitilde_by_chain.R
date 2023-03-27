#' Compute the Mean Conditional Probability of Second-Stage Correct Classification, by First-Stage and True Outcome Across all Subjects for each MCMC Chain
#'
#' @param n_chains An integer specifying the number of MCMC chains to compute over.
#' @param chains_list A numeric list containing the samples from \code{n_chains}
#'   MCMC chains.
#' @param V A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{V}.
#' @param n_cat The number of categorical values that the true outcome, \eqn{Y},
#'   the first-stage observed outcome, \eqn{Y*}, and the second-stage
#'   observed outcome, \eqn{\tilde{Y}}, can take.
#'
#' @return \code{pitilde_by_chain} returns a numeric matrix of the average
#'   conditional probability \eqn{P( \tilde{Y} = j | Y^* = j, Y = j, V)} across all subjects for
#'   each MCMC chain. Rows of the matrix correspond to MCMC chains, up to \code{n_chains}.
#'   The first column contains the conditional probability \eqn{P( \tilde{Y} = 1 | Y^* = 1, Y = 1, V)}.
#'   The second column contains the conditional probability \eqn{P( \tilde{Y} = 2 | Y^* = 2, Y = 2, V)}.
#'
#' @include pistar_compute_for_chains.R
#' @include mean_pistarjj_compute.R
#'
#' @importFrom stats rnorm
#'
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
