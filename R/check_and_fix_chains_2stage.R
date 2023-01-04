#' Check Assumption and Fix Label Switching if Assumption is Broken for a List of MCMC Samples
#'
#' @param n_chains An integer specifying the number of MCMC chains to compute over.
#' @param chains_list A numeric list containing the samples from \code{n_chains}
#'   MCMC chains.
#' @param pistarjj_matrix A numeric matrix of the average
#'   conditional probability \eqn{P(Y^* = j | Y = j, Z)} across all subjects for
#'   each MCMC chain, obtained from the \code{pistar_by_chain} function.
#' @param dim_x The number of columns of the design matrix of the true outcome mechanism, \code{X}.
#' @param dim_z The number of columns of the design matrix of the observation mechanism, \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{check_and_fix_chains} returns a numeric list of the samples from
#'   \code{n_chains} MCMC chains which have been corrected for label switching if
#'   the following assumption is not met: \eqn{P(Y^* = j | Y = j, Z) > 0.50 \forall j}.
#'
#' @include label_switch.R
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 100
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#'
#' beta <- matrix(c(1, 2), ncol = 1)
#' gamma <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = FALSE)
#'
#' my_data <- COMBO_data(sample_size = n,
#'                       x_mu = x_mu, x_sigma = x_sigma,
#'                       z_shape = z_shape,
#'                       beta = beta, gamma = gamma)
#' X <- my_data[["x_design_matrix"]]
#' Z <- my_data[["z_design_matrix"]]
#'
#' example_chain1 <- matrix(c(rnorm(n, mean = beta[1,1]),
#'                            rnorm(n, mean = beta[2,1]),
#'                            rnorm(n, mean = gamma[1,1]),
#'                            rnorm(n, mean = gamma[1,2]),
#'                            rnorm(n, mean = gamma[2,1]),
#'                            rnorm(n, mean = gamma[2,2])),
#'                            nrow = n, byrow = FALSE)
#' colnames(example_chain1) <- c("beta[1,1]", "beta[1,2]",
#'                               "gamma[1,1,1]", "gamma[1,2,1]",
#'                               "gamma[1,1,2]", "gamma[1,2,2]")
#' example_chain2 <- matrix(c(rnorm(n, mean = beta[1,1]),
#'                            rnorm(n, mean = beta[2,1]),
#'                            rnorm(n, mean = gamma[1,1]),
#'                            rnorm(n, mean = gamma[1,2]),
#'                            rnorm(n, mean = gamma[2,1]),
#'                            rnorm(n, mean = gamma[2,2])),
#'                            nrow = n, byrow = FALSE)
#' colnames(example_chain2) <- c("beta[1,1]", "beta[1,2]",
#'                               "gamma[1,1,1]", "gamma[1,2,1]",
#'                               "gamma[1,1,2]", "gamma[1,2,2]")
#' chains_list <- list(example_chain1, example_chain2)
#' pistar_by_chain_matrix <- pistar_by_chain(n_chains = 2,
#'                                           chains_list = chains_list,
#'                                           Z = Z, n = n, n_cat = 2)
#' fixed_chains <- check_and_fix_chains(n_chains = 2, chains_list = chains_list,
#'                                      pistarjj_matrix = pistar_by_chain_matrix,
#'                                      dim_x = ncol(X), dim_z = ncol(Z),
#'                                      n_cat = 2)
#'
#' pistar_by_chain_matrix
#' chains_list[[1]][1:5,]
#' fixed_chains[[1]][1:5,]
#' chains_list[[2]][1:5,]
#' fixed_chains[[2]][1:5,]
#' }
check_and_fix_chains_2stage <- function(n_chains, chains_list,
                                        pistarjj_matrix, pitildejjj_matrix,
                                        dim_x, dim_z, dim_v, n_cat){
  fixed_output_list <- list()
  for(i in 1:n_chains){
    pistar11 = pistarjj_matrix[i, 1]
    pistar22 = pistarjj_matrix[i, 2]
    pitilde111 = pitildejjj_matrix[i, 1]
    pitilde222 = pitildejjj_matrix[i, 2]
    output = if(pistar11 > .50 & pistar22 > .50 &
                pitilde111 > .50 & pitilde222 > .50){
      chains_list[[i]]
    } else {label_switch_2stage(chains_list[[i]],
                                dim_x, dim_z, dim_v, n_cat)}
    fixed_output_list[[i]] = output
  }
  return(fixed_output_list)
}
