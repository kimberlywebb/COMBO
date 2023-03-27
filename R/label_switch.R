#' Fix Label Switching in MCMC Results from a Binary Outcome Misclassification Model
#'
#' @param chain_matrix A numeric matrix containing the posterior samples for all
#'   parameters in a given MCMC chain. \code{chain_matrix} must be a named
#'   object (i.e. each parameter must be named as \code{beta[j, p]} or \code{gamma[k,j,p]}).
#' @param dim_x An integer specifying the number of columns of the design matrix of the true outcome mechanism, \code{X}.
#' @param dim_z An integer specifying the number of columns of the design matrix of the observation mechanism, \code{Z}.
#' @param n_cat An integer specifying the number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*} can take.
#'
#' @return \code{label_switch} returns a named matrix of MCMC posterior samples for
#'   all parameters after performing label switching according the following pattern:
#'   all \eqn{\beta} terms are multiplied by -1, all \eqn{\gamma} terms are "swapped"
#'   with the opposite \code{j} index.
#'
#' @importFrom stats rnorm rmultinom
#'
label_switch <- function(chain_matrix, dim_x, dim_z, n_cat){

  beta_names <- paste0("beta[1,", 1:dim_x, "]")
  gamma_names <- paste0("gamma[1,", rep(1:n_cat, dim_z), ",", rep(1:dim_z, each = n_cat), "]")

  beta_cols <- which(colnames(chain_matrix) %in% beta_names)
  gamma_cols <- which(colnames(chain_matrix) %in% gamma_names)

  n_flip_gammas <- length(gamma_cols) / 2
  gamma_cols_1 <- gamma_cols[c(TRUE, FALSE)]
  gamma_cols_2 <- gamma_cols[c(FALSE, TRUE)]

  return_chain_matrix <- chain_matrix

  return_chain_matrix[,beta_cols] <- chain_matrix[,beta_cols] * -1

  for(i in 1:n_flip_gammas){
    return_chain_matrix[,gamma_cols_1[i]] <- chain_matrix[,gamma_cols_2[i]]
    return_chain_matrix[,gamma_cols_2[i]] <- chain_matrix[,gamma_cols_1[i]]
  }

  return(return_chain_matrix)
}
