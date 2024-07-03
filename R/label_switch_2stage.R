#' Fix Label Switching in MCMC Results from a Binary Outcome Misclassification Model
#'
#' @param chain_matrix A numeric matrix containing the posterior samples for all
#'   parameters in a given MCMC chain. \code{chain_matrix} must be a named
#'   object (i.e. each parameter must be named as \code{beta[j, p]}, \code{gamma[k,j,p]},
#'   or \code{delta[l,k,j,p]}).
#' @param dim_x An integer specifying the number of columns of the design matrix of the true outcome mechanism, \code{X}.
#' @param dim_z An integer specifying the number of columns of the design matrix of the first-stage observation mechanism, \code{Z}.
#' @param dim_v An integer specifying the number of columns of the design matrix of the second-stage observation mechanism, \code{V}.
#' @param n_cat An integer specifying the number of categorical values that the true outcome, \eqn{Y},
#'   the first-stage observed outcome, \eqn{Y^*}, and the second-stage observed
#'   outcome \eqn{\tilde{Y}} can take.
#'
#' @return \code{label_switch_2stage} returns a named matrix of MCMC posterior samples for
#'   all parameters after performing label switching according the following pattern:
#'   all \eqn{\beta} terms are multiplied by -1, all \eqn{\gamma} and \eqn{\delta} terms are "swapped"
#'   with the opposite \code{j} index.
#'
#' @importFrom stats rnorm rmultinom
#'
label_switch_2stage <- function(chain_matrix, dim_x, dim_z, dim_v, n_cat){

  beta_names <- paste0("beta[1,", 1:dim_x, "]")
  gamma_names <- paste0("gamma1[1,", rep(1:n_cat, dim_z), ",", rep(1:dim_z, each = n_cat), "]")
  delta_names <- paste0("gamma2[1,",
                        rep(1:n_cat, dim_v*dim_v), ",",
                        rep(rep(1:n_cat, each = dim_v), dim_v), ",",
                        rep(1:dim_v, each = n_cat * n_cat), "]")

  beta_cols <- which(colnames(chain_matrix) %in% beta_names)
  gamma_cols <- which(colnames(chain_matrix) %in% gamma_names)
  delta_cols <- which(colnames(chain_matrix) %in% delta_names)

  n_flip_gammas <- length(gamma_cols) / 2
  gamma_cols_1 <- gamma_cols[c(TRUE, FALSE)]
  gamma_cols_2 <- gamma_cols[c(FALSE, TRUE)]

  n_flip_deltas <- length(delta_cols) / 2
  delta_cols_1 <- delta_cols[c(TRUE, TRUE, FALSE, FALSE)]
  delta_cols_2 <- delta_cols[c(FALSE, FALSE, TRUE, TRUE)]

  return_chain_matrix <- chain_matrix

  return_chain_matrix[,beta_cols] <- chain_matrix[,beta_cols] * -1

  for(i in 1:n_flip_gammas){
    return_chain_matrix[,gamma_cols_1[i]] <- chain_matrix[,gamma_cols_2[i]]
    return_chain_matrix[,gamma_cols_2[i]] <- chain_matrix[,gamma_cols_1[i]]
  }

  for(i in 1:n_flip_deltas){
    return_chain_matrix[,delta_cols_1[i]] <- chain_matrix[,delta_cols_2[i]]
    return_chain_matrix[,delta_cols_2[i]] <- chain_matrix[,delta_cols_1[i]]
  }

  return(return_chain_matrix)
}
