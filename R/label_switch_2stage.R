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
#' @examples \dontrun{
#' set.seed(123)
#' n <- 1000
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
#'
#' obs_Y = my_data[["obs_Y"]]
#' X = my_data[["x_design_matrix"]]
#' Z = my_data[["z_design_matrix"]]
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
#'
#' label_switch_chain1 <- label_switch(example_chain1, ncol(X), ncol(Z), n_cat = 2)
#'
#' head(example_chain1)
#' head(label_switch_chain1)
#' }
label_switch_2stage <- function(chain_matrix, dim_x, dim_z, dim_v, n_cat){

  beta_names <- paste0("beta[1,", 1:dim_x, "]")
  gamma_names <- paste0("gamma[1,", rep(1:n_cat, dim_z), ",", rep(1:dim_z, each = n_cat), "]")
  delta_names <- paste0("delta[1,",
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
