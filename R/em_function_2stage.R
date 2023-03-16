#' EM-Algorithm Function for Estimation of the Two-Stage Misclassification Model
#'
#' @param param_current A numeric vector of regression parameters, in the order
#'   \eqn{\beta, \gamma, \delta}. The \eqn{\gamma} vector is obtained from the matrix form.
#'   In matrix form, the gamma parameter matrix rows
#'   correspond to parameters for the \code{Y* = 1}
#'   observed outcome, with the dimensions of \code{Z}.
#'   In matrix form, the gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{j = 1, \dots,} \code{n_cat}. The numeric vector \eqn{\gamma} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_v <- c(gamma_matrix)}.
#'   The \eqn{\delta} vector is obtained from the array form. In array form,
#'   the first dimension (matrix rows) of \code{delta}
#'   corresponds to parameters for the \eqn{\tilde{Y} = 1}
#'   second-stage observed outcome, with the dimensions of the \code{V}
#'   The second dimension (matrix columns) correspond to the first-stage
#'   observed outcome categories \eqn{Y^* \in \{1, 2\}}. The third dimension of
#'   \code{delta_start} corresponds to to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}. The numeric vector \eqn{\delta} is obtained by
#'   concatenating the delta array, i.e. \code{delta_vector <- c(delta_array)}.
#' @param obs_Ystar_matrix A numeric matrix of indicator variables (0, 1) for the first-stage observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param obs_Ytilde_matrix A numeric matrix of indicator variables (0, 1) for the second-stage observed
#'   outcome \eqn{\tilde{Y}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param X A numeric design matrix for the true outcome mechanism.
#' @param Z A numeric design matrix for the first-stage observation mechanism.
#' @param V A numeric design matrix for the second-stage observation mechanism.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrices, \code{X}, \code{Z}, and \code{V}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcomes, \code{Y*} and \eqn{\tilde{Y}}, can take.
#'
#' @return \code{em_function_2stage} returns a numeric vector of updated parameter
#'   estimates from one iteration of the EM-algorithm.
#'
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include pitilde_compute.R
#' @include w_j_2stage.R
#' @include q_beta_f.R
#' @include q_gamma_f.R
#' @include q_delta_f.R
#'
#' @importFrom stats rnorm rgamma rmultinom coefficients binomial glm
#'
#' @examples \dontrun{
#'
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#' v_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#' true_delta <- array(c(1.5, 1, .5, .5, -.5, 0, -1, -1), dim = c(2, 2, 2))
#'
#' my_data <- COMBO_data_2stage(sample_size = n,
#'                              x_mu = x_mu, x_sigma = x_sigma,
#'                              z_shape = z_shape, v_shape = v_shape,
#'                              beta = true_beta, gamma = true_gamma, delta = true_delta)
#'
#' obs_Ystar_matrix = my_data[["obs_Ystar_matrix"]]
#' obs_Ytilde_matrix = my_data[["obs_Ytilde_matrix"]]
#' X = my_data[["x_design_matrix"]]
#' Z = my_data[["z_design_matrix"]]
#' V = my_data[["v_design_matrix"]]
#'
#' starting_values <- rnorm(14)
#'
#' new_parameters <- em_function_2stage(starting_values,
#'                                      obs_Ystar_matrix = obs_Ystar_matrix,
#'                                      obs_Ytilde_matrix = obs_Ytilde_matrix,
#'                                      X = X, Z = Z, V = V,
#'                                      sample_size = n, n_cat = 2)
#'
#' new_parameters
#' }
em_function_2stage <- function(param_current,
                               obs_Ystar_matrix, obs_Ytilde_matrix,
                               X, Z, V,
                               sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  gamma_current = matrix(c(param_current[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))]),
                         ncol = n_cat, byrow = FALSE)
  delta_current = array(c(param_current[(ncol(X) + (n_cat * ncol(Z)) + 1):length(param_current)]),
                        dim = c(ncol(V), 2, 2))

  probabilities = matrix(pi_compute(beta_current, X, sample_size, n_cat),
                         ncol = n_cat, byrow = FALSE)
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  conditional_probabilities2 = pitilde_compute(delta_current, V, sample_size, n_cat)

  weights = w_j_2stage(ystar_matrix = obs_Ystar_matrix,
                       ytilde_matrix = obs_Ytilde_matrix,
                       pitilde_array = conditional_probabilities2,
                       pistar_matrix = conditional_probabilities,
                       pi_matrix = probabilities,
                       sample_size = sample_size, n_cat = n_cat)

  Ystar01 = obs_Ystar_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                                             weights = weights[,1],
                                             family = "binomial"(link = "logit")) )
  gamma1_new <- unname(stats::coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                                             weights = weights[,2],
                                             family = "binomial"(link = "logit")) )
  gamma2_new <- unname(stats::coefficients(fit.gamma2))

  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(X),
                                           family = stats::binomial()) )
  beta_new <- unname(stats::coefficients(fit.beta))


  outcome_j1_k1 <- ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,1] == 1,
                          1,
                          ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,1] == 0,
                                 0, NA))

  fit.delta11 <- suppressWarnings( stats::glm(outcome_j1_k1 ~ . + 0, as.data.frame(V),
                                              weights = weights[,1],
                                              family = "binomial"(link = "logit")) )
  delta11_new <- unname(coefficients(fit.delta11))


  outcome_j1_k2 <- ifelse(obs_Ystar_matrix[,1] == 0 & obs_Ytilde_matrix[,1] == 1,
                          1,
                          ifelse(obs_Ystar_matrix[,1] == 0 & obs_Ytilde_matrix[,1] == 0,
                                 0, NA))

  fit.delta21 <- suppressWarnings( stats::glm(outcome_j1_k2 ~ . + 0, as.data.frame(V),
                                              weights = weights[,1],
                                              family = "binomial"(link = "logit")) )
  delta21_new <- unname(coefficients(fit.delta21))


  outcome_j2_k1 <- ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,1] == 1,
                          1,
                          ifelse(obs_Ystar_matrix[,1] == 1 & obs_Ytilde_matrix[,1] == 0,
                                 0, NA))

  fit.delta12 <- suppressWarnings( stats::glm(outcome_j2_k1 ~ . + 0, as.data.frame(V),
                                              weights = weights[,2],
                                              family = "binomial"(link = "logit")) )
  delta12_new <- unname(coefficients(fit.delta12))


  outcome_j2_k2 <- ifelse(obs_Ystar_matrix[,1] == 0 & obs_Ytilde_matrix[,1] == 1,
                          1,
                          ifelse(obs_Ystar_matrix[,1] == 0 & obs_Ytilde_matrix[,1] == 0,
                                 0, NA))

  fit.delta22 <- suppressWarnings( stats::glm(outcome_j2_k2 ~ . + 0, as.data.frame(V),
                                              weights = weights[,2],
                                              family = "binomial"(link = "logit")) )
  delta22_new <- unname(coefficients(fit.delta22))

  delta_new <- c(delta11_new, delta21_new, delta12_new, delta22_new)

  param_new = c(beta_new, gamma1_new, gamma2_new, delta_new)
  param_current = param_new
  return(param_new)

}
