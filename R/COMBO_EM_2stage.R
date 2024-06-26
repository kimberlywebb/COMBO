#' EM-Algorithm Estimation of the Two-Stage Binary Outcome Misclassification Model
#'
#' Jointly estimate \eqn{\beta}, \eqn{\gamma^{(1)}}, \eqn{\gamma^{(2)}} parameters from the true outcome,
#' first-stage observation, and second-stage observation mechanisms, respectively,
#' in a two-stage binary outcome misclassification model.
#'
#' @param Ystar1 A numeric vector of indicator variables (1, 2) for the first-stage observed
#'   outcome \eqn{Y^{*(1)}}. There should be no \code{NA} terms. The reference category is 2.
#' @param Ystar2 A numeric vector of indicator variables (1, 2) for the second-stage
#'   observed outcome \eqn{Y^{*(2)}}. There should be no \code{NA} terms. The
#'   reference category is 2.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z1_matrix A numeric matrix of covariates in the first-stage observation mechanism.
#'   \code{z1_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z2_matrix A numeric matrix of covariates in the second-stage observation mechanism.
#'   \code{z2_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true outcome mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} plus 1.
#' @param gamma1_start A numeric vector or matrix of starting values for the \eqn{\gamma^{(1)}}
#'   parameters in the first-stage observation mechanism. In matrix form, the \code{gamma1_start} matrix rows
#'   correspond to parameters for the \eqn{Y^{*(1)} = 1}
#'   first-stage observed outcome, with the dimensions of \code{z1_matrix} plus 1, and the
#'   parameter matrix columns correspond to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}. A numeric vector for \code{gamma1_start} is
#'   obtained by concatenating the matrix, i.e. \code{gamma1_start <- c(gamma1_matrix)}.
#' @param gamma2_start A numeric array of starting values for the \eqn{\gamma^{(2)}} parameters
#'   in the second-stage observation mechanism. The first dimension (matrix rows)
#'   of \code{gamma2_start} correspond to parameters for the \eqn{Y^{*(2)} = 1}
#'   second-stage observed outcome, with the dimensions of the \code{z2_matrix}
#'   plus 1. The second dimension (matrix columns) correspond to the first-stage
#'   observed outcome categories \eqn{Y^{*(1)} \in \{1, 2\}}. The third dimension of
#'   \code{gamma2_start} corresponds to to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations An integer specifying the maximum number of
#'   iterations of the EM algorithm. The default is \code{1500}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMBO_EM_2stage} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}, followed by the standard
#'   error estimates, \code{SE}. The final column, \code{Convergence}, reports
#'   whether or not the algorithm converged for a given parameter estimate.
#'
#'   Estimates are provided for the two-stage binary misclassification model.
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
#' @include em_function_2stage.R
#' @include loglik_2stage.R
#' @include COMBO_data_2stage.R
#'
#' @importFrom stats rnorm rgamma rmultinom glm
#' @importFrom turboEM turboem
#' @importFrom Matrix nearPD
#'
#' @examples \donttest{
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z1_shape <- 1
#' z2_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma1 <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#' true_gamma2 <- array(c(1.5, 1, .5, .5, -.5, 0, -1, -1), dim = c(2, 2, 2))
#'
#' my_data <- COMBO_data_2stage(sample_size = n,
#'                              x_mu = x_mu, x_sigma = x_sigma,
#'                              z1_shape = z1_shape, z2_shape = z2_shape,
#'                              beta = true_beta, gamma1 = true_gamma1, gamma2 = true_gamma2)
#' table(my_data[["obs_Ystar2"]], my_data[["obs_Ystar1"]], my_data[["true_Y"]])
#'
#' beta_start <- rnorm(length(c(true_beta)))
#' gamma1_start <- rnorm(length(c(true_gamma1)))
#' gamma2_start <- rnorm(length(c(true_gamma2)))
#'
#' EM_results <- COMBO_EM_2stage(Ystar1 = my_data[["obs_Ystar1"]],
#'                               Ystar2 = my_data[["obs_Ystar2"]],
#'                               x_matrix = my_data[["x"]],
#'                               z1_matrix = my_data[["z1"]],
#'                               z2_matrix = my_data[["z2"]],
#'                               beta_start = beta_start,
#'                               gamma1_start = gamma1_start,
#'                               gamma2_start = gamma2_start)
#'
#' EM_results}
COMBO_EM_2stage <- function(Ystar1, Ystar2,
                            x_matrix, z1_matrix, z2_matrix,
                            beta_start, gamma1_start, gamma2_start,
                            tolerance = 1e-7, max_em_iterations = 1500,
                            em_method = "squarem"){

  Ystar <- Ystar1
  Ytilde <- Ystar2
  z_matrix <- z1_matrix
  v_matrix <- z2_matrix
  gamma_start <- gamma1_start
  delta_start <- gamma2_start

  if (is.data.frame(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.numeric(z_matrix))
    stop("'z1_matrix' should be a numeric matrix.")

  if (is.vector(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.matrix(z_matrix))
    stop("'z1_matrix' should be a matrix or data.frame.")

  if (is.vector(v_matrix))
    v_matrix <- as.matrix(v_matrix)
  if (!is.matrix(v_matrix))
    stop("'z2_matrix' should be a matrix or data.frame.")

  if (!is.null(x_matrix)) {
    if (is.data.frame(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.numeric(x_matrix))
      stop("'x_matrix' must be numeric.")
    if (is.vector(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.matrix(x_matrix))
      stop("'x_matrix' must be a data.frame or matrix.")
  }

  if (!is.numeric(Ystar) || !is.vector(Ystar))
    stop("'Ystar1' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ystar))) != 0)
    stop("'Ystar1' must be coded 1/2, where the reference category is 2.")

  if (!is.numeric(Ytilde) || !is.vector(Ytilde))
    stop("'Ystar2' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ytilde))) != 0)
    stop("'Ystar2' must be coded 1/2, where the reference category is 2.")

  n_cat = 2
  sample_size = length(Ystar)
  sample_size_tilde <- length(Ytilde)

  if (sample_size_tilde != sample_size)
    stop("The lengths of 'Ystar' and 'Ytilde' must be the same.")
  if (nrow(z_matrix) != sample_size)
    stop("The number of rows of 'z_matrix' must match the length of 'Ystar'.")
  if (nrow(v_matrix) != sample_size_tilde)
    stop("The number of rows of 'v_matrix' must match the length of 'Ytilde'.")
  if (!is.null(x_matrix) && nrow(x_matrix) != sample_size)
    stop("The number of rows of 'x_matrix' must match the length of 'Ystar'.")

  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  V = matrix(c(rep(1, sample_size), c(v_matrix)),
             byrow = FALSE, nrow = sample_size)

  obs_Ystar_reps = matrix(rep(Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_Ystar_matrix = 1 * (obs_Ystar_reps == category_matrix)

  obs_Ytilde_reps <- matrix(rep(Ytilde, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ytilde_matrix <- 1 * (obs_Ytilde_reps == category_matrix)

  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)

  results = turboEM::turboem(par = c(c(beta_start), c(gamma_start), c(delta_start)),
                             fixptfn = em_function_2stage, objfn = loglik_2stage,
                             method = c(em_method),
                             obs_Ystar_matrix = obs_Ystar_matrix,
                             obs_Ytilde_matrix = obs_Ytilde_matrix,
                             X = X, Z = Z, V = V,
                             sample_size = sample_size, n_cat = n_cat,
                             control.run = control_settings)

  # Naive model
  Ystar_01 <- ifelse(Ystar == 1, 1, 0)
  Ytilde_01 <- ifelse(Ytilde == 1, 1, 0)

  outcome1_model <- glm(Ystar_01 ~ X[,-1], family = "binomial")

  outcome2_ystar1_model <- glm(Ytilde_01[which(Ystar == 1)] ~ V[which(Ystar == 1), -1],
                               family = "binomial")
  outcome2_ystar2_model <- glm(Ytilde_01[which(Ystar == 2)] ~ V[which(Ystar == 2), -1],
                               family = "binomial")

  naive_param <- c(coef(summary(outcome1_model))[,2],
                   coef(summary(outcome2_ystar1_model))[,2],
                   coef(summary(outcome2_ystar2_model))[,2])

  naive_se <- c(coef(summary(outcome1_model))[,3],
                coef(summary(outcome2_ystar1_model))[,3],
                coef(summary(outcome2_ystar2_model))[,3])

  naive_convergence <- c(rep(outcome1_model$converged, nrow(coef(summary(outcome1_model)))),
                         rep(outcome2_ystar1_model$converged, nrow(coef(summary(outcome2_ystar1_model)))),
                         rep(outcome2_ystar2_model$converged, nrow(coef(summary(outcome2_ystar2_model)))))


  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)

  results_i_delta <- array(turboEM::pars(results)[((ncol(X) + (n_cat * ncol(Z))) + 1):length(turboEM::pars(results))],
                           dim = c(ncol(V), n_cat, n_cat))
  results_i_pitilde <- pitilde_compute(results_i_delta, V, sample_size, n_cat)

  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])

  flip_pistar11 <- 1 - pistar_22
  flip_pistar22 <- 1 - pistar_11

  J <- pistar_11 + pistar_22 - 1
  J_flip <- flip_pistar11 + flip_pistar22 - 1

  pitilde_111 <- mean(results_i_pitilde[1:sample_size, 1, 1])
  pitilde_222 <- mean(results_i_pitilde[(sample_size + 1):(2*sample_size), 2, 2])

  estimates_i <- if ((J_flip <= J) |
                     (is.na(pistar_11) & is.na(pistar_22))) {
    # If turboem cannot estimate the parameters they will be NA.
    turboEM::pars(results)
  } else {
    gamma_index = (ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)

    delta_index = ((ncol(X) + (n_cat * ncol(Z))) + 1):length(turboEM::pars(results))
    n_delta_param = length(delta_index) / n_cat
    delta_flip_index = (ncol(X) + (n_cat * ncol(Z))) + c((n_delta_param + 1):length(delta_index), 1:n_delta_param)

    c(-1*turboEM::pars(results)[1:ncol(X)], turboEM::pars(results)[c(gamma_flip_index, delta_flip_index)])

  }

  #sigma_EM = tryCatch(solve(turboEM::hessian(results)[[1]]), silent = TRUE,
  #                    error = function(e) NA)
  #SE_EM = tryCatch(sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
  #                                  row = length(c(c(beta_start), c(gamma_start))),
  #                                  byrow = FALSE))),
  #                 silent = TRUE,
  #                 error = function(e) rep(NA, ncol(X) + (n_cat * ncol(Z))))

  # Do label switching for the SE estimates too.
  sigma_EM = tryCatch(solve(turboEM::hessian(results)[[1]]),
                      silent = TRUE,
                      error = function(e) matrix(NA,
                                                 nrow = length(c(c(beta_start), c(gamma_start), c(delta_start))),
                                                 ncol = length(c(c(beta_start), c(gamma_start), c(delta_start)))))
  #SE_EM = sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
   #                        nrow = length(c(c(beta_start), c(gamma_start), c(delta_start))),
    #                       byrow = FALSE)))

  SE_EM <- if ((J_flip <= J) |
               (is.na(pistar_11) & is.na(pistar_22))) {

    tryCatch(sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
                              nrow = length(c(c(beta_start), c(gamma_start), c(delta_start))),
                              byrow = FALSE))),
             silent = TRUE,
             error = function(e) rep(NA, length = nrow(sigma_EM)))

  } else {
    gamma_index = (ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)

    delta_index = ((ncol(X) + (n_cat * ncol(Z))) + 1):length(turboEM::pars(results))
    n_delta_param = length(delta_index) / n_cat
    delta_flip_index = (ncol(X) + (n_cat * ncol(Z))) + c((n_delta_param + 1):length(delta_index), 1:n_delta_param)

    tryCatch(sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
                              nrow = length(c(c(beta_start), c(gamma_start), c(delta_start))),
                              byrow = FALSE)))[c(1:ncol(X), gamma_flip_index, delta_flip_index)],
             silent = TRUE,
             error = function(e) rep(NA, length = nrow(sigma_EM)))

  }

  beta_param_names <- paste0(rep("beta_", ncol(X)), 1:ncol(X))
  gamma_param_names <- paste0(rep("gamma1_", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  delta_param_names <- paste0(rep("gamma2_", length(c(delta_start))),
                              rep(1:ncol(V), n_cat * n_cat),
                              rep(1, length(c(delta_start))),
                              rep(c(rep(1, ncol(V)), rep(2, ncol(V))), n_cat),
                              c(rep(1, ncol(V) * n_cat), rep(2, ncol(V) * n_cat)))
  naive_param_names <- paste0("naive_", c(beta_param_names,
                                          paste0(rep("gamma2_", (n_cat * ncol(V))),
                                                 rep(1:ncol(V), n_cat),
                                                 rep(1:n_cat, each = ncol(V)))))

  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names,
                                        delta_param_names,
                                        naive_param_names),
                          Estimates = c(c(estimates_i),
                                        c(naive_param)),
                          SE = c(SE_EM, naive_se),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names,
                                                       delta_param_names))),
                                          naive_convergence))

  return(estimates)
}
