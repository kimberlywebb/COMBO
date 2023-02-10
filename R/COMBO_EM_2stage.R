#' EM-Algorithm Estimation of the Two-Stage Binary Outcome Misclassification Model
#'
#' Jointly estimate \eqn{\beta}, \eqn{\gamma}, \eqn{\delta} parameters from the true outcome,
#' first-stage observation, and second-stage observation mechanisms, respectively,
#' in a two-stage binary outcome misclassification model.
#'
#' @param Ystar A numeric vector of indicator variables (1, 2) for the first-stage observed
#'   outcome \code{Y*}. There should be no \code{NA} terms. The reference category is 2.
#' @param Ytilde A numeric vector of indicator variables (1, 2) for the second-stage
#'   observed outcome \eqn{\tilde{Y}}. There should be no \code{NA} terms. The
#'   reference category is 2.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the first-stage observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param v_matrix A numeric matrix of covariates in the second-stage observation mechanism.
#'   \code{v_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true outcome mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} plus 1.
#' @param gamma_start A numeric vector or matrix of starting values for the \eqn{\gamma}
#'   parameters in the first-stage observation mechanism. In matrix form, the \code{gamma_start} matrix rows
#'   correspond to parameters for the \code{Y* = 1}
#'   first-stage observed outcome, with the dimensions of \code{z_matrix} plus 1, and the
#'   gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
#' @param delta_start A numeric array of starting values for the \eqn{\delta} parameters
#'   in the second-stage observation mechanism. The first dimension (matrix rows)
#'   of \code{delta_start} correspond to parameters for the \eqn{\tilde{Y} = 1}
#'   second-stage observed outcome, with the dimensions of the \code{v_matrix}
#'   plus 1. The second dimension (matrix columns) correspond to the first-stage
#'   observed outcome categories \eqn{Y^* \in \{1, 2\}}. The third dimension of
#'   \code{delta_start} corresponds to to the true outcome categories
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
#' @examples
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
#' table(my_data[["obs_Ytilde"]], my_data[["obs_Ystar"]], my_data[["true_Y"]])
#'
#' beta_start <- rnorm(length(c(true_beta)))
#' gamma_start <- rnorm(length(c(true_gamma)))
#' delta_start <- rnorm(length(c(true_delta)))
#'
#' EM_results <- COMBO_EM_2stage(Ystar = my_data[["obs_Ystar"]],
#'                               Ytilde = my_data[["obs_Ytilde"]],
#'                               x_matrix = my_data[["x"]],
#'                               z_matrix = my_data[["z"]],
#'                               v_matrix = my_data[["v"]],
#'                               beta_start = beta_start,
#'                               gamma_start = gamma_start,
#'                               delta_start = delta_start)
#'
#' EM_results
COMBO_EM_2stage <- function(Ystar, Ytilde,
                            x_matrix, z_matrix, v_matrix,
                            beta_start, gamma_start, delta_start,
                            tolerance = 1e-7, max_em_iterations = 1500,
                            em_method = "squarem"){

  if (is.data.frame(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.numeric(z_matrix))
    stop("'z_matrix' should be a numeric matrix.")

  if (is.vector(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.matrix(z_matrix))
    stop("'z_matrix' should be a matrix or data.frame.")

  if (is.vector(v_matrix))
    v_matrix <- as.matrix(v_matrix)
  if (!is.matrix(v_matrix))
    stop("'v_matrix' should be a matrix or data.frame.")

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
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ystar))) != 0)
    stop("'Ystar' must be coded 1/2, where the reference category is 2.")

  if (!is.numeric(Ytilde) || !is.vector(Ytilde))
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ytilde))) != 0)
    stop("'Ystar' must be coded 1/2, where the reference category is 2.")

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
  naive_start_beta <- glm(Ystar_01 ~ X[,-1], family = "binomial")
  naive_start_delta1 <- glm(Ytilde_01[which(Ystar == 1)] ~ V[which(Ystar == 1), -1],
                            family = "binomial")
  naive_start_delta2 <- glm(Ytilde_01[which(Ystar == 2)] ~ V[which(Ystar == 2), -1],
                            family = "binomial")
  naive_results <- optim(par = c(unname(coef(naive_start_beta)),
                                 unname(coef(naive_start_delta1)),
                                 unname(coef(naive_start_delta2))),
                         fn = naive_loglik_2stage,
                         X = X, V = V,
                         obs_Ystar_matrix = obs_Ystar_matrix,
                         obs_Ytilde_matrix = obs_Ytilde_matrix,
                         sample_size = sample_size,
                         n_cat = n_cat,
                         hessian = TRUE,
                         control = list(maxit = max_em_iterations))
  naive_se <- tryCatch(sqrt(diag(solve(naive_results$hessian))),
                       silent = TRUE,
                       error = function(e) rep(NA, length(c(unname(coef(naive_start_beta)),
                                                            unname(coef(naive_start_delta1)),
                                                            unname(coef(naive_start_delta2))))))
  naive_convergence <- ifelse(naive_results$convergence == 1, "maxit reached",
                              ifelse(naive_results$convergence == 10,
                                     "degenerency in nelder-mead simplex",
                                     ifelse(naive_results$convergence == 51,
                                            "warning from L-BFGS-B",
                                            ifelse(naive_results$convergence == 52,
                                                   "error from L-BFGS-B",
                                                   ifelse(naive_results$convergence == 0,
                                                          TRUE, NA)))))

  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)

  results_i_delta <- array(turboEM::pars(results)[((ncol(X) + (n_cat * ncol(Z))) + 1):length(turboEM::pars(results))],
                           dim = c(ncol(V), n_cat, n_cat))
  results_i_pitilde <- pitilde_compute(results_i_delta, V, sample_size, n_cat)

  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])

  pitilde_111 <- mean(results_i_pitilde[1:sample_size, 1, 1])
  pitilde_222 <- mean(results_i_pitilde[(sample_size + 1):(2*sample_size), 2, 2])

  estimates_i <- if ((pistar_11 > .50 | pistar_22 > .50 | pitilde_111 > .50 | pitilde_222 > .50) |
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

  SE_EM <- if ((pistar_11 > .50 | pistar_22 > .50 | pitilde_111 > .50 | pitilde_222 > .50) |
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

  beta_param_names <- paste0(rep("beta", ncol(X)), 1:ncol(X))
  gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  delta_param_names <- paste0(rep("delta", length(c(delta_start))),
                              rep(1:ncol(V), n_cat * n_cat),
                              rep(1, length(c(delta_start))),
                              rep(c(rep(1, ncol(V)), rep(2, ncol(V))), n_cat),
                              c(rep(1, ncol(V) * n_cat), rep(2, ncol(V) * n_cat)))
  naive_param_names <- paste0("naive_", c(beta_param_names,
                                          paste0(rep("delta", (n_cat * ncol(V))),
                                                 rep(1:ncol(V), n_cat),
                                                 rep(1:n_cat, each = ncol(V)))))

  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names,
                                        delta_param_names,
                                        naive_param_names),
                          Estimates = c(c(estimates_i),
                                        c(naive_results$par)),
                          SE = c(SE_EM, naive_se),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names,
                                                       delta_param_names))),
                                          rep(naive_convergence,
                                              length(naive_param_names))))

  return(estimates)
}
