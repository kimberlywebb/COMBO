#' EM-Algorithm Estimation of the Binary Outcome Misclassification Model
#'
#' Jointly estimate \eqn{\beta} and \eqn{\gamma} parameters from the true outcome
#' and observation mechanisms, respectively, in a binary outcome misclassification
#' model.
#'
#' @param Ystar A numeric vector of indicator variables (1, 2) for the observed
#'   outcome \code{Y*}. There should be no \code{NA} terms. The reference category is 2.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true outcome mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} plus 1.
#' @param gamma_start A numeric vector or matrix of starting values for the \eqn{\gamma}
#'   parameters in the observation mechanism. In matrix form, the \code{gamma_start} matrix rows
#'   correspond to parameters for the \code{Y* = 1}
#'   observed outcome, with the dimensions of \code{z_matrix} plus 1, and the
#'   gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{Y \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations An integer specifying the maximum number of
#'   iterations of the EM algorithm. The default is \code{1500}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMBO_EM} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}, followed by the standard
#'   error estimates, \code{SE}. The final column, \code{Convergence}, reports
#'   whether or not the algorithm converged for a given parameter estimate.
#'
#'   Estimates are provided for the binary misclassification model, as well as two
#'   additional cases. The "SAMBA" parameter estimates are from the R Package,
#'   SAMBA, which uses the EM algorithm to estimate a binary outcome misclassification
#'   model that assumes there is perfect specificity. The "PSens" parameter estimates
#'   are estimated using the EM algorithm for the binary outcome misclassification
#'   model that assumes there is perfect sensitivitiy. The "Naive" parameter
#'   estimates are from a simple logistic regression \code{Y* ~ X}.
#'
#' @references Beesley, L. and Mukherjee, B. (2020).
#'   Statistical inference for association studies using electronic health records:
#'   Handling both selection bias and outcome misclassification.
#'   Biometrics, 78, 214-226.
#'
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_j.R
#' @include q_beta_f.R
#' @include q_gamma_f.R
#' @include em_function.R
#' @include loglik.R
#' @include perfect_sensitivity_EM.R
#' @include COMBO_data.R
#'
#' @importFrom stats rnorm rgamma rmultinom glm
#' @importFrom turboEM turboem
#' @importFrom SAMBA obsloglikEM
#' @importFrom Matrix nearPD
#'
#' @examples
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#'
#' true_beta <- matrix(c(1, -2), ncol = 1)
#' true_gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#'
#' x_matrix = matrix(rnorm(n, x_mu, x_sigma), ncol = 1)
#' X = matrix(c(rep(1, n), x_matrix[,1]), ncol = 2, byrow = FALSE)
#' z_matrix = matrix(rgamma(n, z_shape), ncol = 1)
#' Z = matrix(c(rep(1, n), z_matrix[,1]), ncol = 2, byrow = FALSE)
#'
#' exp_xb = exp(X %*% true_beta)
#' pi_result = exp_xb[,1] / (exp_xb[,1] + 1)
#' pi_matrix = matrix(c(pi_result, 1 - pi_result), ncol = 2, byrow = FALSE)
#'
#' true_Y <- rep(NA, n)
#' for(i in 1:n){
#'     true_Y[i] = which(stats::rmultinom(1, 1, pi_matrix[i,]) == 1)
#' }
#'
#' exp_zg = exp(Z %*% true_gamma)
#' pistar_denominator = matrix(c(1 + exp_zg[,1], 1 + exp_zg[,2]), ncol = 2, byrow = FALSE)
#' pistar_result = exp_zg / pistar_denominator
#'
#' pistar_matrix = matrix(c(pistar_result[,1], 1 - pistar_result[,1],
#'                          pistar_result[,2], 1 - pistar_result[,2]),
#'                        ncol = 2, byrow = FALSE)
#'
#' obs_Y <- rep(NA, n)
#' for(i in 1:n){
#'     true_j = true_Y[i]
#'     obs_Y[i] = which(rmultinom(1, 1,
#'                      pistar_matrix[c(i, n + i),
#'                                      true_j]) == 1)
#'  }
#'
#' Ystar <- obs_Y
#'
#' starting_values <- rep(1,6)
#' beta_start <- matrix(starting_values[1:2], ncol = 1)
#' gamma_start <- matrix(starting_values[3:6], ncol = 2, nrow = 2, byrow = FALSE)
#'
#' EM_results <- COMBO_EM(Ystar, x_matrix = x_matrix, z_matrix = z_matrix,
#'                        beta_start = beta_start, gamma_start = gamma_start)
#'
#' EM_results
COMBO_EM <- function(Ystar,
                     x_matrix, z_matrix,
                     beta_start, gamma_start,
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

  n_cat = 2
  sample_size = length(Ystar)

  if (nrow(z_matrix) != sample_size)
    stop("The number of rows of 'z_matrix' must match the length of 'Ystar'.")
  if (!is.null(x_matrix) && nrow(x_matrix) != sample_size)
    stop("The number of rows of 'x_matrix' must match the length of 'Ystar'.")

  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)

  obs_Y_reps = matrix(rep(Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_Y_matrix = 1 * (obs_Y_reps == category_matrix)

  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)

  results = turboEM::turboem(par = c(c(beta_start), c(gamma_start)),
                             fixptfn = em_function, objfn = loglik,
                             method = c(em_method),
                             obs_Y_matrix = obs_Y_matrix,
                             X = X, Z = Z,
                             sample_size = sample_size, n_cat = n_cat,
                             control.run = control_settings)

  Ystar01 = ifelse(Ystar == 1, 1, ifelse(Ystar == 2, 0, NA))
  log_reg = stats::glm(Ystar01 ~ . + 0, as.data.frame(X),
                       family = "binomial"(link = "logit"))

  SAMBA_start <- c(beta_start, c(gamma_start)[1:(1 + ncol(z_matrix))])
  SAMBA_i <- SAMBA::obsloglikEM(Ystar01, Z = x_matrix,
                                X = z_matrix, start = SAMBA_start,
                                tol = tolerance,
                                maxit = max_em_iterations)

  perfect_sens_start <- c(beta_start, c(gamma_start)[(2 + ncol(z_matrix)):length(c(gamma_start))])
  perfect_sens_i <- perfect_sensitivity_EM(Ystar01, Z = x_matrix,
                                           X = z_matrix, start = perfect_sens_start,
                                           tolerance = tolerance,
                                           max_em_iterations = max_em_iterations)

  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
  estimates_i <- if ((pistar_11 > .50 | pistar_22 > .50) |
                     (is.na(pistar_11) & is.na(pistar_22))) {
    # If turboem cannot estimate the parameters they will be NA.
    turboEM::pars(results)
  } else {
    gamma_index = (ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
    c(-1*turboEM::pars(results)[1:ncol(X)], turboEM::pars(results)[gamma_flip_index])
  }

  sigma_EM = tryCatch(solve(turboEM::hessian(results)[[1]]), silent = TRUE,
                      error = function(e) NA)
  SE_EM = tryCatch(sqrt(diag(Matrix::nearPD(sigma_EM)$mat)),
                   silent = TRUE,
                   error = function(e) rep(NA, ncol(X) + (n_cat * ncol(Z))))

  beta_param_names <- paste0(rep("beta", ncol(X)), 1:ncol(X))
  gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  SAMBA_beta_param_names <- paste0(rep("SAMBA_beta", ncol(X)), 1:ncol(X))
  SAMBA_gamma_param_names <- paste0(rep("SAMBA_gamma", (ncol(Z))),
                              rep(1:ncol(Z), 1),
                              rep(1, ncol(Z)))
  PSens_beta_param_names <- paste0(rep("PSens_beta", ncol(X)), 1:ncol(X))
  PSens_gamma_param_names <- paste0(rep("PSens_gamma", (ncol(Z))),
                                    rep(1:ncol(Z), 1),
                                    rep(2, ncol(Z)))
  naive_beta_param_names <- paste0("naive_", beta_param_names)
  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names,
                                        SAMBA_beta_param_names,
                                        SAMBA_gamma_param_names,
                                        PSens_beta_param_names,
                                        PSens_gamma_param_names,
                                        naive_beta_param_names),
                          Estimates = c(estimates_i,
                                        unname(SAMBA_i$param),
                                        unname(perfect_sens_i$param),
                                        unname(log_reg$coefficients)),
                          SE = c(SE_EM, sqrt(diag(SAMBA_i$variance)),
                                 sqrt(diag(perfect_sens_i$variance)),
                                 unname(summary(log_reg)$coefficients[,2])),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names))),
                                          rep(NA, length(c(SAMBA_beta_param_names,
                                                           SAMBA_gamma_param_names))),
                                          rep(NA, length(c(PSens_beta_param_names,
                                                           PSens_gamma_param_names))),
                                          rep(log_reg$converged, ncol(X))))

  return(estimates)
}
